/** \class MuonOITracksForestSelector
 *
 *  \brief XGBoost-forest HighPurity selector for outside-in (OI) muon tracks.
 *
 *  Serves the OI selection of both Phase-2 muon HLT chains (pixel-track chain
 *  and seeds chain): same 22-feature extraction with standalone-muon matching
 *  (RecoMuon/L3TrackFinder/interface/OITrackSelectorFeatures.h), chain-
 *  specific forest (.bin) and working points (set per cfi, from the
 *  training's thresholds.json).
 *
 *  Inference: batched traversal of the compact forest over the event's
 *  tracks (interface/CompactForest.h), score = sigmoid(margin). A track is kept if its score is
 *  >= the working point of its pT bin (ptBinEdges/decisionThresholds; single
 *  decisionThreshold if no bins are configured). Outputs the selected track
 *  collection and the scores of all input tracks (input order).
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoMuon/L3TrackFinder/interface/CompactForest.h"
#include "RecoMuon/L3TrackFinder/interface/OITrackSelectorFeatures.h"

#include <algorithm>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

class MuonOITracksForestSelector : public edm::stream::EDProducer<edm::GlobalCache<muonhp::CompactForest>> {
public:
  MuonOITracksForestSelector(const edm::ParameterSet&, const muonhp::CompactForest*);

  static void fillDescriptions(edm::ConfigurationDescriptions&);
  static std::unique_ptr<muonhp::CompactForest> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const muonhp::CompactForest*) {}

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<reco::TrackCollection> standaloneMuonsToken_;
  const muonhp::BinnedWorkingPoints workingPoints_;
  const bool dumpFeatures_;
};

MuonOITracksForestSelector::MuonOITracksForestSelector(const edm::ParameterSet& iConfig, const muonhp::CompactForest*)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      standaloneMuonsToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("standaloneMuons"))),
      workingPoints_(iConfig.getParameter<double>("decisionThreshold"),
                     iConfig.getParameter<std::vector<double>>("ptBinEdges"),
                     iConfig.getParameter<std::vector<double>>("decisionThresholds")),
      dumpFeatures_(iConfig.getUntrackedParameter<bool>("dumpFeatures")) {
  produces<reco::TrackCollection>();
  produces<std::vector<float>>("scores");
}

std::unique_ptr<muonhp::CompactForest> MuonOITracksForestSelector::initializeGlobalCache(
    const edm::ParameterSet& iConfig) {
  constexpr int kFeatures = muonhp::OITrackFeatures::kSize;
  const int nFeatures = iConfig.getParameter<int>("nFeatures");
  if (nFeatures != kFeatures)
    throw cms::Exception("Configuration") << "MuonOITracksForestSelector: nFeatures = " << nFeatures
                                          << " but the OI extractor provides " << kFeatures << " features";
  const edm::FileInPath modelPath(iConfig.getParameter<edm::FileInPath>("modelPath"));
  auto forest = muonhp::CompactForest::load(modelPath.fullPath(), kFeatures);
  edm::LogInfo("MuonOITracksForestSelector")
      << "Loaded " << modelPath.relativePath() << ": " << forest->nTrees() << " trees, " << forest->nNodes()
      << " nodes, baseLogit " << forest->baseLogit();
  return forest;
}

void MuonOITracksForestSelector::produce(edm::Event& iEvent, const edm::EventSetup&) {
  const auto& tracks = iEvent.get(tracksToken_);
  const auto& standaloneMuons = iEvent.get(standaloneMuonsToken_);
  const auto* forest = globalCache();

  // Features of all tracks, then one batched forest evaluation (CompactForest
  // evaluates tree-major over the tracks: same scores, better cache reuse).
  constexpr size_t kFeatures = muonhp::OITrackFeatures::kSize;
  std::vector<float> features(tracks.size() * kFeatures);
  for (size_t i = 0; i < tracks.size(); ++i) {
    const auto f = muonhp::extractOITrackFeatures(tracks[i], standaloneMuons).toArray();
    std::copy(f.begin(), f.end(), features.begin() + i * kFeatures);
  }
  auto scores = std::make_unique<std::vector<float>>(tracks.size());
  forest->evaluate(features.data(), tracks.size(), kFeatures, scores->data());

  auto selectedTracks = std::make_unique<reco::TrackCollection>();
  for (size_t i = 0; i < tracks.size(); ++i) {
    const float score = (*scores)[i];
    if (score >= workingPoints_.threshold(tracks[i].pt()))
      selectedTracks->push_back(tracks[i]);

    if (dumpFeatures_) {
      // One line per track, for the Python/C++ feature and score cross-check
      // (muonHighPurityTrackSelection/production/features_validation).
      std::ostringstream line;
      line << "MUONHP_FEATURES," << moduleDescription().moduleLabel() << ',' << iEvent.id().run() << ','
           << iEvent.id().luminosityBlock() << ',' << iEvent.id().event() << ',' << i << std::scientific
           << std::setprecision(9);
      for (size_t k = 0; k < kFeatures; ++k)
        line << ',' << features[i * kFeatures + k];
      line << ',' << score;
      edm::LogPrint("MuonOITracksForestSelector") << line.str();
    }
  }

  LogTrace("MuonOITracksForestSelector") << "selected " << selectedTracks->size() << " of " << tracks.size()
                                         << " tracks";
  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(scores), "scores");
}

void MuonOITracksForestSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("hltPhase2L3OIMuCtfWithMaterialTracks"))
      ->setComment("Input OI track collection");
  desc.add<edm::InputTag>("standaloneMuons", edm::InputTag("hltL2MuonsFromL1TkMuon", "UpdatedAtVtx"))
      ->setComment("Standalone muons updated at the vertex (matching features)");
  desc.add<edm::FileInPath>("modelPath")
      ->setComment(
          "Compact forest binary (.bin) exported by the training pipeline (RecoMuon/L3TrackFinder/data/IO|OI)");
  desc.add<double>("decisionThreshold", 0.5)
      ->setComment("Score threshold used for every track when ptBinEdges is empty (training global F2 point)");
  desc.add<std::vector<double>>("ptBinEdges", {})
      ->setComment(
          "Lower pT edges [GeV] of the per-bin working points (last bin open-ended); "
          "training thresholds.json pt_bin_edges. Empty = single threshold.");
  desc.add<std::vector<double>>("decisionThresholds", {})
      ->setComment(
          "Per-pT-bin score thresholds, one per ptBinEdges entry; training thresholds.json "
          "pt_bin_f2_thresholds");
  desc.add<int>("nFeatures", muonhp::OITrackFeatures::kSize)
      ->setComment("Feature count the forest was trained with; must equal the extractor's (22)");
  desc.addUntracked<bool>("dumpFeatures", false)
      ->setComment("Print one MUONHP_FEATURES line per track (features + score) for the training cross-check");
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonOITracksForestSelector);
