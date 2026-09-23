/** \class MuonIOTracksForestSelector
 *
 *  \brief XGBoost-forest HighPurity selector for inside-out (IO) muon tracks.
 *
 *  Serves both IO selectors of the Phase-2 muon HLT: the pixel tracks of the
 *  pixel-track chain (hltPhase2MuonPixelTracks) and the IO tracks built from
 *  LST seeds in the seeds chain (hltPhase2MuonIOTracks). Both feed the same
 *  33-feature extraction (RecoMuon/L3TrackFinder/interface/
 *  IOTrackSelectorFeatures.h); only the forest (.bin) and the working points
 *  differ (set per cfi, from the training's thresholds.json).
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

#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoMuon/L3TrackFinder/interface/CompactForest.h"
#include "RecoMuon/L3TrackFinder/interface/IOTrackSelectorFeatures.h"

#include <algorithm>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

class MuonIOTracksForestSelector : public edm::stream::EDProducer<edm::GlobalCache<muonhp::CompactForest>> {
public:
  MuonIOTracksForestSelector(const edm::ParameterSet&, const muonhp::CompactForest*);

  static void fillDescriptions(edm::ConfigurationDescriptions&);
  static std::unique_ptr<muonhp::CompactForest> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const muonhp::CompactForest*) {}

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuonsToken_;
  const muonhp::BinnedWorkingPoints workingPoints_;
  const bool dumpFeatures_;
};

MuonIOTracksForestSelector::MuonIOTracksForestSelector(const edm::ParameterSet& iConfig, const muonhp::CompactForest*)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      l1TkMuonsToken_(consumes<l1t::TrackerMuonCollection>(iConfig.getParameter<edm::InputTag>("l1TkMuons"))),
      workingPoints_(iConfig.getParameter<double>("decisionThreshold"),
                     iConfig.getParameter<std::vector<double>>("ptBinEdges"),
                     iConfig.getParameter<std::vector<double>>("decisionThresholds")),
      dumpFeatures_(iConfig.getUntrackedParameter<bool>("dumpFeatures")) {
  produces<reco::TrackCollection>();
  produces<std::vector<float>>("scores");
}

std::unique_ptr<muonhp::CompactForest> MuonIOTracksForestSelector::initializeGlobalCache(
    const edm::ParameterSet& iConfig) {
  constexpr int kFeatures = muonhp::IOTrackFeatures::kSize;
  const int nFeatures = iConfig.getParameter<int>("nFeatures");
  if (nFeatures != kFeatures)
    throw cms::Exception("Configuration") << "MuonIOTracksForestSelector: nFeatures = " << nFeatures
                                          << " but the IO extractor provides " << kFeatures << " features";
  const edm::FileInPath modelPath(iConfig.getParameter<edm::FileInPath>("modelPath"));
  auto forest = muonhp::CompactForest::load(modelPath.fullPath(), kFeatures);
  edm::LogInfo("MuonIOTracksForestSelector")
      << "Loaded " << modelPath.relativePath() << ": " << forest->nTrees() << " trees, " << forest->nNodes()
      << " nodes, baseLogit " << forest->baseLogit();
  return forest;
}

void MuonIOTracksForestSelector::produce(edm::Event& iEvent, const edm::EventSetup&) {
  const auto& tracks = iEvent.get(tracksToken_);
  const auto& l1TkMuons = iEvent.get(l1TkMuonsToken_);
  const auto* forest = globalCache();

  // Features of all tracks, then one batched forest evaluation (CompactForest
  // evaluates tree-major over the tracks: same scores, better cache reuse).
  constexpr size_t kFeatures = muonhp::IOTrackFeatures::kSize;
  std::vector<float> features(tracks.size() * kFeatures);
  for (size_t i = 0; i < tracks.size(); ++i) {
    const auto f = muonhp::extractIOTrackFeatures(tracks[i], l1TkMuons).toArray();
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
      edm::LogPrint("MuonIOTracksForestSelector") << line.str();
    }
  }

  LogTrace("MuonIOTracksForestSelector") << "selected " << selectedTracks->size() << " of " << tracks.size()
                                         << " tracks";
  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(scores), "scores");
}

void MuonIOTracksForestSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("hltPhase2MuonPixelTracks"))->setComment("Input track collection");
  desc.add<edm::InputTag>("l1TkMuons", edm::InputTag("l1tTkMuonsGmt"))
      ->setComment("L1 tracker muons (matching and stub features)");
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
  desc.add<int>("nFeatures", muonhp::IOTrackFeatures::kSize)
      ->setComment("Feature count the forest was trained with; must equal the extractor's (33)");
  desc.addUntracked<bool>("dumpFeatures", false)
      ->setComment("Print one MUONHP_FEATURES line per track (features + score) for the training cross-check");
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonIOTracksForestSelector);
