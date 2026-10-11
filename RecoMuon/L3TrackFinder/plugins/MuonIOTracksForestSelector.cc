/** \class MuonIOTracksForestSelector
 *
 *  \brief XGBoost-forest HighPurity selector for inside-out (IO) muon tracks.
 *
 *  Serves both IO selectors of the Phase-2 muon HLT: the pixel tracks of the
 *  pixel-track chain (hltPhase2MuonPixelTracks) and the IO tracks built from
 *  LST seeds in the seeds chain (hltPhase2MuonIOTracks). Both feed the same
 *  feature extraction (RecoMuon/L3TrackFinder/interface/
 *  IOTrackSelectorFeatures.h); only the forests (.bin) and the working points
 *  differ (set per cfi, from the training's thresholds.json).
 *
 *  Features: track parameters, the stub count of the matched L1 tracker muon
 *  and the raw matching quantities of the candidates of every muon station
 *  (DT, CSC and ME0 segments, GE1/1, GE2/1 and RPC rec-hits;
 *  interface/MuonChamberMatch.h).
 *
 *  Inference: batched traversal of the compact forest over the event's tracks
 *  (interface/CompactForest.h), score = sigmoid(margin). A track is
 *  kept if its score is >= the working point of its pT bin
 *  (ptBinEdges/decisionThresholds; single decisionThreshold if no bins are
 *  configured). Outputs the selected track collection and the scores of all
 *  input tracks (input order).
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "RecoMuon/L3TrackFinder/interface/CompactForest.h"
#include "RecoMuon/L3TrackFinder/interface/IOTrackSelectorFeatures.h"
#include "RecoMuon/L3TrackFinder/interface/MuonChamberMatch.h"

#include <algorithm>
#include <cmath>
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
  void fillChamberPoints(const edm::Event&, const edm::EventSetup&);

  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuonsToken_;
  const edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken_;
  const edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken_;
  const edm::EDGetTokenT<GEMSegmentCollection> gemSegmentsToken_;
  const edm::EDGetTokenT<GEMRecHitCollection> gemRecHitsToken_;
  const edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsToken_;
  const edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeometryToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryToken_;
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeometryToken_;
  const edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeometryToken_;
  const muonhp::BinnedWorkingPoints workingPoints_;
  const bool dumpFeatures_;

  // per-event buffer
  std::vector<muonhp::ChamberPoint> points_;
};

MuonIOTracksForestSelector::MuonIOTracksForestSelector(const edm::ParameterSet& iConfig, const muonhp::CompactForest*)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      l1TkMuonsToken_(consumes<l1t::TrackerMuonCollection>(iConfig.getParameter<edm::InputTag>("l1TkMuons"))),
      dtSegmentsToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dtSegments"))),
      cscSegmentsToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"))),
      gemSegmentsToken_(consumes<GEMSegmentCollection>(iConfig.getParameter<edm::InputTag>("gemSegments"))),
      gemRecHitsToken_(consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"))),
      rpcRecHitsToken_(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("rpcRecHits"))),
      dtGeometryToken_(esConsumes()),
      cscGeometryToken_(esConsumes()),
      gemGeometryToken_(esConsumes()),
      rpcGeometryToken_(esConsumes()),
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
  edm::LogInfo("MuonIOTracksForestSelector") << "Loaded " << modelPath.relativePath() << ": " << forest->nTrees()
                                             << " trees, " << forest->nNodes() << " nodes, baseLogit "
                                             << forest->baseLogit();
  return forest;
}

// The event's segments and rec-hits in global coordinates, in the order of the training code: DT, CSC, ME0, GEM,
// RPC barrel, RPC endcap, each in the order of its collection.
void MuonIOTracksForestSelector::fillChamberPoints(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace muonhp;
  points_.clear();

  // loops over a detector collection, as the n-tuple table producers do
  auto forEach = [](const auto& collection, const auto& geometry, auto&& f) {
    for (auto idIt = collection.id_begin(); idIt != collection.id_end(); ++idIt) {
      const auto range = collection.get(*idIt);
      const GeomDet* det = geometry.idToDet(*idIt);
      for (auto it = range.first; it != range.second; ++it)
        f(*idIt, *it, *det);
    }
  };

  auto segment = [&](const auto& s, const GeomDet& det, uint8_t station) {
    const GlobalPoint p = det.toGlobal(s.localPosition());
    const GlobalVector d = det.toGlobal(s.localDirection());
    points_.push_back({chamberEta(p.x(), p.y(), p.z()), float(p.phi().value()), float(d.phi().value()), station, true});
  };

  const auto& dtGeometry = iSetup.getData(dtGeometryToken_);
  forEach(iEvent.get(dtSegmentsToken_),
          dtGeometry,
          [&](const DTChamberId& id, const DTRecSegment4D& s, const GeomDet& det) {
            if (s.hasPhi())
              segment(s, det, id.station());  // kSt1..kSt4 = 1..4
          });

  const auto& cscGeometry = iSetup.getData(cscGeometryToken_);
  forEach(iEvent.get(cscSegmentsToken_),
          cscGeometry,
          [&](const CSCDetId& id, const CSCSegment& s, const GeomDet& det) { segment(s, det, id.station()); });

  const auto& gemGeometry = iSetup.getData(gemGeometryToken_);
  forEach(iEvent.get(gemSegmentsToken_),
          gemGeometry,
          [&](const GEMDetId& id, const GEMSegment& s, const GeomDet& det) {
            if (id.station() == 0)  // ME0
              segment(s, det, kME0);
          });

  forEach(iEvent.get(gemRecHitsToken_), gemGeometry, [&](const GEMDetId& id, const GEMRecHit& h, const GeomDet& det) {
    // GE1/1 and GE2/1 (the ME0 hits are not used: its segments are)
    if ((id.station() != 1 && id.station() != 2) || std::abs(h.BunchX()) > kChamberMaxHitBx)
      return;
    const GlobalPoint p = det.toGlobal(h.localPosition());
    // eta from z and r as stored (float) in the training n-tuple
    const double z = float(p.z()), r = chamberPerp(p.x(), p.y());
    points_.push_back(
        {float(std::asinh(z / r)), float(p.phi().value()), 0.f, uint8_t(id.station() == 1 ? kGE11 : kGE21), false});
  });

  const auto& rpcGeometry = iSetup.getData(rpcGeometryToken_);
  const auto& rpcRecHits = iEvent.get(rpcRecHitsToken_);
  for (const bool barrel : {true, false}) {
    forEach(rpcRecHits, rpcGeometry, [&](const RPCDetId& id, const RPCRecHit& h, const GeomDet& det) {
      if ((id.region() == 0) != barrel || std::abs(h.BunchX()) > kChamberMaxHitBx)
        return;
      const GlobalPoint p = det.toGlobal(h.localPosition());
      points_.push_back(
          {chamberEta(p.x(), p.y(), p.z()), float(p.phi().value()), 0.f, uint8_t(kRPC1 + id.station() - 1), false});
    });
  }
}

void MuonIOTracksForestSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace muonhp;
  const auto& tracks = iEvent.get(tracksToken_);
  const auto& l1TkMuons = iEvent.get(l1TkMuonsToken_);
  const size_t nTracks = tracks.size();

  auto scores = std::make_unique<std::vector<float>>(nTracks);
  auto selectedTracks = std::make_unique<reco::TrackCollection>();
  if (nTracks == 0) {
    iEvent.put(std::move(selectedTracks));
    iEvent.put(std::move(scores), "scores");
    return;
  }

  // Features of all tracks, then one batched forest evaluation (CompactForest
  // evaluates tree-major over the tracks: same scores, better cache reuse).
  fillChamberPoints(iEvent, iSetup);
  constexpr size_t kFeatures = IOTrackFeatures::kSize;
  std::vector<float> features(nTracks * kFeatures);
  for (size_t i = 0; i < nTracks; ++i) {
    const auto f = extractIOTrackFeatures(tracks[i], l1TkMuons, points_).toArray();
    std::copy(f.begin(), f.end(), features.begin() + i * kFeatures);
  }
  globalCache()->evaluate(features.data(), nTracks, kFeatures, scores->data());

  for (size_t i = 0; i < nTracks; ++i) {
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

  LogTrace("MuonIOTracksForestSelector") << "selected " << selectedTracks->size() << " of " << nTracks << " tracks ("
                                         << points_.size() << " segments and hits)";
  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(scores), "scores");
}

void MuonIOTracksForestSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("hltPhase2MuonPixelTracks"))->setComment("Input track collection");
  desc.add<edm::InputTag>("l1TkMuons", edm::InputTag("l1tTkMuonsGmt"))
      ->setComment("L1 tracker muons (stub count of the matched one)");
  desc.add<edm::InputTag>("dtSegments", edm::InputTag("hltDt4DSegments"))->setComment("DT 4D segments");
  desc.add<edm::InputTag>("cscSegments", edm::InputTag("hltCscSegments"))->setComment("CSC segments");
  desc.add<edm::InputTag>("gemSegments", edm::InputTag("hltGemSegments"))
      ->setComment("GEM segments (the ME0 ones are used)");
  desc.add<edm::InputTag>("gemRecHits", edm::InputTag("hltGemRecHits"))->setComment("GEM rec-hits (GE1/1, GE2/1)");
  desc.add<edm::InputTag>("rpcRecHits", edm::InputTag("hltRpcRecHits"))->setComment("RPC rec-hits");
  desc.add<edm::FileInPath>("modelPath")
      ->setComment("Compact forest (.bin) exported by the training pipeline (RecoMuon/L3TrackFinder/data/IO)");
  desc.add<double>("decisionThreshold", 0.5)
      ->setComment("Score threshold used for every track when ptBinEdges is empty");
  desc.add<std::vector<double>>("ptBinEdges", {})
      ->setComment(
          "Lower pT edges [GeV] of the per-bin working points (last bin open-ended); "
          "training thresholds.json pt_bin_edges. Empty = single threshold.");
  desc.add<std::vector<double>>("decisionThresholds", {})
      ->setComment(
          "Per-pT-bin score thresholds, one per ptBinEdges entry; training thresholds.json "
          "pt_bin_f2_thresholds");
  desc.add<int>("nFeatures", muonhp::IOTrackFeatures::kSize)
      ->setComment("Feature count the forest was trained with; must equal the extractor's");
  desc.addUntracked<bool>("dumpFeatures", false)
      ->setComment("Print one MUONHP_FEATURES line per track (features + score) for the training cross-check");
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonIOTracksForestSelector);
