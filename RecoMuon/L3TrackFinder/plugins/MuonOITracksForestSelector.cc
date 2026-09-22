/** \class MuonOITracksForestSelector
 *
 *  \brief XGBoost forest (compact binary) selector for OI muon tracks.
 *
 *  Identical 22-feature extraction as the retired DNN selector, shared through
 *  RecoMuon/L3TrackFinder/interface/OITrackSelectorFeatures.h (named struct
 *  muonhp::OITrackFeatures; toArray() fixes the canonical training order).
 *  Replaces ONNX Runtime DNN inference with a serial tree traversal of the
 *  compact gradient-boosted-tree binary (.bin) exported by the OI XGBoost
 *  trainers: smaller model files and faster at low track multiplicity.
 *
 *  Handles both the pixel-chain and the seeds-chain (general) OI selectors:
 *  same extraction, only the .bin model and thresholds differ (set per cfi).
 *
 *  The 22-feature production set = the original 26-feature layout minus
 *  ptErr, chi2, sigmaPtOverPt and relUncertaintyProduct (pruned in the
 *  round-2 training campaign).
 *
 *  The compact .bin format (identical to MuonIOTracksForestSelector):
 *    int32 nNodes, int32 nTrees, float baseLogit,
 *    int8  feat[nNodes]   (-1 = leaf),
 *    float val[nNodes]    (threshold / leaf value),
 *    int32 left[nNodes], right[nNodes],
 *    int32 roots[nTrees]
 *
 *  Traversal: for each tree, walk from root; at internal nodes go left if
 *  x[feat[node]] < val[node], else right. At leaf, add val[node] to margin.
 *  Score = sigmoid(margin). Selection: score >= threshold (optionally
 *  pT-binned via ptBinEdges/decisionThresholds, thresholds[i] applying to
 *  pT in [edges[i], edges[i+1]) with the last bin open-ended).
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoMuon/L3TrackFinder/interface/OITrackSelectorFeatures.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

// ---------------------------------------------------------------------------
// GlobalCache: compact forest loaded ONCE per process from the .bin file.
// Binary format: int32 nNodes, int32 nTrees, float baseLogit, then
//   int8 feat[nNodes] (-1=leaf), float val[nNodes],
//   int32 left[nNodes], right[nNodes], int32 roots[nTrees].
// (Identical layout to MuonIOTracksForestSelector's ForestCache.)
// ---------------------------------------------------------------------------
struct OIForestCache {
  int nNodes = 0;
  int nTrees = 0;
  float baseLogit = 0.0f;
  std::vector<int8_t> feat;
  std::vector<float> val;
  std::vector<int32_t> left;
  std::vector<int32_t> right;
  std::vector<int32_t> roots;

  static std::unique_ptr<OIForestCache> load(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open())
      throw cms::Exception("MuonOITracksForestSelector") << "Cannot open compact forest binary: " << path;
    int32_t nNodes = 0, nTrees = 0;
    float baseLogit = 0.0f;
    in.read(reinterpret_cast<char*>(&nNodes), 4);
    in.read(reinterpret_cast<char*>(&nTrees), 4);
    in.read(reinterpret_cast<char*>(&baseLogit), 4);
    if (!in)
      throw cms::Exception("MuonOITracksForestSelector") << "compact tree binary header truncated: " << path;

    auto cache = std::make_unique<OIForestCache>();
    cache->nNodes = nNodes;
    cache->nTrees = nTrees;
    cache->baseLogit = baseLogit;
    cache->feat.resize(nNodes);
    cache->val.resize(nNodes);
    cache->left.resize(nNodes);
    cache->right.resize(nNodes);
    cache->roots.resize(nTrees);
    in.read(reinterpret_cast<char*>(cache->feat.data()), nNodes);
    in.read(reinterpret_cast<char*>(cache->val.data()), 4LL * nNodes);
    in.read(reinterpret_cast<char*>(cache->left.data()), 4LL * nNodes);
    in.read(reinterpret_cast<char*>(cache->right.data()), 4LL * nNodes);
    in.read(reinterpret_cast<char*>(cache->roots.data()), 4LL * nTrees);
    if (!in)
      throw cms::Exception("MuonOITracksForestSelector") << "compact tree binary body truncated: " << path;

    edm::LogInfo("MuonOITracksForestSelector")
        << "Loaded compact forest: nNodes=" << nNodes << " nTrees=" << nTrees << " baseLogit=" << baseLogit
        << " from " << path;
    return cache;
  }
};

class MuonOITracksForestSelector : public edm::stream::EDProducer<edm::GlobalCache<OIForestCache>> {
public:
  explicit MuonOITracksForestSelector(const edm::ParameterSet&, const OIForestCache*);
  ~MuonOITracksForestSelector() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);
  static std::unique_ptr<OIForestCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const OIForestCache*) {}

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  // Working point for one track: global threshold, or the threshold of the
  // pT bin the track falls into.
  float thresholdForPt(double pt) const {
    if (ptBinEdges_.empty())
      return decisionThreshold_;
    size_t b = 0;
    while (b + 1 < ptBinEdges_.size() && pt >= ptBinEdges_[b + 1])
      ++b;
    return static_cast<float>(decisionThresholds_[b]);
  }

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> standaloneMuonsToken_;

  // Model parameters
  const float decisionThreshold_;
  // Optional pT-binned working points (same convention as
  // MuonIOTracksForestSelector): bin i applies decisionThresholds_[i] to
  // tracks with pT in [ptBinEdges_[i], ptBinEdges_[i+1]), last bin
  // open-ended; both vectors empty or equal-sized.
  const std::vector<double> ptBinEdges_;
  const std::vector<double> decisionThresholds_;
  const bool useStandaloneMuonFeatures_;
  const int nFeatures_;
  const bool dumpFeatures_;
  unsigned int eventCounter_;


};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

MuonOITracksForestSelector::MuonOITracksForestSelector(const edm::ParameterSet& iConfig, const OIForestCache* cache)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      standaloneMuonsToken_(
          consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("standaloneMuons"))),
      decisionThreshold_(iConfig.getParameter<double>("decisionThreshold")),
      ptBinEdges_(iConfig.getParameter<std::vector<double>>("ptBinEdges")),
      decisionThresholds_(iConfig.getParameter<std::vector<double>>("decisionThresholds")),
      useStandaloneMuonFeatures_(iConfig.getParameter<bool>("useStandaloneMuonFeatures")),
      nFeatures_(iConfig.getParameter<int>("nFeatures")),
      dumpFeatures_(iConfig.getUntrackedParameter<bool>("dumpFeatures")),
      eventCounter_(0) {
  if (!ptBinEdges_.empty()) {
    if (decisionThresholds_.size() != ptBinEdges_.size())
      throw cms::Exception("MuonOITracksForestSelector")
          << "pT-binned working points are inconsistent: " << decisionThresholds_.size()
          << " decisionThresholds but " << ptBinEdges_.size()
          << " ptBinEdges; they must have equal size (thresholds[i] applies to pT in [edges[i], edges[i+1]), "
             "last bin open-ended).";
    for (size_t i = 1; i < ptBinEdges_.size(); ++i)
      if (ptBinEdges_[i] <= ptBinEdges_[i - 1])
        throw cms::Exception("MuonOITracksForestSelector") << "ptBinEdges must be strictly increasing.";
  } else if (!decisionThresholds_.empty()) {
    throw cms::Exception("MuonOITracksForestSelector")
        << "decisionThresholds set without ptBinEdges; set both for pT-binned working points or neither.";
  }
  produces<reco::TrackCollection>();
  produces<std::vector<float>>("scores");
}

std::unique_ptr<OIForestCache> MuonOITracksForestSelector::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  edm::FileInPath modelPath(iConfig.getParameter<std::string>("modelPath"));
  return OIForestCache::load(modelPath.fullPath());
}


void MuonOITracksForestSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonOITracksForestSelector";

  auto selectedTracks = std::make_unique<reco::TrackCollection>();
  auto scores = std::make_unique<std::vector<float>>();

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<reco::TrackCollection> standaloneMuons;
  iEvent.getByToken(standaloneMuonsToken_, standaloneMuons);

  if (!tracks.isValid() || tracks->empty()) {
    iEvent.put(std::move(selectedTracks));
    iEvent.put(std::move(scores), "scores");
    return;
  }

  const auto* cache = globalCache();

  const unsigned int evtIdx = eventCounter_++;
  for (size_t i = 0; i < tracks->size(); ++i) {
    const auto& track = (*tracks)[i];
    const auto featureArray = muonhp::extractOITrackFeatures(track, *standaloneMuons).toArray();
    const auto& features = featureArray;
    if (static_cast<int>(features.size()) != nFeatures_) {
      throw cms::Exception("MuonOITracksForestSelector")
          << "Feature count mismatch: extracted " << features.size() << " features, expected " << nFeatures_
          << ". The deployed 22-feature ABI is fixed; check nFeatures vs the trained forest model.";
    }

    // --- Forest inference: serial tree traversal ---
    float margin = cache->baseLogit;
    for (int t = 0; t < cache->nTrees; ++t) {
      int32_t node = cache->roots[t];
      while (cache->feat[node] >= 0)
        node = (features[cache->feat[node]] < cache->val[node]) ? cache->left[node] : cache->right[node];
      margin += cache->val[node];
    }
    float prob = 1.0f / (1.0f + std::exp(-margin));
    prob = std::clamp(prob, 0.0f, 1.0f);

    scores->push_back(prob);

    // Optional dump for cross-validation
    if (dumpFeatures_) {
      std::ostringstream oss;
      oss << evtIdx << "," << i;
      oss << std::scientific << std::setprecision(9);
      for (float f : features)
        oss << "," << f;
      oss << "," << prob;
      std::cout << oss.str() << "\n";
    }

    if (prob >= thresholdForPt(track.pt())) {
      selectedTracks->push_back(*(reco::TrackRef(tracks, i)));
    }
  }

  LogTrace(metname) << " Selected " << selectedTracks->size() << " out of " << tracks->size() << " tracks";

  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(scores), "scores");
}

void MuonOITracksForestSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks", edm::InputTag("hltPhase2L3OIMuCtfWithMaterialTracks"))
      ->setComment("Input OI track collection");
  desc.add<edm::InputTag>("standaloneMuons", edm::InputTag("hltL2MuonsFromL1TkMuon", "UpdatedAtVtx"))
      ->setComment("Standalone (L2 muon vertex) track collection for matching features");
  desc.add<std::string>("modelPath", "RecoMuon/L3TrackFinder/data/OI_track_selector_forest.bin")
      ->setComment("Path to compact gradient-boosted-tree binary (.bin)");
  desc.add<double>("decisionThreshold", 0.5)
      ->setComment("Probability threshold for track selection (use F2-optimal from training); "
                   "used for every track when ptBinEdges is empty");
  desc.add<std::vector<double>>("ptBinEdges", {})
      ->setComment("Lower pT edges [GeV] of per-bin working points (thresholds[i] applies to pT in "
                   "[edges[i], edges[i+1]), last bin open-ended); from the training thresholds.json "
                   "pt_bin_edges. Empty = single-threshold mode.");
  desc.add<std::vector<double>>("decisionThresholds", {})
      ->setComment("Per-pT-bin probability thresholds (must match ptBinEdges in size and ordering); "
                   "from the training thresholds.json pt_bin_f2_thresholds.");
  desc.add<bool>("useStandaloneMuonFeatures", true)->setComment("Include standalone-muon matching features");
  desc.add<int>("nFeatures", 22)
      ->setComment("Total number of input features the forest was trained with "
                   "(always 22: the extractor emits exactly the OI production set)");
  desc.addUntracked<bool>("dumpFeatures", false)
      ->setComment("Print one CSV-like line per track to stdout, for cross-validation against build_dataset().");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonOITracksForestSelector);
