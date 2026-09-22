/** \class MuonIOTracksForestSelector
 *
 *  \brief XGBoost forest (compact binary) muon track HighPurity selector.
 *
 *  Identical 33-feature extraction as MuonIOTracksDNNSelector, but replaces
 *  ONNX Runtime DNN inference with a serial tree traversal of the compact
 *  gradient-boosted-tree binary (.bin). This is 2.5x smaller than the ONNX
 *  model and ~2x faster at the typical 3-12 muon tracks per event.
 *
 *  Handles both the pixel-track and the IO-track (seeds) selector: both feed
 *  reco::Track + L1TkMu matching features through the same 33-feature
 *  extraction; only the .bin model and threshold differ (set per cfi).
 *
 *  The compact .bin format (identical to PixelTrackForestHighPuritySelector):
 *    int32 nNodes, int32 nTrees, float baseLogit,
 *    int8  feat[nNodes]   (-1 = leaf),
 *    float val[nNodes]    (threshold / leaf value),
 *    int32 left[nNodes], right[nNodes],
 *    int32 roots[nTrees]
 *
 *  Traversal: for each tree, walk from root; at internal nodes go left if
 *  x[feat[node]] < val[node], else right. At leaf, add val[node] to margin.
 *  Score = sigmoid(margin). Selection: score >= threshold.
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
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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
// ---------------------------------------------------------------------------
struct ForestCache {
  int nNodes = 0;
  int nTrees = 0;
  float baseLogit = 0.0f;
  std::vector<int8_t> feat;
  std::vector<float> val;
  std::vector<int32_t> left;
  std::vector<int32_t> right;
  std::vector<int32_t> roots;

  static std::unique_ptr<ForestCache> load(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
      throw cms::Exception("MuonIOTracksForestSelector") << "cannot open compact tree binary: " << path;

    int32_t nNodes = 0, nTrees = 0;
    float baseLogit = 0.0f;
    in.read(reinterpret_cast<char*>(&nNodes), 4);
    in.read(reinterpret_cast<char*>(&nTrees), 4);
    in.read(reinterpret_cast<char*>(&baseLogit), 4);
    if (!in)
      throw cms::Exception("MuonIOTracksForestSelector") << "compact tree binary header truncated: " << path;

    auto cache = std::make_unique<ForestCache>();
    cache->nNodes = nNodes;
    cache->nTrees = nTrees;
    cache->baseLogit = baseLogit;
    cache->feat.resize(nNodes);
    cache->val.resize(nNodes);
    cache->left.resize(nNodes);
    cache->right.resize(nNodes);
    cache->roots.resize(nTrees);

    in.read(reinterpret_cast<char*>(cache->feat.data()), nNodes);
    in.read(reinterpret_cast<char*>(cache->val.data()), nNodes * 4);
    in.read(reinterpret_cast<char*>(cache->left.data()), nNodes * 4);
    in.read(reinterpret_cast<char*>(cache->right.data()), nNodes * 4);
    in.read(reinterpret_cast<char*>(cache->roots.data()), nTrees * 4);
    if (!in)
      throw cms::Exception("MuonIOTracksForestSelector") << "compact tree binary truncated/corrupt: " << path;

    edm::LogInfo("MuonIOTracksForestSelector") << "Loaded compact forest: nNodes=" << nNodes << " nTrees=" << nTrees
                                               << " baseLogit=" << baseLogit << " from " << path;
    return cache;
  }
};

class MuonIOTracksForestSelector : public edm::stream::EDProducer<edm::GlobalCache<ForestCache>> {
public:
  explicit MuonIOTracksForestSelector(const edm::ParameterSet&, const ForestCache*);
  ~MuonIOTracksForestSelector() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);
  static std::unique_ptr<ForestCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ForestCache*) {}

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  std::vector<float> extractFeatures(const reco::Track& track, const l1t::TrackerMuonCollection& l1TkMuons) const;

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuonsToken_;

  // Model parameters
  const float decisionThreshold_;
  const bool useL1TkMuFeatures_;
  const bool useStubFeatures_;
  const int nFeatures_;
  const bool dumpFeatures_;
  unsigned int eventCounter_;

  // L1 Matching parameters
  static constexpr float kMatchDR2Cut = 0.09f;    // 0.3^2
  static constexpr float kMatchChi2PtCut = 9.0f;  // 3 sigma

  // Loose matching window for nCompatible feature
  static constexpr float kLooseDR2Cut = 0.25f;  // 0.5^2
  static constexpr float kLooseChi2PtCut = 25.0f;

  // Sentinel for second-best dR2
  static constexpr float kSentinel = 999.0f;

  // Regularization constants - chosen to exactly mirror modelV6.py
  static constexpr float kEpsilon = 1e-6f;     // generic log/division floor
  static constexpr float kChi2PtEps = 1e-12f;  // matches Python (t_ptErr**2 + 1e-12)
  static constexpr float kDPtNormEps = 1e-9f;  // matches Python (l1_pt + 1e-9)

  // Imputation values for non-matched L1 features
  static constexpr float kImputeDR2 = 0.1f;
  static constexpr float kImputeDPtNorm = 1.0f;
  static constexpr float kImputeChi2Pt = 10.0f;
  static constexpr float kImputeMatchScore = 0.2f;
  static constexpr float kImputeSecondDR2 = 1.0f;
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

MuonIOTracksForestSelector::MuonIOTracksForestSelector(const edm::ParameterSet& iConfig, const ForestCache* cache)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      l1TkMuonsToken_(consumes<l1t::TrackerMuonCollection>(iConfig.getParameter<edm::InputTag>("l1TkMuons"))),
      decisionThreshold_(iConfig.getParameter<double>("decisionThreshold")),
      useL1TkMuFeatures_(iConfig.getParameter<bool>("useL1TkMuFeatures")),
      useStubFeatures_(iConfig.getParameter<bool>("useStubFeatures")),
      nFeatures_(iConfig.getParameter<int>("nFeatures")),
      dumpFeatures_(iConfig.getUntrackedParameter<bool>("dumpFeatures")),
      eventCounter_(0) {
  produces<reco::TrackCollection>();
  produces<std::vector<float>>("scores");
}

std::unique_ptr<ForestCache> MuonIOTracksForestSelector::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  edm::FileInPath modelPath(iConfig.getParameter<std::string>("modelPath"));
  return ForestCache::load(modelPath.fullPath());
}

std::vector<float> MuonIOTracksForestSelector::extractFeatures(const reco::Track& track,
                                                               const l1t::TrackerMuonCollection& l1TkMuons) const {
  // -----------------------------------------------------------------------
  // Pruned 33-feature model (v7).
  // Identical feature extraction to MuonIOTracksDNNSelector -- the only
  // change vs the DNN selector is the inference backend (forest .bin vs ONNX).
  //
  // Drops 11 redundant features from the original 44-feature set:
  //   chi2, normalizedChi2, dszErr, dxyErr, nLostHits,
  //   impactSignificance, chi2PerHit, hitEfficiency,
  //   eta (redundant with absEta), ptErr (derivable),
  //   relUncertaintyProduct (corr 0.99999 with sigmaPtOverPt)
  // -----------------------------------------------------------------------

  std::vector<float> features;
  features.reserve(nFeatures_);

  // Raw track quantities
  const float p = track.p();
  const float pt = track.pt();
  const float ptErr = track.ptError();
  const float eta = track.eta();
  const float etaErr = track.etaError();
  const float phi = track.phi();
  const float phiErr = track.phiError();

  const float dxy = track.dxy();
  const float dxyErr = track.dxyError();
  const float dz = track.dz();
  const float dzErr = track.dzError();
  const float qoverpErr = track.qoverpError();
  const float lambdaErr = track.lambdaError();

  const int nPixelHits = track.hitPattern().numberOfValidPixelHits();
  const int nTrkLays = track.hitPattern().trackerLayersWithMeasurement();
  const int nFoundHits = track.numberOfValidHits();

  // Features 0-6: Log features (dropped: ptErr, chi2, normalizedChi2, dszErr, dxyErr)
  features.push_back(std::log10(std::abs(p) + kEpsilon));          // 0: p (log)
  features.push_back(std::log10(std::abs(pt) + kEpsilon));         // 1: pt (log)
  features.push_back(std::log10(std::abs(etaErr) + kEpsilon));     // 2: etaErr (log)
  features.push_back(std::log10(std::abs(phiErr) + kEpsilon));     // 3: phiErr (log)
  features.push_back(std::log10(std::abs(dzErr) + kEpsilon));      // 4: dzErr (log)
  features.push_back(std::log10(std::abs(qoverpErr) + kEpsilon));  // 5: qoverpErr (log)
  features.push_back(std::log10(std::abs(lambdaErr) + kEpsilon));  // 6: lambdaErr (log)

  // Features 7-9: Plain features (dropped: eta, nLostHits)
  features.push_back(static_cast<float>(nPixelHits));  // 7: nPixelHits
  features.push_back(static_cast<float>(nTrkLays));    // 8: nTrkLays
  features.push_back(static_cast<float>(nFoundHits));  // 9: nFoundHits

  // Features 10-16: Derived features (dropped: impactSignificance, chi2PerHit,
  //   hitEfficiency, relUncertaintyProduct)
  // 10: Impact Parameter 3D (log)
  const float ip3d = dxy * dxy + dz * dz;
  features.push_back(std::log10(ip3d + kEpsilon));

  // 11: SigmaPt / Pt (log)
  const float sigmaPtOverPt = ptErr / std::max(pt, kEpsilon);
  features.push_back(std::log10(sigmaPtOverPt + kEpsilon));

  // 12: Separated 2D impact parameter significance (log)
  const float sip2D = std::abs(dxy) / std::max(dxyErr, kEpsilon);
  features.push_back(std::log10(sip2D + kEpsilon));

  // 13: Longitudinal impact parameter significance (log)
  const float sipZ = std::abs(dz) / std::max(dzErr, kEpsilon);
  features.push_back(std::log10(sipZ + kEpsilon));

  // 14: |dxy| / pT (log)
  const float dxyOverPt = std::abs(dxy) / std::max(pt, kEpsilon);
  features.push_back(std::log10(dxyOverPt + kEpsilon));

  // 15: ptErr / p (log)
  const float ptErrOverP = ptErr / std::max(p, kEpsilon);
  features.push_back(std::log10(ptErrOverP + kEpsilon));

  // 16: |dz| / |dxy| topology ratio (log)
  const float dzOverDxy = std::abs(dz) / (std::abs(dxy) + kEpsilon);
  features.push_back(std::log10(dzOverDxy + kEpsilon));

  // 17: |eta|
  features.push_back(std::abs(eta));

  // ---------------------------------------------------------------------
  // L1TkMuon matching block (features 18-32)
  // ---------------------------------------------------------------------
  if (useL1TkMuFeatures_) {
    // Best L1 match search
    float minDR2 = std::numeric_limits<float>::max();
    float matchedL1Pt = -1.0f;
    int bestIndex = -1;

    // Loose-window count for feature 31
    int nCompatible = 0;

    // PASS 1: Find best match (and count loose-compatible candidates)
    for (size_t l1Idx = 0; l1Idx != l1TkMuons.size(); ++l1Idx) {
      const auto& l1TkMu = l1TkMuons.at(l1Idx);

      // Use propagated muon-system kinematics (phEta/phPhi/phPt) - same as training n-tuple.
      const float l1Eta = l1TkMu.phEta();
      const float l1Phi = l1TkMu.phPhi();
      const float l1Pt = l1TkMu.phPt();

      const float ptDiff = pt - l1Pt;
      const float chi2Pt = (ptDiff * ptDiff) / (ptErr * ptErr + kChi2PtEps);
      const float dR2 = reco::deltaR2(eta, phi, l1Eta, l1Phi);

      // Count loosely compatible L1 candidates (feature 31).
      if (dR2 < kLooseDR2Cut && chi2Pt < kLooseChi2PtCut) {
        nCompatible++;
      }

      // Strict compatibility filter (chi2Pt < 9.0) for best-match search.
      if (chi2Pt >= kMatchChi2PtCut)
        continue;

      if (dR2 < minDR2) {
        minDR2 = dR2;
        matchedL1Pt = l1Pt;
        bestIndex = static_cast<int>(l1Idx);
      }
    }

    // PASS 2: Find second-best dR2 (strictly greater than the best, to mirror the
    // Python semantics where exact ties on dR2 are all marked as "best" and thus
    // excluded from the second-best computation).
    float secondBestDR2 = kSentinel;
    if (bestIndex >= 0) {
      for (size_t l1Idx = 0; l1Idx != l1TkMuons.size(); ++l1Idx) {
        if (static_cast<int>(l1Idx) == bestIndex)
          continue;

        const auto& l1TkMu = l1TkMuons.at(l1Idx);
        const float l1Eta = l1TkMu.phEta();
        const float l1Phi = l1TkMu.phPhi();
        const float l1Pt = l1TkMu.phPt();

        const float ptDiff = pt - l1Pt;
        const float chi2Pt = (ptDiff * ptDiff) / (ptErr * ptErr + kChi2PtEps);
        if (chi2Pt >= kMatchChi2PtCut)
          continue;

        const float dR2 = reco::deltaR2(eta, phi, l1Eta, l1Phi);
        if (dR2 > minDR2 && dR2 < secondBestDR2) {
          secondBestDR2 = dR2;
        }
      }
    }

    // Compute matching quantities (only if we found a strictly compatible L1).
    float dPtNorm = kImputeDPtNorm;
    float chi2PtBest = kImputeChi2Pt;
    float matchingScore = kImputeMatchScore;

    if (bestIndex >= 0) {
      dPtNorm = std::abs(pt - matchedL1Pt) / (matchedL1Pt + kDPtNormEps);
      const float ptDiffBest = pt - matchedL1Pt;
      chi2PtBest = (ptDiffBest * ptDiffBest) / (ptErr * ptErr + kChi2PtEps);
      matchingScore = minDR2 * (1.0f + dPtNorm);
    }

    const bool hasL1Match = (minDR2 < kMatchDR2Cut) && (bestIndex >= 0);

    // --- Stub features (18-24) ---
    if (useStubFeatures_) {
      if (hasL1Match) {
        const auto& bestL1 = l1TkMuons[bestIndex];

        // Count only non-null stubs so nStubsTotal == nEndcap + nBarrel + nOther,
        // matching Python's ak.sum(is_stub_for_l1, axis=2).
        int nStubsTotal = 0;
        int nStubsEndcap = 0;
        int nStubsBarrel = 0;
        int maxStubQuality = 0;
        int minDepthRegion = std::numeric_limits<int>::max();
        int bestStubIndex = -1;

        for (size_t s = 0; s != bestL1.stubs().size(); ++s) {
          const auto stubRef = bestL1.stubs().at(s);
          if (stubRef.isNull())
            continue;

          ++nStubsTotal;
          if (stubRef->type() == 0)
            ++nStubsEndcap;
          else if (stubRef->type() == 1)
            ++nStubsBarrel;

          // Best stub: highest quality; on equal quality, smallest depthRegion.
          // First-seen wins on full ties (matches Python ak.firsts on the masked array).
          if (stubRef->quality() > maxStubQuality ||
              (stubRef->quality() == maxStubQuality && stubRef->depthRegion() < minDepthRegion)) {
            maxStubQuality = stubRef->quality();
            minDepthRegion = stubRef->depthRegion();
            bestStubIndex = static_cast<int>(s);
          }
        }
        features.push_back(static_cast<float>(nStubsTotal));   // 18
        features.push_back(static_cast<float>(nStubsEndcap));  // 19
        features.push_back(static_cast<float>(nStubsBarrel));  // 20

        if (bestStubIndex >= 0) {
          const auto bestStub = bestL1.stubs().at(bestStubIndex);
          features.push_back(static_cast<float>(bestStub->quality()));      // 21
          features.push_back(static_cast<float>(bestStub->etaRegion()));    // 22
          features.push_back(static_cast<float>(bestStub->phiRegion()));    // 23
          features.push_back(static_cast<float>(bestStub->depthRegion()));  // 24
        } else {
          // L1 matched but no usable stubs - Python: maxQual fill_none - 0,
          // best-stub region fill_none - -1.
          features.push_back(0.0f);   // 21 stubQual_max
          features.push_back(-1.0f);  // 22 stubMax_etaRegion
          features.push_back(-1.0f);  // 23 stubMax_phiRegion
          features.push_back(-1.0f);  // 24 stubMax_depthRegion
        }
      } else {
        features.push_back(0.0f);   // 18 nStubs
        features.push_back(0.0f);   // 19 nStubs_Endcap
        features.push_back(0.0f);   // 20 nStubs_Barrel
        features.push_back(0.0f);   // 21 stubQual_max
        features.push_back(-1.0f);  // 22 stubMax_etaRegion
        features.push_back(-1.0f);  // 23 stubMax_phiRegion
        features.push_back(-1.0f);  // 24 stubMax_depthRegion
      }
    }

    // --- L1TkMu matching features (25-29) ---
    features.push_back(hasL1Match ? 1.0f : 0.0f);  // 25: hasMatch

    if (hasL1Match) {
      features.push_back(std::log10(std::abs(minDR2) + kEpsilon));         // 26: dR2min
      features.push_back(std::log10(std::abs(dPtNorm) + kEpsilon));        // 27: dPtNorm
      features.push_back(std::log10(std::abs(chi2PtBest) + kEpsilon));     // 28: chi2Pt
      features.push_back(std::log10(std::abs(matchingScore) + kEpsilon));  // 29: matchingScore
    } else {
      features.push_back(std::log10(std::abs(kImputeDR2) + kEpsilon));         // 26
      features.push_back(std::log10(std::abs(kImputeDPtNorm) + kEpsilon));     // 27
      features.push_back(std::log10(std::abs(kImputeChi2Pt) + kEpsilon));      // 28
      features.push_back(std::log10(std::abs(kImputeMatchScore) + kEpsilon));  // 29
    }

    // --- NEW L1 matching features (30-31) ---
    // 30: nCompatible - number of L1 candidates within loose window
    features.push_back(static_cast<float>(nCompatible));

    // 31: secondBest_dR2 (log, imputed)
    //     hasSecond is true only if we found a real second-best dR2 in the strict window.
    const bool hasSecond = secondBestDR2 < (kSentinel - 1.0f);
    if (hasSecond) {
      features.push_back(std::log10(std::abs(secondBestDR2) + kEpsilon));
    } else {
      features.push_back(std::log10(std::abs(kImputeSecondDR2) + kEpsilon));
    }
  }

  // Feature 32: Low pT indicator
  float exponent = (pt - 5.0f) * 2.0f;
  exponent = std::clamp(exponent, -20.0f, 20.0f);
  const float lowPtIndicator = 1.0f / (1.0f + std::exp(exponent));
  features.push_back(lowPtIndicator);

  return features;
}

void MuonIOTracksForestSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonIOTracksForestSelector";

  auto selectedTracks = std::make_unique<reco::TrackCollection>();
  auto scores = std::make_unique<std::vector<float>>();

  // Get input collections
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<l1t::TrackerMuonCollection> l1TkMuons;
  iEvent.getByToken(l1TkMuonsToken_, l1TkMuons);

  if (!tracks.isValid() || tracks->empty()) {
    iEvent.put(std::move(selectedTracks));
    iEvent.put(std::move(scores), "scores");
    return;
  }

  const auto* cache = globalCache();

  const unsigned int evtIdx = eventCounter_++;
  for (size_t i = 0; i < tracks->size(); ++i) {
    const auto& track = (*tracks)[i];
    auto features = extractFeatures(track, *l1TkMuons);
    if (static_cast<int>(features.size()) != nFeatures_) {
      throw cms::Exception("MuonIOTracksForestSelector")
          << "Feature count mismatch: extracted " << features.size() << " features, expected " << nFeatures_
          << ". Check useStubFeatures/useL1TkMuFeatures vs the trained forest model.";
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

    if (prob >= decisionThreshold_) {
      selectedTracks->push_back(*(reco::TrackRef(tracks, i)));
    }
  }

  std::cout << metname << " Selected " << selectedTracks->size() << " out of " << tracks->size() << " tracks\n";

  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(scores), "scores");
}

void MuonIOTracksForestSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks", edm::InputTag("hltPhase2MuonPixelTracks"))->setComment("Input track collection");
  desc.add<edm::InputTag>("l1TkMuons", edm::InputTag("l1tTkMuonsGmt"))
      ->setComment("L1 Tracker Muon collection for matching features");
  desc.add<std::string>("modelPath", "RecoMuon/L3TrackFinder/data/pixel_track_selector_forest.bin")
      ->setComment("Path to compact gradient-boosted-tree binary (.bin)");
  desc.add<double>("decisionThreshold", 0.5)
      ->setComment("Probability threshold for track selection (use F2-optimal from training)");
  desc.add<bool>("useL1TkMuFeatures", true)->setComment("Include L1 Tracker Muon matching features");
  desc.add<bool>("useStubFeatures", true)->setComment("Include stub-related features (requires stub info in event)");
  desc.add<int>("nFeatures", 33)
      ->setComment(
          "Total number of input features for the pruned 33-feature model "
          "(7 log + 3 plain + 7 derived + 1 absEta + 7 stub + 5 L1 match + 2 new L1 + 1 low-pT = 33)");
  desc.addUntracked<bool>("dumpFeatures", false)
      ->setComment("Print one CSV-like line per track to stdout, for cross-validation against build_dataset().");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonIOTracksForestSelector);
