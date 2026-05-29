/** \class MuonIOTracksDNNSelector.cc
 *
 *  \brief DNN-based muon (pixel) tracks selection
 *
 *  This module implements a DNN-based selector for muon tracks, designed to replace
 *  traditional cut-based selections with a more powerful multivariate approach.
 *  The DNN model is expected to be trained on a comprehensive set of track features,
 *  including kinematic variables, track quality metrics, and optionally L1Tk Muon
 *  matching information and stub features.
 *
 *  \author Luca Ferragina (INFN BO), 2026
 *
 * Feature order (44 total):
 *
 * Log features (log10(|x| + 1e-6)):
 *   0: track_p
 *   1: track_pt
 *   2: track_ptErr
 *   3: track_chi2
 *   4: track_normalizedChi2
 *   5: track_etaErr
 *   6: track_phiErr
 *   7: track_dszErr
 *   8: track_dxyErr
 *   9: track_dzErr
 *  10: track_qoverpErr
 *  11: track_lambdaErr
 *
 * Plain features (no transform):
 *  12: track_eta
 *  13: track_nPixelHits
 *  14: track_nTrkLays
 *  15: track_nFoundHits
 *  16: track_nLostHits
 *
 * Derived features (original):
 *  17: track_impact3D              (log10(dxy^2 + dz^2 + eps))
 *  18: track_impactSignificance    (log10(sqrt((dxy/dxyErr)^2 + (dz/dzErr)^2) + eps))
 *  19: track_chi2PerHit            (log10(chi2 / max(nFoundHits, 1) + eps))
 *  20: track_hitEfficiency         (nFoundHits / max(nFoundHits + nLostHits, 1))
 *  21: track_sigmaPtOverPt         (log10(ptErr / pt + eps))
 *  22: track_relUncertaintyProduct (log10((ptErr/pt) * (qoverpErr/|qoverp|) + eps))
 *
 * NEW derived features:
 *  23: track_sip2D                 (log10(|dxy| / dxyErr + eps))
 *  24: track_sipZ                  (log10(|dz| / dzErr + eps))
 *  25: track_dxyOverPt             (log10(|dxy| / pt + eps))
 *  26: track_ptErrOverP            (log10(ptErr / p + eps))
 *  27: track_dzOverDxy             (log10(|dz| / (|dxy| + eps) + eps))
 *  28: track_absEta                (|eta|)
 *
 * L1TkMuon stub features (if useStubFeatures=true):
 *  29: nStubs
 *  30: nStubs_Endcap
 *  31: nStubs_Barrel
 *  32: stubQual_max
 *  33: stubMax_etaRegion
 *  34: stubMax_phiRegion
 *  35: stubMax_depthRegion
 *
 * L1TkMuon matching features (if useL1TkMuFeatures=true):
 *  36: L1TkMu_hasMatch             (1.0 if matched, 0.0 otherwise)
 *  37: L1TkMu_dR2min               (log, imputed with 0.1)
 *  38: L1TkMu_dPtNorm              (log, imputed with 1.0)
 *  39: L1TkMu_chi2Pt               (log, imputed with 10.0)
 *  40: L1TkMu_matchingScore        (log, imputed with 0.2)
 *
 * NEW L1TkMuon matching features:
 *  41: L1TkMu_nCompatible          (count of L1 candidates within loose window)
 *  42: L1TkMu_secondBest_dR2       (log, imputed with 1.0)
 *
 * Low pT indicator:
 *  43: is_low_pt                   (1 / (1 + exp(clip((pt - 5.0) * 2.0, -20, 20))))
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

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include <cmath>
#include <algorithm>
#include <limits>

class MuonIOTracksDNNSelector : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>> {
public:
  explicit MuonIOTracksDNNSelector(const edm::ParameterSet&, const cms::Ort::ONNXRuntime*);
  ~MuonIOTracksDNNSelector() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);
  static std::unique_ptr<cms::Ort::ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const cms::Ort::ONNXRuntime*);

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
  static constexpr float kEpsilon = 1e-6f;     // generic log/division floor (matches Python 1e-6)
  static constexpr float kChi2PtEps = 1e-12f;  // matches Python (t_ptErr**2 + 1e-12)
  static constexpr float kDPtNormEps = 1e-9f;  // matches Python (l1_pt + 1e-9)

  // Imputation values for non-matched L1 features
  static constexpr float kImputeDR2 = 0.1f;
  static constexpr float kImputeDPtNorm = 1.0f;
  static constexpr float kImputeChi2Pt = 10.0f;
  static constexpr float kImputeMatchScore = 0.2f;
  static constexpr float kImputeSecondDR2 = 1.0f;
};

MuonIOTracksDNNSelector::MuonIOTracksDNNSelector(const edm::ParameterSet& iConfig, const cms::Ort::ONNXRuntime* cache)
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

std::unique_ptr<cms::Ort::ONNXRuntime> MuonIOTracksDNNSelector::initializeGlobalCache(const edm::ParameterSet& iConfig) {
  edm::FileInPath modelPath(iConfig.getParameter<std::string>("modelPath"));
  return std::make_unique<cms::Ort::ONNXRuntime>(modelPath.fullPath());
}

void MuonIOTracksDNNSelector::globalEndJob(const cms::Ort::ONNXRuntime* cache) {}

std::vector<float> MuonIOTracksDNNSelector::extractFeatures(const reco::Track& track,
                                                            const l1t::TrackerMuonCollection& l1TkMuons) const {
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
  const float chi2 = track.chi2();
  const float normalizedChi2 = track.normalizedChi2();

  const float dszErr = track.dszError();
  const float dxy = track.dxy();
  const float dxyErr = track.dxyError();
  const float dz = track.dz();
  const float dzErr = track.dzError();
  const float qoverp = track.qoverp();
  const float qoverpErr = track.qoverpError();
  const float lambdaErr = track.lambdaError();

  const int nPixelHits = track.hitPattern().numberOfValidPixelHits();
  const int nTrkLays = track.hitPattern().trackerLayersWithMeasurement();
  const int nFoundHits = track.numberOfValidHits();
  const int nLostHits = track.numberOfLostHits();

  // Features 0-11: Log features
  features.push_back(std::log10(std::abs(p) + kEpsilon));               // 0
  features.push_back(std::log10(std::abs(pt) + kEpsilon));              // 1
  features.push_back(std::log10(std::abs(ptErr) + kEpsilon));           // 2
  features.push_back(std::log10(std::abs(chi2) + kEpsilon));            // 3
  features.push_back(std::log10(std::abs(normalizedChi2) + kEpsilon));  // 4
  features.push_back(std::log10(std::abs(etaErr) + kEpsilon));          // 5
  features.push_back(std::log10(std::abs(phiErr) + kEpsilon));          // 6
  features.push_back(std::log10(std::abs(dszErr) + kEpsilon));          // 7
  features.push_back(std::log10(std::abs(dxyErr) + kEpsilon));          // 8
  features.push_back(std::log10(std::abs(dzErr) + kEpsilon));           // 9
  features.push_back(std::log10(std::abs(qoverpErr) + kEpsilon));       // 10
  features.push_back(std::log10(std::abs(lambdaErr) + kEpsilon));       // 11

  // Features 12-16: Plain features
  features.push_back(eta);                             // 12
  features.push_back(static_cast<float>(nPixelHits));  // 13
  features.push_back(static_cast<float>(nTrkLays));    // 14
  features.push_back(static_cast<float>(nFoundHits));  // 15
  features.push_back(static_cast<float>(nLostHits));   // 16

  // Features 17-22: Original derived features
  // 17: Impact Parameter 3D (log)
  const float ip3d = dxy * dxy + dz * dz;
  features.push_back(std::log10(ip3d + kEpsilon));

  // 18: Combined Impact Significance (log)
  const float sipDxy = dxy / std::max(dxyErr, kEpsilon);
  const float sipDz = dz / std::max(dzErr, kEpsilon);
  const float sipCombined = std::sqrt(sipDxy * sipDxy + sipDz * sipDz);
  features.push_back(std::log10(sipCombined + kEpsilon));

  // 19: Chi2 per hit (log)
  const float chi2PerHit = chi2 / std::max(nFoundHits, 1);
  features.push_back(std::log10(chi2PerHit + kEpsilon));

  // 20: Hit Efficiency (linear)
  const float hitEfficiency = static_cast<float>(nFoundHits) / std::max(nFoundHits + nLostHits, 1);
  features.push_back(hitEfficiency);

  // 21: SigmaPt / Pt (log)
  const float sigmaPtOverPt = ptErr / std::max(pt, kEpsilon);
  features.push_back(std::log10(sigmaPtOverPt + kEpsilon));

  // 22: Relative Uncertainty Product (log)
  const float relUncertaintyProduct = sigmaPtOverPt * (qoverpErr / std::max(std::abs(qoverp), kEpsilon));
  features.push_back(std::log10(relUncertaintyProduct + kEpsilon));

  // Features 23-28: NEW derived features
  // 23: Separated 2D impact parameter significance (log)
  const float sip2D = std::abs(dxy) / std::max(dxyErr, kEpsilon);
  features.push_back(std::log10(sip2D + kEpsilon));

  // 24: Longitudinal impact parameter significance (log)
  const float sipZ = std::abs(dz) / std::max(dzErr, kEpsilon);
  features.push_back(std::log10(sipZ + kEpsilon));

  // 25: |dxy| / pT (log)
  const float dxyOverPt = std::abs(dxy) / std::max(pt, kEpsilon);
  features.push_back(std::log10(dxyOverPt + kEpsilon));

  // 26: ptErr / p (log) - at large |eta|, p >> pt, so ptErr/p reveals forward-track quality
  const float ptErrOverP = ptErr / std::max(p, kEpsilon);
  features.push_back(std::log10(ptErrOverP + kEpsilon));

  // 27: |dz| / |dxy| topology ratio (log)
  //     Note: Python uses (|dxy| + eps) in the denominator, NOT max(|dxy|, eps).
  //     For |dxy| ~ 0 this gives 1/eps = 1e6, floored by the outer log10(... + eps).
  const float dzOverDxy = std::abs(dz) / (std::abs(dxy) + kEpsilon);
  features.push_back(std::log10(dzOverDxy + kEpsilon));

  // 28: |eta|
  features.push_back(std::abs(eta));

  // ---------------------------------------------------------------------
  // L1TkMuon matching block (features 29-42)
  // ---------------------------------------------------------------------
  if (useL1TkMuFeatures_) {
    // Best L1 match search
    float minDR2 = std::numeric_limits<float>::max();
    float matchedL1Pt = -1.0f;
    int bestIndex = -1;

    // Loose-window count for feature 41
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

      // Count loosely compatible L1 candidates (feature 41).
      // In the Python training is_loose_compatible is computed on ALL L1 candidates
      // without any pre-filtering. The strict chi2Pt < 9 cut only applies to the
      // best-match search and the second-best dR2 below.
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
        // Strict ">" on minDR2: any L1 sharing the exact best dR2 is excluded
        // from being second-best, matching Python's is_best_match treatment of ties.
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

    // --- Stub features (29-35) ---
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

        features.push_back(static_cast<float>(nStubsTotal));   // 29
        features.push_back(static_cast<float>(nStubsEndcap));  // 30
        features.push_back(static_cast<float>(nStubsBarrel));  // 31

        if (bestStubIndex >= 0) {
          const auto bestStub = bestL1.stubs().at(bestStubIndex);
          features.push_back(static_cast<float>(bestStub->quality()));      // 32
          features.push_back(static_cast<float>(bestStub->etaRegion()));    // 33
          features.push_back(static_cast<float>(bestStub->phiRegion()));    // 34
          features.push_back(static_cast<float>(bestStub->depthRegion()));  // 35
        } else {
          // L1 matched but no usable stubs - Python: maxQual fill_none - 0,
          // best-stub region fill_none - -1.
          features.push_back(0.0f);   // 32 stubQual_max
          features.push_back(-1.0f);  // 33 stubMax_etaRegion
          features.push_back(-1.0f);  // 34 stubMax_phiRegion
          features.push_back(-1.0f);  // 35 stubMax_depthRegion
        }
      } else {
        features.push_back(0.0f);   // 29 nStubs
        features.push_back(0.0f);   // 30 nStubs_Endcap
        features.push_back(0.0f);   // 31 nStubs_Barrel
        features.push_back(0.0f);   // 32 stubQual_max
        features.push_back(-1.0f);  // 33 stubMax_etaRegion
        features.push_back(-1.0f);  // 34 stubMax_phiRegion
        features.push_back(-1.0f);  // 35 stubMax_depthRegion
      }
    }

    // --- L1TkMu matching features (36-40) ---
    features.push_back(hasL1Match ? 1.0f : 0.0f);  // 36: hasMatch

    if (hasL1Match) {
      features.push_back(std::log10(std::abs(minDR2) + kEpsilon));         // 37: dR2min
      features.push_back(std::log10(std::abs(dPtNorm) + kEpsilon));        // 38: dPtNorm
      features.push_back(std::log10(std::abs(chi2PtBest) + kEpsilon));     // 39: chi2Pt
      features.push_back(std::log10(std::abs(matchingScore) + kEpsilon));  // 40: matchingScore
    } else {
      features.push_back(std::log10(std::abs(kImputeDR2) + kEpsilon));         // 37
      features.push_back(std::log10(std::abs(kImputeDPtNorm) + kEpsilon));     // 38
      features.push_back(std::log10(std::abs(kImputeChi2Pt) + kEpsilon));      // 39
      features.push_back(std::log10(std::abs(kImputeMatchScore) + kEpsilon));  // 40
    }

    // --- NEW L1 matching features (41-42) ---

    // 41: nCompatible - number of L1 candidates within loose window
    features.push_back(static_cast<float>(nCompatible));

    // 42: secondBest_dR2 (log, imputed)
    //     hasSecond is true only if we found a real second-best dR2 in the strict window.
    const bool hasSecond = secondBestDR2 < (kSentinel - 1.0f);
    if (hasSecond) {
      features.push_back(std::log10(std::abs(secondBestDR2) + kEpsilon));
    } else {
      features.push_back(std::log10(std::abs(kImputeSecondDR2) + kEpsilon));
    }
  }

  // Feature 43: Low pT indicator
  float exponent = (pt - 5.0f) * 2.0f;
  exponent = std::clamp(exponent, -20.0f, 20.0f);
  const float lowPtIndicator = 1.0f / (1.0f + std::exp(exponent));
  features.push_back(lowPtIndicator);

  return features;
}

void MuonIOTracksDNNSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonIOTracksDNNSelector";

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

  // Prepare batch input
  std::vector<float> inputData;
  inputData.reserve(tracks->size() * nFeatures_);

  const unsigned int evtIdx = eventCounter_++;
  for (size_t i = 0; i < tracks->size(); ++i) {
    const auto& track = (*tracks)[i];
    auto features = extractFeatures(track, *l1TkMuons);
    if (static_cast<int>(features.size()) != nFeatures_) {
      throw cms::Exception("MuonIOTracksDNNSelector")
          << "Feature count mismatch: extracted " << features.size() << " features, expected " << nFeatures_
          << ". Check useStubFeatures/useL1TkMuFeatures vs the trained ONNX model.";
    }
    inputData.insert(inputData.end(), features.begin(), features.end());

    // Optional dump for cross-validation
    if (dumpFeatures_) {
      std::ostringstream oss;
      oss << evtIdx << "," << i;
      oss << std::scientific << std::setprecision(9);
      for (float f : features) {
        oss << "," << f;
      }
      std::cout << oss.str() << "\n";
    }
  }

  // Run inference
  std::vector<std::vector<int64_t>> inputShapes = {
      {static_cast<int64_t>(tracks->size()), static_cast<int64_t>(nFeatures_)}};
  cms::Ort::FloatArrays inputTensor({inputData});

  auto outputs = globalCache()->run({"input"}, inputTensor, inputShapes, {"output"}, 1);

  // Model outputs probabilities directly (sigmoid is in the network)
  const auto& probs = outputs[0];

  LogDebug(metname) << "Processing " << tracks->size() << " tracks with threshold " << decisionThreshold_;

  for (size_t i = 0; i < tracks->size(); ++i) {
    float prob = probs[i];
    LogDebug(metname) << "  Track " << i << ": DNN score = " << prob;

    // Clamp probability to valid range (safety check)
    prob = std::clamp(prob, 0.0f, 1.0f);

    scores->push_back(prob);

    if (prob >= decisionThreshold_) {
      LogDebug(metname) << "    -> Selected";
      selectedTracks->push_back(*(reco::TrackRef(tracks, i)));
    }
  }

  std::cout << "Selected " << selectedTracks->size() << " out of " << tracks->size() << " tracks\n";

  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(scores), "scores");
}

void MuonIOTracksDNNSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("tracks", edm::InputTag("hltPhase2MuonPixelTracks"))->setComment("Input track collection");
  desc.add<edm::InputTag>("l1TkMuons", edm::InputTag("l1tTkMuonsGmt"))
      ->setComment("L1 Tracker Muon collection for matching features");
  desc.add<std::string>("modelPath", "RecoMuon/L3TrackFinder/data/pixel_track_selector.onnx")
      ->setComment("Path to ONNX model file (expects raw unscaled inputs, scaler fused)");
  desc.add<double>("decisionThreshold", 0.5)->setComment("Probability threshold for track selection");
  desc.add<bool>("useL1TkMuFeatures", true)->setComment("Include L1 Tracker Muon matching features");
  desc.add<bool>("useStubFeatures", true)->setComment("Include stub-related features (requires stub info in event)");
  desc.add<int>("nFeatures", 44)
      ->setComment(
          "Total number of input features "
          "(17 base + 6 original derived + 6 new derived + 7 stub + 5 L1 match + 2 new L1 + 1 low-pT + 1 raw-pT = 45)");
  desc.addUntracked<bool>("dumpFeatures", false)
      ->setComment("Print one CSV-like line per track to stdout, for cross-validation against build_dataset().");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonIOTracksDNNSelector);
