#ifndef RecoMuon_L3TrackFinder_OITrackSelectorFeatures_h
#define RecoMuon_L3TrackFinder_OITrackSelectorFeatures_h

// Feature extraction for the muon outside-in (OI) track high-purity selectors
// (pixel-chain and seeds/general-chain deployments), used by
// MuonOITracksForestSelector; the feature set is identical to the retired
// DNN selector it replaces.
//
// The feature set is the 22-feature OI production set (round-2 pruned):
// ptErr, chi2, sigmaPtOverPt and relUncertaintyProduct of the original
// 26-feature layout are NOT part of the trained ABI and are not extracted.
// The struct gives every feature a name; toArray() returns them in the
// canonical training order used by the compact-forest/ONNX inputs.
// Training-side counterpart: the build_dataset() feature extraction of the
// OI XGBoost trainers.
//
// Matching semantics (must equal the python extraction, pixel-level):
// wrapped deltaPhi; score = (chi2Eta + chi2Phi + chi2Pt + chi2Dz) / 9
// (kMatchChi2 normalization); hasMatch = bestScore < 25; matchingScore =
// log10(bestScore + eps) when matched, else kImputeMatchScore in RAW space.

#include <array>
#include <cmath>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoMuon/L3TrackFinder/interface/IOTrackSelectorFeatures.h"  // muonhp::kEpsilon, kLowPtCut

namespace muonhp {

  // Matching/imputation constants - must equal the Python training-side
  // feature extraction.
  inline constexpr float kOIMatchChi2 = 9.0f;
  inline constexpr float kOINoMatchBestScore = 25.0f;   // 5 sigma^2 in /9 units
  inline constexpr float kOIImputeMatchScore = 10.0f;   // RAW space (not log10)
  inline constexpr float kOIStandaloneChi2Eps = 1e-12f;

  // The 22-feature OI production set, in canonical (training) order.
  struct OITrackFeatures {
    // Log-compressed track parameters
    float logP = 0.f;
    float logPt = 0.f;
    float logNormalizedChi2 = 0.f;
    float logEtaErr = 0.f;
    float logPhiErr = 0.f;
    float logDszErr = 0.f;
    float logDxyErr = 0.f;
    float logDzErr = 0.f;
    float logQoverpErr = 0.f;
    float logLambdaErr = 0.f;
    // Plain kinematics + hit content
    float eta = 0.f;
    float nPixelHits = 0.f;
    float nTrkLays = 0.f;
    float nFoundHits = 0.f;
    float nLostHits = 0.f;
    // Derived track-quality features
    float logImpact3DSq = 0.f;         // log10(dxy^2 + dz^2 + eps)
    float logImpactSignificance = 0.f;
    float logChi2PerHit = 0.f;
    float hitEfficiency = 0.f;
    // Standalone (L2 muon vertex) matching
    float hasStandaloneMatch = 0.f;
    float logStandaloneMatchScore = 0.f;
    // Regime indicator (soft sigmoid around kLowPtCut, from the IO header)
    float lowPtSigmoid = 0.f;

    static constexpr size_t kSize = 22;

    // Canonical ordering consumed by the deployed models. The order below is
    // the training ABI and must not change without retraining.
    std::array<float, kSize> toArray() const {
      return {// 0-9 log track parameters
              logP, logPt, logNormalizedChi2, logEtaErr, logPhiErr, logDszErr, logDxyErr, logDzErr,
              logQoverpErr, logLambdaErr,
              // 10-14 plain kinematics + hit content
              eta, nPixelHits, nTrkLays, nFoundHits, nLostHits,
              // 15-18 derived track quality
              logImpact3DSq, logImpactSignificance, logChi2PerHit, hitEfficiency,
              // 19-20 standalone matching
              hasStandaloneMatch, logStandaloneMatchScore,
              // 21 regime indicator
              lowPtSigmoid};
    }
  };

  inline OITrackFeatures extractOITrackFeatures(const reco::Track& track,
                                                const reco::TrackCollection& standaloneMuons) {
    OITrackFeatures f;

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
    const float qoverpErr = track.qoverpError();
    const float lambdaErr = track.lambdaError();

    f.logP = std::log10(std::abs(p) + kEpsilon);
    f.logPt = std::log10(std::abs(pt) + kEpsilon);
    f.logNormalizedChi2 = std::log10(std::abs(normalizedChi2) + kEpsilon);
    f.logEtaErr = std::log10(std::abs(etaErr) + kEpsilon);
    f.logPhiErr = std::log10(std::abs(phiErr) + kEpsilon);
    f.logDszErr = std::log10(std::abs(dszErr) + kEpsilon);
    f.logDxyErr = std::log10(std::abs(dxyErr) + kEpsilon);
    f.logDzErr = std::log10(std::abs(dzErr) + kEpsilon);
    f.logQoverpErr = std::log10(std::abs(qoverpErr) + kEpsilon);
    f.logLambdaErr = std::log10(std::abs(lambdaErr) + kEpsilon);

    f.eta = eta;
    f.nPixelHits = static_cast<float>(track.hitPattern().numberOfValidPixelHits());
    f.nTrkLays = static_cast<float>(track.hitPattern().trackerLayersWithMeasurement());
    f.nFoundHits = static_cast<float>(track.numberOfValidHits());
    f.nLostHits = static_cast<float>(track.numberOfLostHits());

    // Derived track-quality features
    f.logImpact3DSq = std::log10(dxy * dxy + dz * dz + kEpsilon);
    const float dxySignificance = dxy / std::max(dxyErr, kEpsilon);
    const float dzSignificance = dz / std::max(dzErr, kEpsilon);
    f.logImpactSignificance =
        std::log10(std::sqrt(dxySignificance * dxySignificance + dzSignificance * dzSignificance) + kEpsilon);
    const int nFound = track.numberOfValidHits();
    const int nLost = track.numberOfLostHits();
    f.logChi2PerHit = std::log10(chi2 / std::max(nFound, 1) + kEpsilon);
    f.hitEfficiency = static_cast<float>(nFound) / std::max(nFound + nLost, 1);

    // Standalone (L2 muon vertex) matching
    float bestScore = kOINoMatchBestScore;
    for (const auto& muon : standaloneMuons) {
      const float muEtaErr = muon.etaError();
      const float muPhiErr = muon.phiError();
      const float muPtErr = muon.ptError();
      const float muDzErr = muon.dzError();

      const float chi2Eta = (eta - muon.eta()) * (eta - muon.eta()) /
                          (etaErr * etaErr + muEtaErr * muEtaErr + kOIStandaloneChi2Eps);
      const float dPhi = reco::deltaPhi(phi, muon.phi());
      const float chi2Phi = dPhi * dPhi / (phiErr * phiErr + muPhiErr * muPhiErr + kOIStandaloneChi2Eps);
      const float chi2Pt =
          (pt - muon.pt()) * (pt - muon.pt()) / (ptErr * ptErr + muPtErr * muPtErr + kOIStandaloneChi2Eps);
      const float chi2Dz =
          (dz - muon.dz()) * (dz - muon.dz()) / (dzErr * dzErr + muDzErr * muDzErr + kOIStandaloneChi2Eps);

      const float score = (chi2Eta + chi2Phi + chi2Pt + chi2Dz) / kOIMatchChi2;
      if (score < bestScore)
        bestScore = score;
    }
    if (bestScore < kOINoMatchBestScore) {
      f.hasStandaloneMatch = 1.0f;
      f.logStandaloneMatchScore = std::log10(bestScore + kEpsilon);
    } else {
      f.hasStandaloneMatch = 0.0f;
      f.logStandaloneMatchScore = kOIImputeMatchScore;  // RAW space (matches training)
    }

    // Regime indicator
    const float exponent = std::clamp((pt - kLowPtCut) * 2.0f, -20.0f, 20.0f);
    f.lowPtSigmoid = 1.0f / (1.0f + std::exp(exponent));

    return f;
  }

}  // namespace muonhp

#endif
