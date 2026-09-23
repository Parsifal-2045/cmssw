#ifndef RecoMuon_L3TrackFinder_OITrackSelectorFeatures_h
#define RecoMuon_L3TrackFinder_OITrackSelectorFeatures_h

// Feature extraction for the muon outside-in (OI) track high-purity selectors
// (pixel-chain and seeds/general-chain deployments), used by
// MuonOITracksForestSelector.
//
// The feature set is the 22-feature OI production set (round-2 pruned):
// ptErr, chi2, sigmaPtOverPt and relUncertaintyProduct of the original
// 26-feature layout are NOT part of the trained ABI and are not extracted.
// The struct gives every feature a name; toArray() returns them in the
// canonical training order used by the compact-forest/ONNX inputs.
// Training-side counterpart: build_dataset() + OI_PRODUCTION_FEATURES in
// muonHighPurityTrackSelection/production/oi/OI_features.py (the training
// pipeline asserts that its kept features follow this order).
//
// Matching semantics (must equal the python extraction, pixel-level):
// wrapped deltaPhi; score = (chi2Eta + chi2Phi + chi2Pt + chi2Dz) / 9
// (kMatchChi2 normalization); hasMatch = bestScore < 25; matchingScore =
// log10(bestScore + eps) when matched, else kImputeMatchScore in RAW space.
//
// Numeric convention as in IOTrackSelectorFeatures.h: raw inputs rounded to
// float (as stored in the training n-tuple), double-precision arithmetic,
// each feature rounded to float once.

#include <algorithm>
#include <array>
#include <cmath>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoMuon/L3TrackFinder/interface/IOTrackSelectorFeatures.h"  // kEpsilon, kLowPtCut, asStored, logFloor

namespace muonhp {

  // Matching/imputation constants - must equal the Python training-side
  // feature extraction.
  inline constexpr double kOIMatchChi2 = 9.0;
  inline constexpr double kOINoMatchBestScore = 25.0;  // 5 sigma^2 in /9 units
  inline constexpr float kOIImputeMatchScore = 10.0f;  // RAW space (not log10)
  inline constexpr double kOIStandaloneChi2Eps = 1e-12;

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
    float logImpact3DSq = 0.f;  // log10(dxy^2 + dz^2 + eps)
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
              logP,
              logPt,
              logNormalizedChi2,
              logEtaErr,
              logPhiErr,
              logDszErr,
              logDxyErr,
              logDzErr,
              logQoverpErr,
              logLambdaErr,
              // 10-14 plain kinematics + hit content
              eta,
              nPixelHits,
              nTrkLays,
              nFoundHits,
              nLostHits,
              // 15-18 derived track quality
              logImpact3DSq,
              logImpactSignificance,
              logChi2PerHit,
              hitEfficiency,
              // 19-20 standalone matching
              hasStandaloneMatch,
              logStandaloneMatchScore,
              // 21 regime indicator
              lowPtSigmoid};
    }
  };

  inline OITrackFeatures extractOITrackFeatures(const reco::Track& track,
                                                const reco::TrackCollection& standaloneMuons) {
    OITrackFeatures f;

    const double p = asStored(track.p());
    const double pt = asStored(track.pt());
    const double ptErr = asStored(track.ptError());
    const double eta = asStored(track.eta());
    const double etaErr = asStored(track.etaError());
    const double phi = asStored(track.phi());
    const double phiErr = asStored(track.phiError());
    const double chi2 = asStored(track.chi2());
    const double normalizedChi2 = asStored(track.normalizedChi2());
    const double dszErr = asStored(track.dszError());
    const double dxy = asStored(track.dxy());
    const double dxyErr = asStored(track.dxyError());
    const double dz = asStored(track.dz());
    const double dzErr = asStored(track.dzError());
    const double qoverpErr = asStored(track.qoverpError());
    const double lambdaErr = asStored(track.lambdaError());

    f.logP = logFloor(p);
    f.logPt = logFloor(pt);
    f.logNormalizedChi2 = logFloor(normalizedChi2);
    f.logEtaErr = logFloor(etaErr);
    f.logPhiErr = logFloor(phiErr);
    f.logDszErr = logFloor(dszErr);
    f.logDxyErr = logFloor(dxyErr);
    f.logDzErr = logFloor(dzErr);
    f.logQoverpErr = logFloor(qoverpErr);
    f.logLambdaErr = logFloor(lambdaErr);

    f.eta = eta;
    f.nPixelHits = track.hitPattern().numberOfValidPixelHits();
    f.nTrkLays = track.hitPattern().trackerLayersWithMeasurement();
    const int nFound = track.numberOfValidHits();
    const int nLost = track.numberOfLostHits();
    f.nFoundHits = nFound;
    f.nLostHits = nLost;

    // Derived track-quality features
    f.logImpact3DSq = std::log10(dxy * dxy + dz * dz + kEpsilon);
    const double dxySignificance = dxy / std::max(dxyErr, kEpsilon);
    const double dzSignificance = dz / std::max(dzErr, kEpsilon);
    f.logImpactSignificance =
        std::log10(std::sqrt(dxySignificance * dxySignificance + dzSignificance * dzSignificance) + kEpsilon);
    f.logChi2PerHit = std::log10(chi2 / std::max(nFound, 1) + kEpsilon);
    f.hitEfficiency = static_cast<double>(nFound) / std::max(nFound + nLost, 1);

    // Standalone (L2 muon vertex) matching
    double bestScore = kOINoMatchBestScore;
    for (const auto& muon : standaloneMuons) {
      const double muEta = asStored(muon.eta());
      const double muEtaErr = asStored(muon.etaError());
      const double muPhiErr = asStored(muon.phiError());
      const double muPt = asStored(muon.pt());
      const double muPtErr = asStored(muon.ptError());
      const double muDz = asStored(muon.dz());
      const double muDzErr = asStored(muon.dzError());

      const double chi2Eta =
          (eta - muEta) * (eta - muEta) / (etaErr * etaErr + muEtaErr * muEtaErr + kOIStandaloneChi2Eps);
      const double dPhi = reco::deltaPhi(phi, asStored(muon.phi()));
      const double chi2Phi = dPhi * dPhi / (phiErr * phiErr + muPhiErr * muPhiErr + kOIStandaloneChi2Eps);
      const double chi2Pt = (pt - muPt) * (pt - muPt) / (ptErr * ptErr + muPtErr * muPtErr + kOIStandaloneChi2Eps);
      const double chi2Dz = (dz - muDz) * (dz - muDz) / (dzErr * dzErr + muDzErr * muDzErr + kOIStandaloneChi2Eps);

      const double score = (chi2Eta + chi2Phi + chi2Pt + chi2Dz) / kOIMatchChi2;
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
    const double exponent = std::clamp((pt - kLowPtCut) * 2.0, -20.0, 20.0);
    f.lowPtSigmoid = 1.0 / (1.0 + std::exp(exponent));

    return f;
  }

}  // namespace muonhp

#endif
