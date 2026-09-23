#ifndef RecoMuon_L3TrackFinder_IOTrackSelectorFeatures_h
#define RecoMuon_L3TrackFinder_IOTrackSelectorFeatures_h

// Feature extraction for the muon inside-out (IO) track high-purity selectors
// (pixel tracks and IO/seeds tracks), used by MuonIOTracksForestSelector.
//
// The feature set is the 33-feature production set the deployed forests were
// trained with. The struct gives every feature a name; toArray() returns them
// in the canonical training order used by the compact-forest/ONNX inputs.
// Training-side counterpart: build_dataset() + io_production_features() in
// muonHighPurityTrackSelection/production/pixel_features.py (the training
// pipeline asserts that its kept features follow this order).
//
// Numeric convention (identical in the training extraction): every raw input
// is rounded to float, as the training n-tuple stores it; all arithmetic is
// done in double precision; each feature is rounded to float once. The
// forests split on exact float values (constant imputed features sit exactly
// on split points), so float-vs-double differences between training and
// inference are not harmless: they must not exist.

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace muonhp {

  // Matching/imputation constants - must equal the Python feature extraction.
  inline constexpr double kEpsilon = 1e-6;           // generic log/division floor
  inline constexpr double kMatchDR2Cut = 0.3 * 0.3;  // as 0.3**2 in the training code
  inline constexpr double kMatchChi2PtCut = 9.0;     // 3 sigma
  inline constexpr double kLooseDR2Cut = 0.5 * 0.5;
  inline constexpr double kLooseChi2PtCut = 25.0;
  inline constexpr double kChi2PtEps = 1e-12;
  inline constexpr double kDPtNormEps = 1e-9;
  inline constexpr double kSentinel = 999.0;
  inline constexpr double kImputeDR2 = 0.1;
  inline constexpr double kImputeDPtNorm = 1.0;
  inline constexpr double kImputeChi2Pt = 10.0;
  inline constexpr double kImputeMatchScore = 0.2;
  inline constexpr double kImputeSecondDR2 = 1.0;
  inline constexpr double kLowPtCut = 5.0;  // GeV (soft sigmoid centre)

  // Raw input as stored in the training n-tuple (float), promoted to double.
  inline double asStored(double x) { return static_cast<float>(x); }
  // log10(|x| + eps), the log compression of the training extraction.
  inline double logFloor(double x) { return std::log10(std::abs(x) + kEpsilon); }

  // The 33-feature IO production set, in canonical (training) order.
  struct IOTrackFeatures {
    // Log-compressed track parameters
    float logP = 0.f;
    float logPt = 0.f;
    float logEtaErr = 0.f;
    float logPhiErr = 0.f;
    float logDzErr = 0.f;
    float logQoverpErr = 0.f;
    float logLambdaErr = 0.f;
    // Hit content
    float nPixelHits = 0.f;
    float nTrkLays = 0.f;
    float nFoundHits = 0.f;
    // Derived track-quality features
    float logImpact3DSq = 0.f;  // log10(dxy^2 + dz^2 + eps)
    float logSigmaPtOverPt = 0.f;
    float logSip2D = 0.f;
    float logSipZ = 0.f;
    float logDxyOverPt = 0.f;
    float logPtErrOverP = 0.f;
    float logDzOverDxy = 0.f;
    float absEta = 0.f;
    // Stub summary of the best-matched L1 tracker muon
    float nStubs = 0.f;
    float nStubsEndcap = 0.f;
    float nStubsBarrel = 0.f;
    float stubQualMax = 0.f;
    float stubMaxEtaRegion = 0.f;
    float stubMaxPhiRegion = 0.f;
    float stubMaxDepthRegion = 0.f;
    // L1 tracker-muon matching
    float hasL1Match = 0.f;
    float logDR2Min = 0.f;
    float logDPtNorm = 0.f;
    float logChi2Pt = 0.f;
    float logMatchingScore = 0.f;
    float nCompatible = 0.f;
    float logSecondBestDR2 = 0.f;
    // Regime indicator (soft sigmoid around kLowPtCut)
    float lowPtSigmoid = 0.f;

    static constexpr size_t kSize = 33;

    // Canonical ordering consumed by the deployed models. The order below is
    // the training ABI and must not change without retraining.
    std::array<float, kSize> toArray() const {
      return {// 0-6 log track parameters
              logP,
              logPt,
              logEtaErr,
              logPhiErr,
              logDzErr,
              logQoverpErr,
              logLambdaErr,
              // 7-9 hit content
              nPixelHits,
              nTrkLays,
              nFoundHits,
              // 10-17 derived track quality
              logImpact3DSq,
              logSigmaPtOverPt,
              logSip2D,
              logSipZ,
              logDxyOverPt,
              logPtErrOverP,
              logDzOverDxy,
              absEta,
              // 18-24 stub summary
              nStubs,
              nStubsEndcap,
              nStubsBarrel,
              stubQualMax,
              stubMaxEtaRegion,
              stubMaxPhiRegion,
              stubMaxDepthRegion,
              // 25-31 L1 matching
              hasL1Match,
              logDR2Min,
              logDPtNorm,
              logChi2Pt,
              logMatchingScore,
              nCompatible,
              logSecondBestDR2,
              // 32 regime indicator
              lowPtSigmoid};
    }
  };

  inline IOTrackFeatures extractIOTrackFeatures(const reco::Track& track, const l1t::TrackerMuonCollection& l1TkMuons) {
    IOTrackFeatures f;

    // Raw track quantities
    const double p = asStored(track.p());
    const double pt = asStored(track.pt());
    const double ptErr = asStored(track.ptError());
    const double eta = asStored(track.eta());
    const double etaErr = asStored(track.etaError());
    const double phi = asStored(track.phi());
    const double phiErr = asStored(track.phiError());
    const double dxy = asStored(track.dxy());
    const double dxyErr = asStored(track.dxyError());
    const double dz = asStored(track.dz());
    const double dzErr = asStored(track.dzError());
    const double qoverpErr = asStored(track.qoverpError());
    const double lambdaErr = asStored(track.lambdaError());

    // Log-compressed track parameters
    f.logP = logFloor(p);
    f.logPt = logFloor(pt);
    f.logEtaErr = logFloor(etaErr);
    f.logPhiErr = logFloor(phiErr);
    f.logDzErr = logFloor(dzErr);
    f.logQoverpErr = logFloor(qoverpErr);
    f.logLambdaErr = logFloor(lambdaErr);

    // Hit content
    f.nPixelHits = track.hitPattern().numberOfValidPixelHits();
    f.nTrkLays = track.hitPattern().trackerLayersWithMeasurement();
    f.nFoundHits = track.numberOfValidHits();

    // Derived track-quality features
    f.logImpact3DSq = std::log10(dxy * dxy + dz * dz + kEpsilon);
    f.logSigmaPtOverPt = std::log10(ptErr / std::max(pt, kEpsilon) + kEpsilon);
    f.logSip2D = std::log10(std::abs(dxy) / std::max(dxyErr, kEpsilon) + kEpsilon);
    f.logSipZ = std::log10(std::abs(dz) / std::max(dzErr, kEpsilon) + kEpsilon);
    f.logDxyOverPt = std::log10(std::abs(dxy) / std::max(pt, kEpsilon) + kEpsilon);
    f.logPtErrOverP = std::log10(ptErr / std::max(p, kEpsilon) + kEpsilon);
    f.logDzOverDxy = std::log10(std::abs(dz) / (std::abs(dxy) + kEpsilon) + kEpsilon);
    f.absEta = std::abs(eta);

    // ------------------------------------------------------------------
    // L1 tracker-muon matching and stub summary
    // ------------------------------------------------------------------
    double minDR2 = std::numeric_limits<double>::max();
    double matchedL1Pt = -1.0;
    int bestIndex = -1;
    int nCompatible = 0;

    // Pass 1: best match (chi2Pt-compatible) and loose-window count.
    for (size_t l1Idx = 0; l1Idx != l1TkMuons.size(); ++l1Idx) {
      const auto& l1TkMu = l1TkMuons[l1Idx];
      // Propagated muon-system kinematics (L1TkMu_pt/eta/phi of the n-tuple).
      const double l1Pt = asStored(l1TkMu.phPt());
      const double ptDiff = pt - l1Pt;
      const double chi2Pt = (ptDiff * ptDiff) / (ptErr * ptErr + kChi2PtEps);
      const double dR2 = reco::deltaR2(eta, phi, asStored(l1TkMu.phEta()), asStored(l1TkMu.phPhi()));

      if (dR2 < kLooseDR2Cut && chi2Pt < kLooseChi2PtCut)
        ++nCompatible;

      if (chi2Pt >= kMatchChi2PtCut)
        continue;

      if (dR2 < minDR2) {
        minDR2 = dR2;
        matchedL1Pt = l1Pt;
        bestIndex = static_cast<int>(l1Idx);
      }
    }

    // Pass 2: second-best dR2 strictly greater than the best (exact ties are
    // all marked as best in the training code, hence excluded).
    double secondBestDR2 = kSentinel;
    if (bestIndex >= 0) {
      for (size_t l1Idx = 0; l1Idx != l1TkMuons.size(); ++l1Idx) {
        if (static_cast<int>(l1Idx) == bestIndex)
          continue;
        const auto& l1TkMu = l1TkMuons[l1Idx];
        const double ptDiff = pt - asStored(l1TkMu.phPt());
        const double chi2Pt = (ptDiff * ptDiff) / (ptErr * ptErr + kChi2PtEps);
        if (chi2Pt >= kMatchChi2PtCut)
          continue;
        const double dR2 = reco::deltaR2(eta, phi, asStored(l1TkMu.phEta()), asStored(l1TkMu.phPhi()));
        if (dR2 > minDR2 && dR2 < secondBestDR2)
          secondBestDR2 = dR2;
      }
    }

    double dPtNorm = kImputeDPtNorm;
    double chi2PtBest = kImputeChi2Pt;
    double matchingScore = kImputeMatchScore;
    if (bestIndex >= 0) {
      dPtNorm = std::abs(pt - matchedL1Pt) / (matchedL1Pt + kDPtNormEps);
      const double ptDiffBest = pt - matchedL1Pt;
      chi2PtBest = (ptDiffBest * ptDiffBest) / (ptErr * ptErr + kChi2PtEps);
      matchingScore = minDR2 * (1.0 + dPtNorm);
    }

    const bool hasL1Match = (minDR2 < kMatchDR2Cut) && (bestIndex >= 0);

    // Stub summary of the best-matched L1 tracker muon. Best stub = highest
    // quality, ties broken by smallest depthRegion, first-seen wins on full
    // ties (the training code takes ak.firsts of the masked stub arrays).
    f.stubMaxEtaRegion = -1.0f;
    f.stubMaxPhiRegion = -1.0f;
    f.stubMaxDepthRegion = -1.0f;
    if (hasL1Match) {
      const auto& bestL1 = l1TkMuons[bestIndex];
      int nStubsTotal = 0, nStubsEndcap = 0, nStubsBarrel = 0;
      int maxStubQuality = 0;
      int minDepthRegion = std::numeric_limits<int>::max();
      int bestStubIndex = -1;

      for (size_t s = 0; s != bestL1.stubs().size(); ++s) {
        const auto stubRef = bestL1.stubs()[s];
        if (stubRef.isNull())
          continue;
        ++nStubsTotal;
        if (stubRef->type() == 0)
          ++nStubsEndcap;
        else if (stubRef->type() == 1)
          ++nStubsBarrel;

        if (stubRef->quality() > maxStubQuality ||
            (stubRef->quality() == maxStubQuality && stubRef->depthRegion() < minDepthRegion)) {
          maxStubQuality = stubRef->quality();
          minDepthRegion = stubRef->depthRegion();
          bestStubIndex = static_cast<int>(s);
        }
      }

      f.nStubs = nStubsTotal;
      f.nStubsEndcap = nStubsEndcap;
      f.nStubsBarrel = nStubsBarrel;
      // L1 matched but no usable stubs: training fill values (maxQual 0,
      // best-stub regions -1) are kept.
      if (bestStubIndex >= 0) {
        const auto bestStub = bestL1.stubs()[bestStubIndex];
        f.stubQualMax = bestStub->quality();
        f.stubMaxEtaRegion = bestStub->etaRegion();
        f.stubMaxPhiRegion = bestStub->phiRegion();
        f.stubMaxDepthRegion = bestStub->depthRegion();
      }
    }

    f.hasL1Match = hasL1Match ? 1.0f : 0.0f;
    f.logDR2Min = logFloor(hasL1Match ? minDR2 : kImputeDR2);
    f.logDPtNorm = logFloor(hasL1Match ? dPtNorm : kImputeDPtNorm);
    f.logChi2Pt = logFloor(hasL1Match ? chi2PtBest : kImputeChi2Pt);
    f.logMatchingScore = logFloor(hasL1Match ? matchingScore : kImputeMatchScore);
    f.nCompatible = nCompatible;
    const bool hasSecond = secondBestDR2 < (kSentinel - 1.0);
    f.logSecondBestDR2 = logFloor(hasSecond ? secondBestDR2 : kImputeSecondDR2);

    // Regime indicator
    const double exponent = std::clamp((pt - kLowPtCut) * 2.0, -20.0, 20.0);
    f.lowPtSigmoid = 1.0 / (1.0 + std::exp(exponent));

    return f;
  }

}  // namespace muonhp

#endif
