#ifndef RecoMuon_L3TrackFinder_IOTrackSelectorFeatures_h
#define RecoMuon_L3TrackFinder_IOTrackSelectorFeatures_h

// Feature extraction for the muon inside-out (IO) track high-purity selectors
// (pixel tracks and IO/seeds tracks), used by MuonIOTracksForestSelector.
// The feature set is identical to the retired DNN selector it replaces.
//
// The feature set is the 33-feature production set the deployed forests were
// trained with. The struct gives every
// feature a name; toArray() returns them in the canonical training order used
// by the compact-forest/ONNX inputs. Training-side counterpart:
// the build_dataset() feature extraction of the XGBoost trainers (33 features).

#include <array>
#include <cmath>
#include <limits>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

namespace muonhp {

  // Matching/imputation constants - must equal the Python feature extraction.
  inline constexpr float kEpsilon = 1e-6f;        // generic log/division floor
  inline constexpr float kMatchDR2Cut = 0.09f;    // 0.3^2
  inline constexpr float kMatchChi2PtCut = 9.0f;  // 3 sigma
  inline constexpr float kLooseDR2Cut = 0.25f;    // 0.5^2
  inline constexpr float kLooseChi2PtCut = 25.0f;
  inline constexpr float kChi2PtEps = 1e-12f;
  inline constexpr float kDPtNormEps = 1e-9f;
  inline constexpr float kSentinel = 999.0f;
  inline constexpr float kImputeDR2 = 0.1f;
  inline constexpr float kImputeDPtNorm = 1.0f;
  inline constexpr float kImputeChi2Pt = 10.0f;
  inline constexpr float kImputeMatchScore = 0.2f;
  inline constexpr float kImputeSecondDR2 = 1.0f;
  inline constexpr float kLowPtCut = 5.0f;  // GeV (soft sigmoid centre)

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
    float logImpact3DSq = 0.f;   // log10(dxy^2 + dz^2 + eps)
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
              logP, logPt, logEtaErr, logPhiErr, logDzErr, logQoverpErr, logLambdaErr,
              // 7-9 hit content
              nPixelHits, nTrkLays, nFoundHits,
              // 10-17 derived track quality
              logImpact3DSq, logSigmaPtOverPt, logSip2D, logSipZ, logDxyOverPt, logPtErrOverP, logDzOverDxy, absEta,
              // 18-24 stub summary
              nStubs, nStubsEndcap, nStubsBarrel, stubQualMax, stubMaxEtaRegion, stubMaxPhiRegion,
              stubMaxDepthRegion,
              // 25-31 L1 matching
              hasL1Match, logDR2Min, logDPtNorm, logChi2Pt, logMatchingScore, nCompatible, logSecondBestDR2,
              // 32 regime indicator
              lowPtSigmoid};
    }
  };

  inline IOTrackFeatures extractIOTrackFeatures(const reco::Track& track,
                                                const l1t::TrackerMuonCollection& l1TkMuons) {
    IOTrackFeatures f;

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

    // Log-compressed track parameters
    f.logP = std::log10(std::abs(p) + kEpsilon);
    f.logPt = std::log10(std::abs(pt) + kEpsilon);
    f.logEtaErr = std::log10(std::abs(etaErr) + kEpsilon);
    f.logPhiErr = std::log10(std::abs(phiErr) + kEpsilon);
    f.logDzErr = std::log10(std::abs(dzErr) + kEpsilon);
    f.logQoverpErr = std::log10(std::abs(qoverpErr) + kEpsilon);
    f.logLambdaErr = std::log10(std::abs(lambdaErr) + kEpsilon);

    // Hit content
    f.nPixelHits = static_cast<float>(track.hitPattern().numberOfValidPixelHits());
    f.nTrkLays = static_cast<float>(track.hitPattern().trackerLayersWithMeasurement());
    f.nFoundHits = static_cast<float>(track.numberOfValidHits());

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
    float minDR2 = std::numeric_limits<float>::max();
    float matchedL1Pt = -1.0f;
    int bestIndex = -1;
    int nCompatible = 0;

    // Pass 1: best match (chi2Pt-compatible) and loose-window count.
    for (size_t l1Idx = 0; l1Idx != l1TkMuons.size(); ++l1Idx) {
      const auto& l1TkMu = l1TkMuons.at(l1Idx);
      // Propagated muon-system kinematics (matches the training n-tuple).
      const float l1Eta = l1TkMu.phEta();
      const float l1Phi = l1TkMu.phPhi();
      const float l1Pt = l1TkMu.phPt();

      const float ptDiff = pt - l1Pt;
      const float chi2Pt = (ptDiff * ptDiff) / (ptErr * ptErr + kChi2PtEps);
      const float dR2 = reco::deltaR2(eta, phi, l1Eta, l1Phi);

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

    // Pass 2: second-best dR2 strictly greater than the best (mirrors the
    // Python semantics where exact ties are all marked as best and thus
    // excluded from the second-best computation).
    float secondBestDR2 = kSentinel;
    if (bestIndex >= 0) {
      for (size_t l1Idx = 0; l1Idx != l1TkMuons.size(); ++l1Idx) {
        if (static_cast<int>(l1Idx) == bestIndex)
          continue;
        const auto& l1TkMu = l1TkMuons.at(l1Idx);
        const float ptDiff = pt - l1TkMu.phPt();
        const float chi2Pt = (ptDiff * ptDiff) / (ptErr * ptErr + kChi2PtEps);
        if (chi2Pt >= kMatchChi2PtCut)
          continue;
        const float dR2 = reco::deltaR2(eta, phi, l1TkMu.phEta(), l1TkMu.phPhi());
        if (dR2 > minDR2 && dR2 < secondBestDR2)
          secondBestDR2 = dR2;
      }
    }

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

    // Stub summary of the best-matched L1 tracker muon. Best stub = highest
    // quality, ties broken by smallest depthRegion, first-seen wins on full
    // ties (matches the Python ak.firsts on the masked stub arrays).
    if (hasL1Match) {
      const auto& bestL1 = l1TkMuons[bestIndex];
      int nStubsTotal = 0, nStubsEndcap = 0, nStubsBarrel = 0;
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

        if (stubRef->quality() > maxStubQuality ||
            (stubRef->quality() == maxStubQuality && stubRef->depthRegion() < minDepthRegion)) {
          maxStubQuality = stubRef->quality();
          minDepthRegion = stubRef->depthRegion();
          bestStubIndex = static_cast<int>(s);
        }
      }

      f.nStubs = static_cast<float>(nStubsTotal);
      f.nStubsEndcap = static_cast<float>(nStubsEndcap);
      f.nStubsBarrel = static_cast<float>(nStubsBarrel);
      if (bestStubIndex >= 0) {
        const auto bestStub = bestL1.stubs().at(bestStubIndex);
        f.stubQualMax = static_cast<float>(bestStub->quality());
        f.stubMaxEtaRegion = static_cast<float>(bestStub->etaRegion());
        f.stubMaxPhiRegion = static_cast<float>(bestStub->phiRegion());
        f.stubMaxDepthRegion = static_cast<float>(bestStub->depthRegion());
      } else {
        // L1 matched but no usable stubs (Python: maxQual fill_none=0,
        // best-stub regions fill_none=-1).
        f.stubMaxEtaRegion = -1.0f;
        f.stubMaxPhiRegion = -1.0f;
        f.stubMaxDepthRegion = -1.0f;
      }
    } else {
      f.stubMaxEtaRegion = -1.0f;
      f.stubMaxPhiRegion = -1.0f;
      f.stubMaxDepthRegion = -1.0f;
    }

    f.hasL1Match = hasL1Match ? 1.0f : 0.0f;
    if (hasL1Match) {
      f.logDR2Min = std::log10(std::abs(minDR2) + kEpsilon);
      f.logDPtNorm = std::log10(std::abs(dPtNorm) + kEpsilon);
      f.logChi2Pt = std::log10(std::abs(chi2PtBest) + kEpsilon);
      f.logMatchingScore = std::log10(std::abs(matchingScore) + kEpsilon);
    } else {
      f.logDR2Min = std::log10(std::abs(kImputeDR2) + kEpsilon);
      f.logDPtNorm = std::log10(std::abs(kImputeDPtNorm) + kEpsilon);
      f.logChi2Pt = std::log10(std::abs(kImputeChi2Pt) + kEpsilon);
      f.logMatchingScore = std::log10(std::abs(kImputeMatchScore) + kEpsilon);
    }
    f.nCompatible = static_cast<float>(nCompatible);
    const bool hasSecond = secondBestDR2 < (kSentinel - 1.0f);
    f.logSecondBestDR2 = std::log10(std::abs(hasSecond ? secondBestDR2 : kImputeSecondDR2) + kEpsilon);

    // Regime indicator
    const float exponent = std::clamp((pt - kLowPtCut) * 2.0f, -20.0f, 20.0f);
    f.lowPtSigmoid = 1.0f / (1.0f + std::exp(exponent));

    return f;
  }

}  // namespace muonhp

#endif
