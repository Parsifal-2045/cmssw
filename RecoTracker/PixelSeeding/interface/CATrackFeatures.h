#ifndef RecoTracker_PixelSeeding_interface_CATrackFeatures_h
#define RecoTracker_PixelSeeding_interface_CATrackFeatures_h

// Feature vector of the displaced-track classifier (CATrackDNN.h), in the order the model was trained
// with (train_disp_nano.py FEATS). Stage 1 is the in-kernel loose->tight gate on the 12 features;
// Stage 2 is PixelTrackTorchHighPuritySelector, which adds the rzKappaOut extras.
// fill() is shared by three call sites that must produce bit-identical values: Kernel_classifyTracks
// (Stage 1), PixelTrackTorchHighPuritySelector (Stage 2) and the host CA-features table producer.

#include <cmath>
#include <cstdint>

#include "DataFormats/TrackingRecHitSoA/interface/OTRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/StubsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoTracker/PixelSeeding/interface/CircleEq.h"
#include "RecoTracker/PixelSeeding/interface/OTHitTag.h"

namespace caTrackFeatures {

  inline constexpr int kNFeat = 12;

  // Shared inverse-variance stub-curvature kernel: given rg2, bend d=dPhiDr and error
  // s=dPhiDrError, return den = 1 + rg2*d^2 and weight w = den^3 / s^2. rg2 is caller-passed:
  // the xg^2+yg^2 and rGlobal()^2 forms differ at the bit level, so each call site keeps its own
  // to stay ABI-exact with the trained models.
  ALPAKA_FN_HOST_ACC inline void stubDenWeight(float rg2, float d, float s, float &den, float &w) {
    den = 1.f + rg2 * d * d;
    w = den * den * den / (s * s);
  }

  // Feature order == trained ABI (train_disp_nano.py FEATS):
  //   0 fitChi2  1 psFrac  2 r0  3 nPS  4 nh  5 spanZ
  //   6 nStubs   7 nl      8 logChi2Stub  9 kErr  10 dcaEst  11 nBarrel
  //
  // HitIter: forward iterator over merged-hit indices. Returns false on empty/corrupt list.
  //
  // otView: raw OT-rechit SoA for extension-walk hits; a hit id with kOTHitTag set indexes THIS
  // view (otIdx(id)), not hh. Its global position enters the geometry features like a pixel hit;
  // it is "not a stub" (OT SoA has no dPhiDr/stubFlags). Null view with tagged ids present -> false.
  // rGlobal() is stored on the merged SoA but not on the OT SoA -> derived inline.
  //
  // rzKappaOut (Stage-2): when non-null, the SAME hit walk produces out[0]=rzChi2 (reduced chi2 of
  // a straight line z=a+b*r; -1 undefined), out[1]=meanStubKappa, out[2]=leverArm (rMax-r0),
  // out[3]=rMax. Null -> not computed (Stage-1 path unaffected). Caller MUST size out[] for 4 floats.
  template <typename HitIter>
  ALPAKA_FN_HOST_ACC inline bool fill(HitIter hitBegin,
                                      HitIter hitEnd,
                                      ::reco::TrackingRecHitConstView const &hh,
                                      int nHitsTot,
                                      float nLayers,
                                      float chi2,
                                      float *feat,
                                      float *rzKappaOut = nullptr,
                                      ::reco::OTRecHitsConstView const *otView = nullptr) {
    int nh = 0;
    for (auto ph = hitBegin; ph != hitEnd; ++ph)
      ++nh;
    if (nh < 3)
      return false;

    float x0 = 0.f, y0 = 0.f, z0 = 0.f, r0 = 0.f;
    float xm = 0.f, ym = 0.f;
    float xN = 0.f, yN = 0.f, zN = 0.f;
    int nStubs = 0, nPS = 0, nBarrel = 0, nStubK = 0;
    float sumW = 0.f, sumWK = 0.f, sumWK2 = 0.f;
    // Stage-2 extras accumulators, touched only when rzKappaOut != nullptr (dead code in the Stage-1
    // instantiation that passes nullptr, so its register pressure is unchanged).
    float sumKw = 0.f, sumKwk = 0.f;
    float rMaxE = 0.f, r0E = -1.f;  // radial extent (rGlobal-based, matches the host table)
    double Sr = 0, Sz = 0, Srr = 0, Srz = 0, Szz = 0;
    const int iMid = nh / 2;
    int i = 0;
    for (auto ph = hitBegin; ph != hitEnd; ++ph, ++i) {
      const uint32_t h = *ph;
      float xg, yg, zg;
      bool otHit = false;
      if (caOTHitTag::isOTId(h)) {
        // Raw OT-rechit extra: the id indexes the OT SoA, not hh. Without otView -> bail.
        if (!otView)
          return false;
        const uint32_t o = caOTHitTag::otIdx(h);
        // Corruption guard, symmetric to the h >= nHitsTot check below: a container-content
        // overflow can produce garbage words with the tag bit set; never index the OT view out of range.
        if (o >= uint32_t(otView->metadata().size()))
          return false;
        xg = (*otView)[o].xGlobal();
        yg = (*otView)[o].yGlobal();
        zg = (*otView)[o].zGlobal();
        otHit = true;
      } else {
        if (h >= static_cast<uint32_t>(nHitsTot))
          return false;  // content overflow corruption guard
        xg = hh[h].xGlobal();
        yg = hh[h].yGlobal();
        zg = hh[h].zGlobal();
      }
      if (i == 0) {
        x0 = xg, y0 = yg, z0 = zg, r0 = std::sqrt(xg * xg + yg * yg);
      }
      if (i == iMid)
        xm = xg, ym = yg;
      xN = xg, yN = yg, zN = zg;
      // r-z linearity sums (rGlobal()-based, matches the host table; OT SoA derives rGlobal from x/y).
      float rg = 0.f;
      if (rzKappaOut) {
        rg = otHit ? std::sqrt(xg * xg + yg * yg) : hh[h].rGlobal();
        if (r0E < 0.f)
          r0E = rg;  // first hit's radius (host: if (r0<0) r0=rg)
        rMaxE = std::max(rMaxE, rg);
        Sr += rg;
        Sz += zg;
        Srr += double(rg) * rg;
        Srz += double(rg) * zg;
        Szz += double(zg) * zg;
      }
      // OT extras are never stubs (no dPhiDr/stubFlags columns) -> contribute as a pixel hit does.
      if (!otHit && ::reco::isStub(hh, h)) {
        ++nStubs;
        const auto flags = hh[h].stubFlags();
        if (::reco::StubFlags::isPS(flags))
          ++nPS;
        if (::reco::StubFlags::isBarrel(flags))
          ++nBarrel;
        const float s = hh[h].dPhiDrError();
        if (s > 0.f) {
          const float d = hh[h].dPhiDr();
          float den, w;
          stubDenWeight(xg * xg + yg * yg, d, s, den, w);
          const float k = d / std::sqrt(den);
          sumW += w;
          sumWK += w * k;
          sumWK2 += w * k * k;
          ++nStubK;
          if (rzKappaOut) {
            // meanStubKappa weight uses rGlobal()^2 (NOT xg^2+yg^2), matching the host table bit for bit.
            float denR, wR;
            stubDenWeight(rg * rg, d, s, denR, wR);
            sumKw += wR;
            sumKwk += wR * (d / std::sqrt(denR));
          }
        }
      }
    }

    const float kErr = (sumW > 0.f) ? 1.f / std::sqrt(sumW) : -1.f;
    const float chi2Stub = (nStubK >= 3 && sumW > 0.f) ? (sumWK2 - sumWK * sumWK / sumW) / float(nStubK - 1) : -1.f;
    CircleEq<float> eq(x0, y0, xm, ym, xN, yN);
    const float curv3 = eq.curvature();
    const float dcaEst = (std::abs(curv3) > 0.f) ? std::abs(eq.dca0()) / std::abs(curv3) : 0.f;

    // fitChi2: sanitized (non-finite/negative/huge raw chi2 -> worst cap) and log1p-compressed
    // (mirrors feat[8]); the trainer applies the identical transform so mean/scale stay consistent.
    float fitChi2Safe = chi2;
    if (!(fitChi2Safe >= 0.f) || fitChi2Safe > 1.0e4f)  // catches NaN, negative, +inf, huge
      fitChi2Safe = 1.0e4f;
    feat[0] = std::log1p(fitChi2Safe);                  // fitChi2 (sanitized, log)
    feat[1] = float(nPS) / float(std::max(nStubs, 1));  // psFrac
    feat[2] = r0;
    feat[3] = float(nPS);
    feat[4] = float(nh);
    feat[5] = std::abs(zN - z0);  // spanZ
    feat[6] = float(nStubs);
    feat[7] = nLayers;                              // nl
    feat[8] = std::log1p(std::max(chi2Stub, 0.f));  // logChi2Stub
    feat[9] = kErr;
    feat[10] = dcaEst;
    feat[11] = float(nBarrel);
    // Stage-2 extras, finalized from the single walk above (same arithmetic as the host producer).
    if (rzKappaOut) {
      rzKappaOut[1] = (sumKw > 0.f) ? sumKwk / sumKw : 0.f;  // meanStubKappa
      float rz = -1.f;                                       // rzChi2 (-1 = undefined)
      // nh hits were all validated above (fill returns false on a bad index before reaching here).
      if (nh >= 3) {
        const double D = nh * Srr - Sr * Sr;
        if (std::abs(D) > 0.0) {
          const double b = (nh * Srz - Sr * Sz) / D, a = (Sz - b * Sr) / nh;
          rz = float(std::max(0.0, (Szz - a * Sz - b * Srz)) / std::max(1, nh - 2));
        }
      }
      rzKappaOut[0] = rz;
      // Radial-extent pair (host CATrackFeaturesTableProducer: leverArm = rMax - r0, rMax = rMax).
      rzKappaOut[2] = rMaxE - (r0E < 0.f ? 0.f : r0E);  // leverArm
      rzKappaOut[3] = rMaxE;                            // rMax
    }
    return true;
  }

}  // namespace caTrackFeatures

#endif  // RecoTracker_PixelSeeding_interface_CATrackFeatures_h
