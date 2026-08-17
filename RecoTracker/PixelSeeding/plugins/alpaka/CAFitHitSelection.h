#ifndef RecoTracker_PixelSeeding_plugins_alpaka_CAFitHitSelection_h
#define RecoTracker_PixelSeeding_plugins_alpaka_CAFitHitSelection_h

#include <cstdint>

#include "RecoTracker/PixelSeeding/interface/OTHitTag.h"  // caOTHitTag::isOTId (tagged-OT-extra skip)

// Compile-time selection of which hits the BrokenLine fit uses (and how tracks are
// binned by multiplicity). Track FINDING is unaffected -- it stays IT+OT, so the usual
// pixel-track collection is produced with full (hermetic) efficiency; only the FIT (and
// the multiplicity binning that feeds it) is restricted to the selected hit subset.
//
//   All    : use every hit on the track (default)
//   ITonly : use only IT/pixel hits  (drop OT stubs from the fit)
//   OTonly : use only OT/stub hits   (drop IT pixel hits from the fit)
//
// Edit kMode + recompile to switch. Tracks with < 3 selected hits cannot be fit and are skipped.
namespace caFitHitSel {

  enum class Mode { All, ITonly, OTonly };

  constexpr Mode kMode = Mode::All;

  // Whether a hit (given its isStub flag) is used by the fit in the current mode.
  // hasStubs guards the case where the OT collections aren't present (offsetStubs < 0):
  // there, every hit is a pixel hit and we keep everything.
  constexpr inline bool useHit(bool isStub, bool hasStubs) {
    if (kMode == Mode::All || !hasStubs)
      return true;
    return (kMode == Mode::ITonly) ? !isStub : isStub;
  }

  // --- Overlap-module (same-layer) PIXEL hit de-duplication for the BL fit ---------------------
  constexpr float kDedupDsMin = 0.1f;  // cm; transverse merge radius (pixel-pixel only)
  constexpr float kDedupDsMin2 = kDedupDsMin * kDedupDsMin;

  // Single walk used by BOTH the multiplicity counter and the fit's hit selection, so they can never
  // disagree (a mismatch would trip the nSel in [nHitsL,nHitsH] assertions).
  //   k <  0 -> return the number of hits the fit will use (deduped, kMode-filtered).
  //   k >= 0 -> return the position j in [0,nhits) of the k-th kept hit.
  // excludeTaggedOT: when set, a tagged OT-rechit extra id (bit30) is SKIPPED rather than treated as an
  // overflow terminator. A tagged id indexes the raw OT source, not this merged hit SoA, and is attached
  // to the track by the incremental-verify KF (not by the capacity-limited BL fit), so it must NOT be
  // part of this fit-eligibility / multiplicity count nor truncate it. Used only by the extension's
  // verify guard; every other caller leaves it false and sees a tagged id terminate the walk.
  template <typename TupleCont, typename HitsView>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE uint32_t dedupWalk(TupleCont const* __restrict__ foundNtuplets,
                                                    uint32_t it,
                                                    HitsView hh,
                                                    bool hasStubs,
                                                    int k,
                                                    bool excludeTaggedOT = false) {
    auto const* hitId = foundNtuplets->begin(it);
    auto const nhits = foundNtuplets->size(it);
    auto const nTot = hh.metadata().size();
    uint32_t nkept = 0;
    int lastKeptJ = -1;
    for (uint32_t j = 0; j < uint32_t(nhits); ++j) {
      auto const h = hitId[j];
      if (h >= static_cast<unsigned int>(nTot)) {
        if (excludeTaggedOT && caOTHitTag::isOTId(h))
          continue;  // tagged OT extra -> not a merged-fit member; skip, keep walking
        break;       // content buffer overflow guard
      }
      // Stub-ness comes from reco::isStub(view, idx): stubs are the merged-collection hits at index
      // >= hh.offsetStubs().
      if (!useHit(isStub(hh, int32_t(h)), hasStubs))
        continue;
      // Merge only two consecutive PIXEL hits that are near-coincident in the transverse plane.
      if (hasStubs && lastKeptJ >= 0 && !isStub(hh, int32_t(h)) && !isStub(hh, int32_t(hitId[lastKeptJ]))) {
        auto const hp = hitId[lastKeptJ];
        float const dx = float(hh[h].xGlobal()) - float(hh[hp].xGlobal());
        float const dy = float(hh[h].yGlobal()) - float(hh[hp].yGlobal());
        if (dx * dx + dy * dy < kDedupDsMin2)
          continue;  // overlap partner of the previous kept pixel -> merged (keep the earlier hit)
      }
      if (int(nkept) == k)
        return j;  // k-th kept hit
      ++nkept;
      lastKeptJ = int(j);
    }
    return nkept;  // k < 0 (or out of range): the kept count
  }

}  // namespace caFitHitSel

#endif  // RecoTracker_PixelSeeding_plugins_alpaka_CAFitHitSelection_h
