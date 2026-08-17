#ifndef RecoTracker_PixelSeeding_plugins_alpaka_OTStubFormationVectorHitStyleKernels_h
#define RecoTracker_PixelSeeding_plugins_alpaka_OTStubFormationVectorHitStyleKernels_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include "DataFormats/Math/interface/approx_atan2.h"
#include "DataFormats/TrackingRecHitSoA/interface/OTRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/StubsSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoTracker/PixelSeeding/interface/StackedModuleGeometrySoA.h"

// Opt-in per-module stub-formation diagnostics (device printf). When defined, the counting kernel
// runs a sequential single-lane path and DiagnosticSummaryKernel becomes available.
// #define STUB_DIAGNOSTIC_COUNTERS

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace otStubFormationVectorHitStyle {
    // CA module offset parameters, filled by the producer
    struct CAModuleOffsets {
      static constexpr uint32_t nPixelModules = 4000;
      uint32_t barrelStartGeom;
      uint32_t backwardStartGeom;
      uint32_t forwardStartGeom;
      uint32_t nBarrelModules;
      uint32_t nBackwardModules;
      uint32_t nForwardModules;
    };

    // Parallax correction along the lower sensor's local-x axis (VectorHit convention). The
    // origin->hit vector is projected on the stack frame (zLoc = lower->upper stack direction,
    // xLoc = lower local-x made orthogonal to zLoc) and the local slope lVx / lVz is multiplied by
    // the sensor separation. Returns 0 for degenerate inputs.
    template <typename TAcc>
    ALPAKA_FN_ACC ALPAKA_FN_INLINE float computeParallaxCorrection(
        TAcc const& acc,
        float gx_lower,
        float gy_lower,
        float gz_lower,
        // unit vector (global) from the lower to the upper sensor centre
        float globalLowUpNormX,
        float globalLowUpNormY,
        float globalLowUpNormZ,
        // direction of the lower sensor's local-x axis in global coordinates (unit)
        float localXInGlobalX,
        float localXInGlobalY,
        float localXInGlobalZ,
        // sensor separation, in the same length unit as the coordinates
        float separation,
        bool debug = false,
        uint32_t moduleId = 0) {
      constexpr float eps = 1e-6f;

      // gV: vector from the origin to the lower hit (global). The hit coordinates are already
      // beam-spot corrected by the converter, so the origin is the beam spot.
      float gVx = gx_lower;
      float gVy = gy_lower;
      float gVz = gz_lower;

      // zLoc direction: along the stack from LOWER -> UPPER (global, assumed normalized already)
      float zLocX = globalLowUpNormX;
      float zLocY = globalLowUpNormY;
      float zLocZ = globalLowUpNormZ;

      // Safety: ensure zLoc is valid (normalized-ish and non-zero)
      float zLocN2 = zLocX * zLocX + zLocY * zLocY + zLocZ * zLocZ;
      if (zLocN2 < eps)
        return 0.0f;

      // If not perfectly normalized, normalize it
      if (alpaka::math::abs(acc, zLocN2 - 1.0f) > 1e-3f) {
        float invN = 1.0f / alpaka::math::sqrt(acc, zLocN2);
        zLocX *= invN;
        zLocY *= invN;
        zLocZ *= invN;
      }

      // xAxis (global) is the direction of LOWER local-x, but it may not be exactly orthogonal to zLoc.
      // Project it into the plane orthogonal to zLoc to get a proper in-plane x direction.
      float x0X = localXInGlobalX;
      float x0Y = localXInGlobalY;
      float x0Z = localXInGlobalZ;

      float x0N2 = x0X * x0X + x0Y * x0Y + x0Z * x0Z;
      if (x0N2 < eps)
        return 0.0f;

      if (alpaka::math::abs(acc, x0N2 - 1.0f) > 1e-3f) {
        float invN = 1.0f / alpaka::math::sqrt(acc, x0N2);
        x0X *= invN;
        x0Y *= invN;
        x0Z *= invN;
      }

      // xPerp = x0 - (x0*zLoc) zLoc
      float x0DotZ = x0X * zLocX + x0Y * zLocY + x0Z * zLocZ;
      float xPerpX = x0X - x0DotZ * zLocX;
      float xPerpY = x0Y - x0DotZ * zLocY;
      float xPerpZ = x0Z - x0DotZ * zLocZ;

      float xPerpN2 = xPerpX * xPerpX + xPerpY * xPerpY + xPerpZ * xPerpZ;
      if (xPerpN2 < eps) {
        // Degenerate: provided localX is (nearly) parallel to stack direction
        if (debug) {
          printf("    [ParallaxDebug Module %u] WARNING: xPerp too small (localX || stack), returning 0\n", moduleId);
        }
        return 0.0f;
      }

      float invXPerpN = 1.0f / alpaka::math::sqrt(acc, xPerpN2);
      float xLocX = xPerpX * invXPerpN;
      float xLocY = xPerpY * invXPerpN;
      float xLocZ = xPerpZ * invXPerpN;

      // Local components needed for VectorHit-style correction:
      //   lVx = gV * xLoc
      //   lVz = gV * zLoc
      float lVx = gVx * xLocX + gVy * xLocY + gVz * xLocZ;
      float lVz = gVx * zLocX + gVy * zLocY + gVz * zLocZ;

      if (debug) {
        printf("    [ParallaxDebug Module %u] Inputs: gPos=(%.4f, %.4f, %.4f) sep=%.4f\n",
               moduleId,
               gx_lower,
               gy_lower,
               gz_lower,
               separation);
        printf("    [ParallaxDebug Module %u] zLoc(L->U)=(%.6f, %.6f, %.6f)\n", moduleId, zLocX, zLocY, zLocZ);
        printf(
            "    [ParallaxDebug Module %u] x0(localX in global)=(%.6f, %.6f, %.6f)  xLoc(projected)=(%.6f, %.6f, "
            "%.6f)\n",
            moduleId,
            localXInGlobalX,
            localXInGlobalY,
            localXInGlobalZ,
            xLocX,
            xLocY,
            xLocZ);
        printf("    [ParallaxDebug Module %u] Local comps: lVx=%.6f lVz=%.6f\n", moduleId, lVx, lVz);
      }

      // Protect against near-parallel-to-plane direction (would blow up slope)
      if (alpaka::math::abs(acc, lVz) < eps) {
        if (debug) {
          printf("    [ParallaxDebug Module %u] WARNING: |lVz| too small, returning 0\n", moduleId);
        }
        return 0.0f;
      }

      // Normalized slope and correction
      float slope_x = lVx / lVz;
      float parallCorr = slope_x * separation;

      if (debug) {
        printf("    [ParallaxDebug Module %u] slope_x=%.6f parallCorr=%.6f\n", moduleId, slope_x, parallCorr);
      }

      if (!alpaka::math::isfinite(acc, parallCorr)) {
        if (debug) {
          printf("    [ParallaxDebug Module %u] WARNING: parallCorr not finite, returning 0\n", moduleId);
        }
        return 0.0f;
      }

      return parallCorr;  // along LOWER-sensor local_x, transporting LOWER->UPPER along globalLowUpNorm
    }

    // Counting kernel: Count stubs per module using VectorHits cuts
    // Work division: Acc2D with Y indexing modules and X indexing warp lanes
    // (warpSize threads per module). All lanes of a warp cooperate on the
    // (iLower, iUpper) pair loop of a single module.
    class CountStubsVectorHitStyleKernel {
    public:
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::reco::OTRecHitsConstView hits,
                                    ::reco::OTHitModuleConstView moduleView,
                                    ::reco::StackedModuleGeometryConstView geometry,
                                    float const* maxWidthBarrelFlat,
                                    float const* maxWidthBarrelTilted,
                                    float const* maxWidthEndcap,
                                    int32_t const* barrelFlatMaxCSDiff,
                                    int32_t const* barrelTiltedMaxCSDiff,
                                    int32_t const* endcapMaxCSDiff,
                                    int32_t const* barrelFlatMaxCS,
                                    int32_t const* barrelTiltedMaxCS,
                                    int32_t const* endcapMaxCS,
                                    int32_t const* barrelFlatMaxCSSum,
                                    int32_t const* barrelTiltedMaxCSSum,
                                    int32_t const* endcapMaxCSSum,
                                    uint32_t* stubCounts,
                                    uint32_t nModules) const {
        // Mark as maybe_unused to suppress warning when diagnostic counters are enabled (single thread path)
        [[maybe_unused]] const int32_t warpSize = alpaka::warp::getSize(acc);
        const int32_t laneId = static_cast<int32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[1u]);

        // One warp per module
        for (uint32_t iModule : cms::alpakatools::uniform_elements_y(acc, nModules)) {
          // Get hit ranges (uniform across the warp)
          uint32_t hitStart = moduleView[iModule].moduleStart();
          uint32_t hitEnd = moduleView[iModule + 1].moduleStart();
          uint32_t upperStart = moduleView[iModule].upperSensorStart();

          uint32_t nLower = upperStart - hitStart;
          uint32_t nUpper = hitEnd - upperStart;

          // Handle modules with no upper hits (P-hit recovery for PS modules)
          if (nUpper == 0) {
            // For PS modules: recover P-hits that have no matching S-hits.
            // These become PHitOnly stub entries (no bend measurement).
            uint8_t moduleType = geometry[iModule].moduleType();
            bool isPS = (moduleType == 0 || moduleType == 1);
            if (laneId == 0) {
              if (isPS && nLower > 0) {
                // Count each P-hit as a PHitOnly stub entry
                stubCounts[iModule + 1] = nLower;
#ifdef STUB_DIAGNOSTIC_COUNTERS
                printf(
                    "STUB_DIAG_VH module=%u type=PS totalPHits=%u pHitOnlyEntries=%u "
                    "reason=noUpperHit\n",
                    iModule,
                    nLower,
                    nLower);
#endif
              } else {
                // SS modules or no lower hits: no stubs
                stubCounts[iModule + 1] = 0;
              }
            }
            continue;
          }

          // Skip if no lower hits
          if (nLower == 0) {
            if (laneId == 0)
              stubCounts[iModule + 1] = 0;
            continue;
          }

          // Get geometry parameters (uniform across the warp)
          float separation = geometry[iModule].sensorSeparation();  // in mm
          bool isBarrel = geometry[iModule].isBarrel();
          bool isFlat = geometry[iModule].isFlat();
          [[maybe_unused]] bool isFlipped = geometry[iModule].isFlipped();
          uint8_t layer = geometry[iModule].layer();
          float globalLowUpNormX = geometry[iModule].globalLowUpNormX();
          float globalLowUpNormY = geometry[iModule].globalLowUpNormY();
          float globalLowUpNormZ = geometry[iModule].globalLowUpNormZ();
          float localXInGlobalX = geometry[iModule].localXInGlobalX();
          float localXInGlobalY = geometry[iModule].localXInGlobalY();
          float localXInGlobalZ = geometry[iModule].localXInGlobalZ();

          // Get layer-dependent cut (separate cuts for flat and tilted barrel modules)
          float cut =
              isBarrel ? (isFlat ? maxWidthBarrelFlat[layer] : maxWidthBarrelTilted[layer]) : maxWidthEndcap[layer];

          // Per-layer cluster size cuts
          int32_t maxCSDiff =
              isBarrel ? (isFlat ? barrelFlatMaxCSDiff[layer] : barrelTiltedMaxCSDiff[layer]) : endcapMaxCSDiff[layer];
          int32_t maxCS = isBarrel ? (isFlat ? barrelFlatMaxCS[layer] : barrelTiltedMaxCS[layer]) : endcapMaxCS[layer];
          int32_t maxCSSum =
              isBarrel ? (isFlat ? barrelFlatMaxCSSum[layer] : barrelTiltedMaxCSSum[layer]) : endcapMaxCSSum[layer];

          // Topological lower->upper normal (flipped sign for flipped modules; see note in
          // formation kernel below)
          float flipSign = isFlipped ? -1.0f : 1.0f;
          float topoLowUpNormX = flipSign * globalLowUpNormX;
          float topoLowUpNormY = flipSign * globalLowUpNormY;
          float topoLowUpNormZ = flipSign * globalLowUpNormZ;
          float separation_cm = separation * 0.1f;

#ifdef STUB_DIAGNOSTIC_COUNTERS
          // Diagnostic counters require sequential per-lower-hit accounting:
          // run the original 1D algorithm on lane 0 only.
          if (laneId == 0) {
            uint8_t moduleType = geometry[iModule].moduleType();
            bool isPS = (moduleType == 0 || moduleType == 1);
            uint32_t pHitsWithStub = 0;
            uint32_t pHitsLostWidthCut = 0;
            uint32_t maxStubsPerPHit = 0;
            uint32_t count = 0;

            for (uint32_t iLower = hitStart; iLower < upperStart; ++iLower) {
              auto const& lowerHit = hits[iLower];
              uint32_t stubsForThisLowerHit = 0;
              for (uint32_t iUpper = upperStart; iUpper < hitEnd; ++iUpper) {
                auto const& upperHit = hits[iUpper];
                if (lowerHit.yLocal() * upperHit.yLocal() < 0.f)
                  continue;
                if (alpaka::math::abs(acc, int32_t(lowerHit.clusterSize()) - int32_t(upperHit.clusterSize())) >
                    maxCSDiff)
                  continue;
                if (int32_t(lowerHit.clusterSize()) > maxCS || int32_t(upperHit.clusterSize()) > maxCS)
                  continue;
                if (int32_t(lowerHit.clusterSize()) + int32_t(upperHit.clusterSize()) > maxCSSum)
                  continue;
                float lx_lower = lowerHit.xLocal();
                float lx_upper = upperHit.xLocal();
                float pC = computeParallaxCorrection(acc,
                                                     lowerHit.xGlobal(),
                                                     lowerHit.yGlobal(),
                                                     lowerHit.zGlobal(),
                                                     topoLowUpNormX,
                                                     topoLowUpNormY,
                                                     topoLowUpNormZ,
                                                     localXInGlobalX,
                                                     localXInGlobalY,
                                                     localXInGlobalZ,
                                                     separation_cm);
                float lpos_lower_corr = 0.0f;
                float lpos_upper_corr = 0.0f;
                if (lx_upper > lx_lower) {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
                  } else {
                    lpos_lower_corr = lx_lower + alpaka::math::abs(acc, pC);
                    lpos_upper_corr = lx_upper;
                  }
                } else if (lx_upper < lx_lower) {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower - alpaka::math::abs(acc, pC);
                    lpos_upper_corr = lx_upper;
                  } else {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
                  }
                } else {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
                  } else {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
                  }
                }
                float width = lpos_lower_corr - lpos_upper_corr;
                if (alpaka::math::abs(acc, width) < cut) {
                  ++count;
                  ++stubsForThisLowerHit;
                }
              }
              if (stubsForThisLowerHit > 0) {
                ++pHitsWithStub;
                if (stubsForThisLowerHit > maxStubsPerPHit)
                  maxStubsPerPHit = stubsForThisLowerHit;
              } else {
                ++pHitsLostWidthCut;
              }
            }

            stubCounts[iModule + 1] = count;
            const char* modType = isPS ? "PS" : "SS";
            printf(
                "STUB_DIAG_VH module=%u type=%s totalPHits=%u stubsFormed=%u pHitsWithStub=%u pHitsLost=%u "
                "lostReason=widthCut maxStubsPerPHit=%u\n",
                iModule,
                modType,
                nLower,
                count,
                pHitsWithStub,
                pHitsLostWidthCut,
                maxStubsPerPHit);
          }
          continue;
#else
          // Warp-cooperative count: lanes stride across the flattened pair index
          // p = (iLower - hitStart) * nUpper + (iUpper - upperStart).
          // Output ordering is irrelevant here -- we only need the total per module.
          uint32_t nPairs = nLower * nUpper;
          int32_t localCount = 0;
          for (uint32_t p = static_cast<uint32_t>(laneId); p < nPairs; p += warpSize) {
            uint32_t iLower = hitStart + p / nUpper;
            uint32_t iUpper = upperStart + p % nUpper;

            auto const& lowerHit = hits[iLower];
            auto const& upperHit = hits[iUpper];

            // Same-sign local-y cut
            if (lowerHit.yLocal() * upperHit.yLocal() < 0.f)
              continue;
            // Cluster size compatibility cuts (per layer)
            if (alpaka::math::abs(acc, int32_t(lowerHit.clusterSize()) - int32_t(upperHit.clusterSize())) > maxCSDiff)
              continue;
            if (int32_t(lowerHit.clusterSize()) > maxCS || int32_t(upperHit.clusterSize()) > maxCS)
              continue;
            if (int32_t(lowerHit.clusterSize()) + int32_t(upperHit.clusterSize()) > maxCSSum)
              continue;

            float lx_lower = lowerHit.xLocal();
            float lx_upper = upperHit.xLocal();

            float pC = computeParallaxCorrection(acc,
                                                 lowerHit.xGlobal(),
                                                 lowerHit.yGlobal(),
                                                 lowerHit.zGlobal(),
                                                 topoLowUpNormX,
                                                 topoLowUpNormY,
                                                 topoLowUpNormZ,
                                                 localXInGlobalX,
                                                 localXInGlobalY,
                                                 localXInGlobalZ,
                                                 separation_cm);

            float lpos_lower_corr = 0.0f;
            float lpos_upper_corr = 0.0f;
            if (lx_upper > lx_lower) {
              if (lx_upper > 0) {
                lpos_lower_corr = lx_lower;
                lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
              } else {
                lpos_lower_corr = lx_lower + alpaka::math::abs(acc, pC);
                lpos_upper_corr = lx_upper;
              }
            } else if (lx_upper < lx_lower) {
              if (lx_upper > 0) {
                lpos_lower_corr = lx_lower - alpaka::math::abs(acc, pC);
                lpos_upper_corr = lx_upper;
              } else {
                lpos_lower_corr = lx_lower;
                lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
              }
            } else {
              if (lx_upper > 0) {
                lpos_lower_corr = lx_lower;
                lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
              } else {
                lpos_lower_corr = lx_lower;
                lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
              }
            }

            float width = lpos_lower_corr - lpos_upper_corr;
            if (alpaka::math::abs(acc, width) < cut)
              ++localCount;
          }

          // Warp-reduce localCount across lanes (xor-butterfly).
          // All lanes must be active for the shuffles to work, so no early return above.
          for (int32_t off = 1; off < warpSize; off <<= 1) {
            localCount += alpaka::warp::shfl_xor(acc, localCount, off);
          }
          if (laneId == 0) {
            stubCounts[iModule + 1] = static_cast<uint32_t>(localCount);
          }
#endif  // STUB_DIAGNOSTIC_COUNTERS
        }
      }
    };

    // Formation kernel for VectorHits style.
    // Work division: Acc2D with Y indexing modules and X indexing warp lanes
    // (warpSize threads per module). All lanes of a warp cooperate on the
    // (iLower, iUpper) pair loop of a single module. Lanes consume pairs in their
    // natural lexicographic order; a warp exclusive prefix scan assigns each
    // passing lane an output slot, so the per-module stubs[] write order is
    // the lexicographic pair order, independent of how the lanes are scheduled.
    class FormStubsVectorHitStyleKernel {
    public:
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::reco::OTRecHitsConstView hits,
                                    ::reco::OTHitModuleConstView moduleView,
                                    ::reco::StackedModuleGeometryConstView geometry,
                                    float const* maxWidthBarrelFlat,
                                    float const* maxWidthBarrelTilted,
                                    float const* maxWidthEndcap,
                                    int32_t const* barrelFlatMaxCSDiff,
                                    int32_t const* barrelTiltedMaxCSDiff,
                                    int32_t const* endcapMaxCSDiff,
                                    int32_t const* barrelFlatMaxCS,
                                    int32_t const* barrelTiltedMaxCS,
                                    int32_t const* endcapMaxCS,
                                    int32_t const* barrelFlatMaxCSSum,
                                    int32_t const* barrelTiltedMaxCSSum,
                                    int32_t const* endcapMaxCSSum,
                                    uint32_t const* stubOffsets,
                                    ::reco::StubsView stubs,
                                    uint32_t nModules) const {
        const int32_t warpSize = alpaka::warp::getSize(acc);
        const int32_t laneId = static_cast<int32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[1u]);

        // One warp per module
        for (uint32_t iModule : cms::alpakatools::uniform_elements_y(acc, nModules)) {
          // Get hit ranges (uniform across the warp)
          uint32_t hitStart = moduleView[iModule].moduleStart();
          uint32_t hitEnd = moduleView[iModule + 1].moduleStart();
          uint32_t upperStart = moduleView[iModule].upperSensorStart();

          uint32_t nLower = upperStart - hitStart;
          uint32_t nUpper = hitEnd - upperStart;

          // Get geometry parameters (uniform across the warp)
          float separation = geometry[iModule].sensorSeparation();  // in mm
          uint8_t moduleType = geometry[iModule].moduleType();
          bool isPS = (moduleType == 0 || moduleType == 1);
          bool isBarrel = geometry[iModule].isBarrel();
          bool isFlat = geometry[iModule].isFlat();
          bool isFlipped = geometry[iModule].isFlipped();
          uint8_t layer = geometry[iModule].layer();
          float globalLowUpNormX = geometry[iModule].globalLowUpNormX();
          float globalLowUpNormY = geometry[iModule].globalLowUpNormY();
          float globalLowUpNormZ = geometry[iModule].globalLowUpNormZ();
          float localXInGlobalX = geometry[iModule].localXInGlobalX();
          float localXInGlobalY = geometry[iModule].localXInGlobalY();
          float localXInGlobalZ = geometry[iModule].localXInGlobalZ();

          float cut =
              isBarrel ? (isFlat ? maxWidthBarrelFlat[layer] : maxWidthBarrelTilted[layer]) : maxWidthEndcap[layer];

          // Per-layer cluster size cuts
          int32_t maxCSDiff =
              isBarrel ? (isFlat ? barrelFlatMaxCSDiff[layer] : barrelTiltedMaxCSDiff[layer]) : endcapMaxCSDiff[layer];
          int32_t maxCS = isBarrel ? (isFlat ? barrelFlatMaxCS[layer] : barrelTiltedMaxCS[layer]) : endcapMaxCS[layer];
          int32_t maxCSSum =
              isBarrel ? (isFlat ? barrelFlatMaxCSSum[layer] : barrelTiltedMaxCSSum[layer]) : endcapMaxCSSum[layer];

          uint32_t stubBase = stubOffsets[iModule];
          uint32_t maxStubIdx = stubOffsets[iModule + 1];  // End offset for this module

          // Handle PHitOnly recovery for PS modules with P-hits but no S-hits
          if (nUpper == 0) {
            if (isPS && nLower > 0) {
              // Create PHitOnly stub entries for each P-hit (lanes stride across P-hits)
              // These have position but no bend measurement.
              for (uint32_t k = static_cast<uint32_t>(laneId); k < nLower; k += warpSize) {
                uint32_t out = stubBase + k;
                if (out >= maxStubIdx)
                  continue;
                uint32_t iLower = hitStart + k;
                auto const& pHit = hits[iLower];

                float xg = pHit.xGlobal();
                float yg = pHit.yGlobal();
                stubs[out].iphi() = unsafe_atan2s<7>(yg, xg);

                // Stub-specific fields set to invalid values (no bend measurement).
                // Negative dPhiDrError tags PHitOnly entries (see reco::isStub).
                stubs[out].dPhiDr() = 0.0f;
                stubs[out].dPhiDrError() = -1.0f;
                stubs[out].dPhiDrErrorPrec() = -1.0f;  // same sentinel: no bend measurement

                stubs[out].lowerHitIdx() = iLower;
                stubs[out].upperHitIdx() = UINT32_MAX;  // Invalid marker

                stubs[out].flags() = ::reco::StubFlags::makeFlags(isBarrel, isFlat, true, layer, isPS);
              }
            }
            continue;
          }

          // Skip if no lower hits
          if (nLower == 0)
            continue;

          // Topological lower->upper normal: globalLowUpNorm is the physical inner->outer
          // direction, while the width calculation needs the topological lower->upper direction,
          // so the sign is flipped for flipped modules.
          float flipSign = isFlipped ? -1.0f : 1.0f;
          float topoLowUpNormX = flipSign * globalLowUpNormX;
          float topoLowUpNormY = flipSign * globalLowUpNormY;
          float topoLowUpNormZ = flipSign * globalLowUpNormZ;
          float separation_cm = separation * 0.1f;

          // Tilt parameters and sensor frames (uniform across the warp)
          float sinTilt = geometry[iModule].sinTilt();
          float cosTilt = geometry[iModule].cosTilt();
          auto const innerFrame =
              isFlipped ? geometry[iModule].upperSensorFrame() : geometry[iModule].lowerSensorFrame();
          auto const outerFrame =
              isFlipped ? geometry[iModule].lowerSensorFrame() : geometry[iModule].upperSensorFrame();

          // Warp-cooperative formation. Lanes stride across the flattened pair index
          // p = (iLower - hitStart) * nUpper + (iUpper - upperStart) in their natural
          // (iLower, iUpper) order. A warp prefix scan over the per-lane
          // result assigns each passing lane an output slot, preserving the original
          // ordering of stubs[] writes.
          uint32_t nPairs = nLower * nUpper;
          uint32_t writeIdx = stubBase;  // running output cursor, uniform across the warp

          for (uint32_t pBase = 0; pBase < nPairs; pBase += warpSize) {
            uint32_t p = pBase + static_cast<uint32_t>(laneId);
            bool predicate = false;
            uint32_t iLower = 0;
            uint32_t iUpper = 0;

            if (p < nPairs) {
              iLower = hitStart + p / nUpper;
              iUpper = upperStart + p % nUpper;

              auto const& lowerHit = hits[iLower];
              auto const& upperHit = hits[iUpper];

              // Evaluate the VectorHit-style predicate. Use a do/while(false) to
              // simulate labelled break out of the predicate block without exiting
              // the warp-collective code below.
              do {
                if (lowerHit.yLocal() * upperHit.yLocal() < 0.f)
                  break;
                if (alpaka::math::abs(acc, int32_t(lowerHit.clusterSize()) - int32_t(upperHit.clusterSize())) >
                    maxCSDiff)
                  break;
                if (int32_t(lowerHit.clusterSize()) > maxCS || int32_t(upperHit.clusterSize()) > maxCS)
                  break;
                if (int32_t(lowerHit.clusterSize()) + int32_t(upperHit.clusterSize()) > maxCSSum)
                  break;

                float lx_lower = lowerHit.xLocal();
                float lx_upper = upperHit.xLocal();

                float pC = computeParallaxCorrection(acc,
                                                     lowerHit.xGlobal(),
                                                     lowerHit.yGlobal(),
                                                     lowerHit.zGlobal(),
                                                     topoLowUpNormX,
                                                     topoLowUpNormY,
                                                     topoLowUpNormZ,
                                                     localXInGlobalX,
                                                     localXInGlobalY,
                                                     localXInGlobalZ,
                                                     separation_cm);

                float lpos_lower_corr = 0.0f;
                float lpos_upper_corr = 0.0f;
                if (lx_upper > lx_lower) {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
                  } else {
                    lpos_lower_corr = lx_lower + alpaka::math::abs(acc, pC);
                    lpos_upper_corr = lx_upper;
                  }
                } else if (lx_upper < lx_lower) {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower - alpaka::math::abs(acc, pC);
                    lpos_upper_corr = lx_upper;
                  } else {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
                  }
                } else {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
                  } else {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
                  }
                }
                float width = lpos_lower_corr - lpos_upper_corr;
                if (alpaka::math::abs(acc, width) >= cut)
                  break;

                predicate = true;
              } while (false);
            }

            // Warp exclusive prefix scan of (predicate ? 1 : 0).
            // All lanes must execute the shuffles regardless of their predicate or the
            // bound check, so the warp stays converged through this section.
            int32_t myCount = predicate ? 1 : 0;
            int32_t inclusive = myCount;
            for (int32_t off = 1; off < warpSize; off <<= 1) {
              int32_t y = alpaka::warp::shfl(acc, inclusive, laneId - off);
              if (laneId >= off)
                inclusive += y;
            }
            int32_t myOffset = inclusive - myCount;
            int32_t totalThisIter = alpaka::warp::shfl(acc, inclusive, warpSize - 1);

            if (predicate) {
              uint32_t out = writeIdx + static_cast<uint32_t>(myOffset);
              if (out < maxStubIdx) {
                // Determine physical inner/outer for stub properties (position, hit
                // indices, errors).
                uint32_t innerIdx = isFlipped ? iUpper : iLower;
                uint32_t outerIdx = isFlipped ? iLower : iUpper;

                float xi = hits[innerIdx].xGlobal();
                float yi = hits[innerIdx].yGlobal();
                float xo = hits[outerIdx].xGlobal();
                float yo = hits[outerIdx].yGlobal();
                float ri2 = xi * xi + yi * yi;
                float ro2 = xo * xo + yo * yo;

                // Contract: a stub stores exactly one sensor's position+error.
                //   PS modules: always the PIXEL sensor (better spatial precision).
                //     PSP (moduleType == 0): pixel is topologically lower.
                //     PSS (moduleType == 1): pixel is topologically upper.
                //   SS modules: the physically-inner sensor.
                uint32_t pickedIdx = isPS ? ((moduleType == 0) ? iLower : iUpper) : innerIdx;
                float xg = hits[pickedIdx].xGlobal();
                float yg = hits[pickedIdx].yGlobal();
                float zg = hits[pickedIdx].zGlobal();
#ifdef PS_STUB_POSITION_DEBUG
                if (isPS) {
                  float stub_r = alpaka::math::sqrt(acc, xg * xg + yg * yg);
                  printf(
                      "PS_STUB_VH module=%u stubIdx=%u pickedIdx=%u moduleType=%u: "
                      "x=%.6f y=%.6f z=%.6f r=%.4f "
                      "xerrLocal=%.4e yerrLocal=%.4e sensorDetId=%u\n",
                      iModule,
                      out,
                      pickedIdx,
                      unsigned(moduleType),
                      xg,
                      yg,
                      zg,
                      stub_r,
                      hits[pickedIdx].xerrLocal(),
                      hits[pickedIdx].yerrLocal(),
                      hits[pickedIdx].sensorDetId());
                }
#endif
                float rg = alpaka::math::sqrt(acc, xg * xg + yg * yg);

                // iPhi - same encoding as pixel hits (16-bit signed, full circle = [-32768, 32767])
                stubs[out].iphi() = unsafe_atan2s<7>(yg, xg);

                // dPhiDr from the inner/outer global positions, with the parallax correction
                // applied for tilted barrel and endcap modules (in the flat barrel local-x is
                // along phi and the correction vanishes).
                float phi_inner = alpaka::math::atan2(acc, yi, xi);
                float phi_outer = alpaka::math::atan2(acc, yo, xo);
                float ri = alpaka::math::sqrt(acc, ri2);
                float ro = alpaka::math::sqrt(acc, ro2);

                float dphi_raw = phi_outer - phi_inner;
                if (dphi_raw > M_PI)
                  dphi_raw -= 2.0f * M_PI;
                if (dphi_raw < -M_PI)
                  dphi_raw += 2.0f * M_PI;

                float dphi;
                bool applyDphiParallax = !isFlat || !isBarrel;  // tilted barrel or endcap
                if (applyDphiParallax) {
                  float zi = hits[innerIdx].zGlobal();
                  float pC_dphi = computeParallaxCorrection(acc,
                                                            xi,
                                                            yi,
                                                            zi,
                                                            globalLowUpNormX,
                                                            globalLowUpNormY,
                                                            globalLowUpNormZ,
                                                            localXInGlobalX,
                                                            localXInGlobalY,
                                                            localXInGlobalZ,
                                                            separation_cm);
                  float dphi_parallax = (ri > 1e-6f) ? pC_dphi / ri : 0.0f;
                  dphi = dphi_raw - dphi_parallax;
                } else {
                  dphi = dphi_raw;
                }

                float dr_geometric = ro - ri;
                float denominator = cosTilt + sinTilt * zg / rg;

                float dr_effective;
                if (alpaka::math::abs(acc, denominator) > 1e-6f) {
                  dr_effective = separation_cm / denominator;
                } else {
                  dr_effective = dr_geometric;
                }

                stubs[out].dPhiDr() = (alpaka::math::abs(acc, dr_effective) > 1e-6f) ? dphi / dr_effective : 0.0f;

                // dPhiDr error: propagate per-sensor local errors to global via the
                // sensor frames, then form the global phi variance.
                float ge_i[6], ge_o[6];
                innerFrame.toGlobal(hits[innerIdx].xerrLocal(), 0.0f, hits[innerIdx].yerrLocal(), ge_i);
                outerFrame.toGlobal(hits[outerIdx].xerrLocal(), 0.0f, hits[outerIdx].yerrLocal(), ge_o);

                float sinphi_i = alpaka::math::sin(acc, phi_inner);
                float cosphi_i = alpaka::math::cos(acc, phi_inner);
                float sig2phi_i = (ri2 > 0.0f) ? (sinphi_i * sinphi_i * ge_i[0] - 2.0f * sinphi_i * cosphi_i * ge_i[1] +
                                                  cosphi_i * cosphi_i * ge_i[2]) /
                                                     ri2
                                               : 0.0f;
                float sinphi_o = alpaka::math::sin(acc, phi_outer);
                float cosphi_o = alpaka::math::cos(acc, phi_outer);
                float sig2phi_o = (ro2 > 0.0f) ? (sinphi_o * sinphi_o * ge_o[0] - 2.0f * sinphi_o * cosphi_o * ge_o[1] +
                                                  cosphi_o * cosphi_o * ge_o[2]) /
                                                     ro2
                                               : 0.0f;
                float dphi_err = alpaka::math::sqrt(acc, sig2phi_i + sig2phi_o);
                stubs[out].dPhiDrError() = (alpaka::math::abs(acc, dr_effective) > 1e-6f)
                                               ? dphi_err / alpaka::math::abs(acc, dr_effective)
                                               : 0.0f;

                // Precision-only azimuthal bend error (dPhiDrErrorPrec). dPhiDrError above adds the
                // full local covariance of both sensors as independent, so the along-strip (local-y)
                // term enters twice. In a 2S stack the two sensors carry identical, aligned strips:
                // the reported along-strip position is the same physical offset in both sensors, it
                // is common mode in dphi = phi_outer - phi_inner, and only the difference of the two
                // projections survives. The term leaks into the azimuth only in the endcap, where the
                // strips are radial just at the module centre; in the barrel local-y is z-hat and the
                // correction is a no-op. Form: sigma2 = ux_i^2 vx_i + ux_o^2 vx_o + (sy_i uy_i - sy_o uy_o)^2
                // with ux, uy the azimuthal angular sensitivity to a unit displacement along local-x,
                // local-y, taken from the sensor frame rotation (same convention as SOAFrame::toGlobal):
                // an in-plane displacement (dX, dY) moves the azimuth by (cos(phi) dY - sin(phi) dX) / r.
                // dPhiDrError is kept as is because the CA stub-curvature cuts are tuned on it.
                {
                  auto const& rotI = innerFrame.rotation();
                  auto const& rotO = outerFrame.rotation();
                  const float riInv = (ri > 1e-6f) ? 1.0f / ri : 0.0f;
                  const float roInv = (ro > 1e-6f) ? 1.0f / ro : 0.0f;
                  // azimuthal angular sensitivity to a unit displacement along each local axis
                  const float uxI = (cosphi_i * rotI.xy() - sinphi_i * rotI.xx()) * riInv;
                  const float uxO = (cosphi_o * rotO.xy() - sinphi_o * rotO.xx()) * roInv;
                  const float uyI = (cosphi_i * rotI.yy() - sinphi_i * rotI.yx()) * riInv;
                  float uyO = (cosphi_o * rotO.yy() - sinphi_o * rotO.yx()) * roInv;
                  // orient the outer strip direction onto the inner one (common-mode requires a
                  // common global sense; the local frames' y-axes may be flipped relative to each
                  // other without any physical meaning)
                  const float dotY = rotI.yx() * rotO.yx() + rotI.yy() * rotO.yy() + rotI.yz() * rotO.yz();
                  if (dotY < 0.0f)
                    uyO = -uyO;
                  const float vxI = alpaka::math::max(acc, hits[innerIdx].xerrLocal(), 0.0f);
                  const float vxO = alpaka::math::max(acc, hits[outerIdx].xerrLocal(), 0.0f);
                  const float syI = alpaka::math::sqrt(acc, alpaka::math::max(acc, hits[innerIdx].yerrLocal(), 0.0f));
                  const float syO = alpaka::math::sqrt(acc, alpaka::math::max(acc, hits[outerIdx].yerrLocal(), 0.0f));
                  // The common-mode cancellation holds only when both sensors report the same
                  // along-strip segmentation, i.e. a 2S/2S stack with identical, aligned strips. In a
                  // PS stack a macro-pixel column faces a 2.4 cm strip, the two quantisations are
                  // independent and the plain sum is the correct treatment. The test compares the two
                  // along-strip variances (the strip CPE returns L^2/12, a per-sensor-type constant,
                  // so identical sensors give bit-equal values); 1e-6 is a float tolerance.
                  const bool sameSegmentation =
                      alpaka::math::abs(acc, hits[innerIdx].yerrLocal() - hits[outerIdx].yerrLocal()) <=
                      1.0e-6f * (hits[innerIdx].yerrLocal() + hits[outerIdx].yerrLocal());
                  const float alongVar = sameSegmentation ? (syI * uyI - syO * uyO) * (syI * uyI - syO * uyO)
                                                          : (syI * uyI * syI * uyI + syO * uyO * syO * uyO);
                  const float sig2prec = vxI * uxI * uxI + vxO * uxO * uxO + alongVar;
                  const float errPrec = alpaka::math::sqrt(acc, alpaka::math::max(acc, sig2prec, 0.0f));
                  stubs[out].dPhiDrErrorPrec() = (alpaka::math::abs(acc, dr_effective) > 1e-6f)
                                                     ? errPrec / alpaka::math::abs(acc, dr_effective)
                                                     : 0.0f;
                }

                stubs[out].lowerHitIdx() = iLower;
                stubs[out].upperHitIdx() = iUpper;
                stubs[out].flags() = ::reco::StubFlags::makeFlags(isBarrel, isFlat, true, layer, isPS);
              }
            }

            writeIdx += static_cast<uint32_t>(totalThisIter);
          }
        }
      }
    };

    // Kernel to fill moduleStart array from stubOffsets
    class FillModuleStartKernel {
    public:
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    uint32_t const* stubOffsets,
                                    ::reco::StubModuleView stubModuleView,
                                    uint32_t nModules) const {
        // One thread per module (+ 1 for the final total)
        for (auto i : cms::alpakatools::uniform_elements(acc, nModules + 1)) {
          stubModuleView[i].moduleStart() = stubOffsets[i];
        }
      }
    };

#ifdef STUB_DIAGNOSTIC_COUNTERS
    // Diagnostic kernel: Compute event-level summary statistics for stub formation (VectorHit style)
    // This kernel runs after stub formation to aggregate statistics across all modules
    class DiagnosticSummaryKernel {
    public:
      template <typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::reco::OTRecHitsConstView hits,
                                    ::reco::OTHitModuleConstView moduleView,
                                    ::reco::StackedModuleGeometryConstView geometry,
                                    float const* maxWidthBarrelFlat,
                                    float const* maxWidthBarrelTilted,
                                    float const* maxWidthEndcap,
                                    uint32_t nModules) const {
        // Single thread computes event-level summary
        // This is inefficient but acceptable for diagnostic purposes
        for (auto idx : cms::alpakatools::uniform_elements(acc, 1u)) {
          (void)idx;  // Unused, we just need one thread

          uint32_t totalPHits_PS = 0;
          uint32_t totalPHits_SS = 0;
          uint32_t pHitsWithStub_PS = 0;
          uint32_t pHitsWithStub_SS = 0;
          uint32_t pHitsLostNoUpper_PS = 0;
          uint32_t pHitsLostNoUpper_SS = 0;
          uint32_t pHitsLostWidthCut_PS = 0;
          uint32_t pHitsLostWidthCut_SS = 0;
          uint32_t totalStubs = 0;
          uint32_t maxStubsPerPHit_PS = 0;
          uint32_t maxStubsPerPHit_SS = 0;
          uint32_t pHitsMultiStub_PS = 0;  // P-hits that formed >1 stub
          uint32_t pHitsMultiStub_SS = 0;

          for (uint32_t iModule = 0; iModule < nModules; ++iModule) {
            uint32_t hitStart = moduleView[iModule].moduleStart();
            uint32_t hitEnd = moduleView[iModule + 1].moduleStart();
            uint32_t upperStart = moduleView[iModule].upperSensorStart();

            uint32_t nLower = upperStart - hitStart;
            uint32_t nUpper = hitEnd - upperStart;

            uint8_t moduleType = geometry[iModule].moduleType();
            bool isPS = (moduleType == 0 || moduleType == 1);
            bool isBarrel = geometry[iModule].isBarrel();
            bool isFlipped = geometry[iModule].isFlipped();
            uint8_t layer = geometry[iModule].layer();
            float separation = geometry[iModule].sensorSeparation();
            float globalLowUpNormX = geometry[iModule].globalLowUpNormX();
            float globalLowUpNormY = geometry[iModule].globalLowUpNormY();
            float globalLowUpNormZ = geometry[iModule].globalLowUpNormZ();
            float localXInGlobalX = geometry[iModule].localXInGlobalX();
            float localXInGlobalY = geometry[iModule].localXInGlobalY();
            float localXInGlobalZ = geometry[iModule].localXInGlobalZ();
            bool isFlat = geometry[iModule].isFlat();
            float cut =
                isBarrel ? (isFlat ? maxWidthBarrelFlat[layer] : maxWidthBarrelTilted[layer]) : maxWidthEndcap[layer];

            if (isPS) {
              totalPHits_PS += nLower;
            } else {
              totalPHits_SS += nLower;
            }

            if (nLower == 0 || nUpper == 0) {
              // Lower hits lost due to no upper hits
              if (nLower > 0 && nUpper == 0) {
                if (isPS) {
                  pHitsLostNoUpper_PS += nLower;
                } else {
                  pHitsLostNoUpper_SS += nLower;
                }
              }
              continue;
            }

            // Count stubs and track per-lower-hit statistics
            for (uint32_t iLower = hitStart; iLower < upperStart; ++iLower) {
              auto const& lowerHit = hits[iLower];
              uint32_t stubsForThisLowerHit = 0;

              for (uint32_t iUpper = upperStart; iUpper < hitEnd; ++iUpper) {
                auto const& upperHit = hits[iUpper];

                // Same-sign local-y cut: reject if hits are on opposite strip sides
                if (lowerHit.yLocal() * upperHit.yLocal() < 0.f)
                  continue;

                // VectorHit style width calculation (simplified for summary)
                float lx_lower = lowerHit.xLocal();
                float lx_upper = upperHit.xLocal();

                // Parallax correction
                float flipSign = isFlipped ? -1.0f : 1.0f;
                float topoLowUpNormX = flipSign * globalLowUpNormX;
                float topoLowUpNormY = flipSign * globalLowUpNormY;
                float topoLowUpNormZ = flipSign * globalLowUpNormZ;

                float pC = computeParallaxCorrection(acc,
                                                     lowerHit.xGlobal(),
                                                     lowerHit.yGlobal(),
                                                     lowerHit.zGlobal(),
                                                     topoLowUpNormX,
                                                     topoLowUpNormY,
                                                     topoLowUpNormZ,
                                                     localXInGlobalX,
                                                     localXInGlobalY,
                                                     localXInGlobalZ,
                                                     separation * 0.1f);

                // Apply parallax correction
                float lpos_lower_corr = 0.0f;
                float lpos_upper_corr = 0.0f;

                if (lx_upper > lx_lower) {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
                  } else {
                    lpos_lower_corr = lx_lower + alpaka::math::abs(acc, pC);
                    lpos_upper_corr = lx_upper;
                  }
                } else if (lx_upper < lx_lower) {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower - alpaka::math::abs(acc, pC);
                    lpos_upper_corr = lx_upper;
                  } else {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
                  }
                } else {
                  if (lx_upper > 0) {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper - alpaka::math::abs(acc, pC);
                  } else {
                    lpos_lower_corr = lx_lower;
                    lpos_upper_corr = lx_upper + alpaka::math::abs(acc, pC);
                  }
                }

                float width = lpos_lower_corr - lpos_upper_corr;

                if (alpaka::math::abs(acc, width) < cut) {
                  stubsForThisLowerHit++;
                  totalStubs++;
                }
              }

              if (stubsForThisLowerHit > 0) {
                if (isPS) {
                  pHitsWithStub_PS++;
                  if (stubsForThisLowerHit > maxStubsPerPHit_PS) {
                    maxStubsPerPHit_PS = stubsForThisLowerHit;
                  }
                  if (stubsForThisLowerHit > 1) {
                    pHitsMultiStub_PS++;
                  }
                } else {
                  pHitsWithStub_SS++;
                  if (stubsForThisLowerHit > maxStubsPerPHit_SS) {
                    maxStubsPerPHit_SS = stubsForThisLowerHit;
                  }
                  if (stubsForThisLowerHit > 1) {
                    pHitsMultiStub_SS++;
                  }
                }
              } else {
                // Lost to width cut
                if (isPS) {
                  pHitsLostWidthCut_PS++;
                } else {
                  pHitsLostWidthCut_SS++;
                }
              }
            }
          }

          // Print event-level summary
          printf(
              "STUB_DIAG_EVENT_VH type=PS totalPHits=%u pHitsWithStub=%u pHitsLost=%u "
              "lostNoUpper=%u lostWidthCut=%u totalStubs=%u maxStubsPerPHit=%u pHitsMultiStub=%u\n",
              totalPHits_PS,
              pHitsWithStub_PS,
              pHitsLostNoUpper_PS + pHitsLostWidthCut_PS,
              pHitsLostNoUpper_PS,
              pHitsLostWidthCut_PS,
              totalStubs,
              maxStubsPerPHit_PS,
              pHitsMultiStub_PS);
          printf(
              "STUB_DIAG_EVENT_VH type=SS totalPHits=%u pHitsWithStub=%u pHitsLost=%u "
              "lostNoUpper=%u lostWidthCut=%u totalStubs=%u maxStubsPerPHit=%u pHitsMultiStub=%u\n",
              totalPHits_SS,
              pHitsWithStub_SS,
              pHitsLostNoUpper_SS + pHitsLostWidthCut_SS,
              pHitsLostNoUpper_SS,
              pHitsLostWidthCut_SS,
              totalStubs,
              maxStubsPerPHit_SS,
              pHitsMultiStub_SS);
          printf("STUB_DIAG_EVENT_VH type=ALL totalPHits=%u pHitsWithStub=%u pHitsLost=%u totalStubs=%u\n",
                 totalPHits_PS + totalPHits_SS,
                 pHitsWithStub_PS + pHitsWithStub_SS,
                 pHitsLostNoUpper_PS + pHitsLostWidthCut_PS + pHitsLostNoUpper_SS + pHitsLostWidthCut_SS,
                 totalStubs);
        }
      }
    };
#endif  // STUB_DIAGNOSTIC_COUNTERS

  }  // namespace otStubFormationVectorHitStyle

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoTracker_PixelSeeding_plugins_alpaka_OTStubFormationVectorHitStyleKernels_h
