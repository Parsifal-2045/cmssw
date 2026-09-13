#ifndef RecoTracker_PixelTrackFitting_test_gblTestMaterial_h
#define RecoTracker_PixelTrackFitting_test_gblTestMaterial_h

// Host copy of the Geant4-material-map walk that the GBL fit runs on the device, for gblReplay: the host
// fit (interface/BrokenLine.h) carries no material walk, and the device one is ALPAKA_FN_ACC code a
// plain host tool cannot call. The functions below transcribe the device ones with
// `alpaka::math::sqrt(acc, x)` replaced by `std::sqrt(x)`; same lattice, same cell boundaries and the same
// accumulation order, and testGblReplayDevice asserts the two give the same doubles on every fixture and
// backend.
// The table itself is not duplicated: blMaterialMapData() is the one compiled-in copy.

#include <cmath>

#include "RecoTracker/PixelTrackFitting/interface/BLMaterialMap.h"

namespace gblTestMaterial {

  // Exact cell walk of the material map along the straight (r,z) chord (r0,z0)->(r1,z1): chord length L
  // and the moments about the ARRIVAL end, W = int rho dl, S1 = int rho d dl, S2 = int rho d^2 dl, with
  // `path3D` > 0 rescaling them to the 3-D path (k, k^2, k^3 with k = path3D/L).
  // cf. ALPAKA_ACCELERATOR_NAMESPACE::brokenline::segmentWalk.
  inline void segmentWalk(
      double r0, double z0, double r1, double z1, double path3D, double& L, double& W, double& S1, double& S2) {
    const float* rho = blMaterialMap::blMaterialMapData();
    const double dr = r1 - r0, dz = z1 - z0;
    L = std::sqrt(dr * dr + dz * dz);
    W = S1 = S2 = 0.;
    if (!(L > 0.))
      return;
    double tR = 2., dtR = 1.;
    if (dr != 0.) {
      const double kNext = std::floor(r0 / double(blMaterialMap::kDR)) + (dr > 0. ? 1. : 0.);
      tR = (kNext * double(blMaterialMap::kDR) - r0) / dr;
      dtR = double(blMaterialMap::kDR) / std::abs(dr);
    }
    double tZ = 2., dtZ = 1.;
    if (dz != 0.) {
      const double kNext = std::floor(z0 / double(blMaterialMap::kDZ)) + (dz > 0. ? 1. : 0.);
      tZ = (kNext * double(blMaterialMap::kDZ) - z0) / dz;
      dtZ = double(blMaterialMap::kDZ) / std::abs(dz);
    }
    constexpr int kMaxCells = 2048;
    double t = 0.;
    for (int cell = 0; cell < kMaxCells && t < 1.; ++cell) {
      double tn = tR < tZ ? tR : tZ;
      if (tn > 1.)
        tn = 1.;
      if (tn > t) {
        const double tm = 0.5 * (t + tn);
        const double q = blMaterialMap::rhoAt(rho, float(r0 + tm * dr), float(z0 + tm * dz));
        if (q > 0.f) {
          const double a = 1. - t, c = 1. - tn;
          W += q * (a - c) * L;
          S1 += q * (a * a - c * c) * 0.5 * L * L;
          S2 += q * (a * a * a - c * c * c) * (1. / 3.) * L * L * L;
        }
      }
      t = tn;
      if (tR <= t)
        tR += dtR;
      if (tZ <= t)
        tZ += dtZ;
    }
    if (path3D > 0.) {
      const double k = path3D / L;
      W *= k;
      S1 *= k * k;
      S2 *= k * k * k;
      L = path3D;
    }
  }

  // Two-equivalent-thin-scatterer split of a segment: W, the interior scatterer's path distance from the
  // arrival end d1 = S2/S1 and its share of the variance w1 = S1^2/(S2 W).
  // cf. ALPAKA_ACCELERATOR_NAMESPACE::brokenline::segmentXX0Moments.
  inline double segmentXX0Moments(
      double r0, double z0, double r1, double z1, double& d1, double& w1, double path3D = 0.) {
    double L, W, S1, S2;
    segmentWalk(r0, z0, r1, z1, path3D, L, W, S1, S2);
    d1 = 0.;
    w1 = 0.;
    if (W > 0. && S1 > 0. && S2 > 0.) {
      d1 = S2 / S1;
      w1 = S1 * S1 / (S2 * W);
    }
    return W;
  }

  // The material rows of PreparedGblData<n>.
  template <int N>
  struct MatData {
    double matXX0[N] = {};  // slot g = the WHOLE of gap g->g+1 (slot N-1 is zero)
    double gapD1[N] = {};   // gap g's interior equivalent-scatterer path distance from its arrival hit [cm]
    double gapW1[N] = {};   // gap g's interior equivalent-scatterer share of the variance, in (0,1]
    double innerXX0 = 0.;   // beamline (z = 0) -> first hit, beam pipe + upstream material
    double innerD1 = 0.;    // the same two-thin split for the upstream segment
    double innerW1 = 0.;
  };

  // The material section of brokenline::prepareGblFitData (alpaka/BrokenLine.h), on the host.
  // hits is any 3xN Eigen-like object indexable as hits(row, col): rows 0,1,2 = x,y,z [cm]; sTotal is
  // that track's 3-D arc length, which sets the path each gap's material is scaled to. The upstream term
  // is always integrated, as in the device walk, so OT-only stub tracks carry the upstream material too.
  template <int N, typename M3xN, typename VN>
  inline void fillMatData(const M3xN& hits, const VN& sTotal, MatData<N>& md) {
    auto rOf = [&](int j) {
      return std::sqrt(double(hits(0, j)) * double(hits(0, j)) + double(hits(1, j)) * double(hits(1, j)));
    };
    for (int i = 0; i < N; ++i) {
      md.matXX0[i] = 0.;
      md.gapD1[i] = 0.;
      md.gapW1[i] = 0.;
    }
    for (int g = 0; g + 1 < N; ++g) {
      const double path = std::abs(sTotal(g + 1) - sTotal(g));
      md.matXX0[g] = segmentXX0Moments(rOf(g), hits(2, g), rOf(g + 1), hits(2, g + 1), md.gapD1[g], md.gapW1[g], path);
    }
    md.innerXX0 = segmentXX0Moments(0., 0., rOf(0), hits(2, 0), md.innerD1, md.innerW1);
  }

}  // namespace gblTestMaterial

#endif
