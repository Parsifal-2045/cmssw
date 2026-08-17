#ifndef RecoTracker_PixelSeeding_interface_StackedModuleGeometrySoA_h
#define RecoTracker_PixelSeeding_interface_StackedModuleGeometrySoA_h

#include <cstdint>

#include "DataFormats/GeometrySurface/interface/SOARotation.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace reco {

  // Stacked Module Geometry SoA: Phase-2 OT stacked-module geometry for stub formation.
  // Each module has two closely-spaced sensors (inner/outer). Module types:
  //   PS (Pixel-Strip): inner=macro-pixel, outer=strip, precise z; SS (Strip-Strip): no precise z.
  GENERATE_SOA_LAYOUT(
      StackedModuleGeometryLayout,
      // DetId as uint32_t
      SOA_COLUMN(uint32_t, detId),
      // Partner DetId as uint32_t
      SOA_COLUMN(uint32_t, partnerDetId),
      // Stacked DetId as uint32_t
      SOA_COLUMN(uint32_t, stackedDetId),
      // Original index in the TrackerGeometry detUnit vector
      SOA_COLUMN(uint32_t, geomIndex),
      // kind of module: 0-P, 1-SinPS, 2-SS
      SOA_COLUMN(uint8_t, moduleType),
      // Sensor separation: distance between inner and outer sensors (mm)
      // Typical values: 2.0-4.0mm for PS, 1.6-4.0mm for SS
      SOA_COLUMN(float, sensorSeparation),

      // Mean radius of the module (cm)
      SOA_COLUMN(float, meanRadius),

      // Tilt angle (rad): module axis from radial direction, atan2(dz_phys, dr_phys) using the
      // PHYSICAL inner->outer direction (NOT topological lower->upper) so dPhiDr has a consistent
      // sign regardless of flipped status. Barrel flat ~0, tilted ~nonzero (toward IP), endcap ~+/-pi/2.
      SOA_COLUMN(float, tiltAngle),

      // Precomputed sin and cos of tilt angle for effective dr calculation
      // Used in dr_effective = separation / (cosTilt + sinTilt * z/r)
      // Note: tiltAngle = atan2(dz, dr) measured from radial axis, so sinTilt > 0 for +z tilt
      SOA_COLUMN(float, sinTilt),
      SOA_COLUMN(float, cosTilt),

      // Module type: true if PS (pixel-strip), false if SS (strip-strip)
      SOA_COLUMN(bool, isPS),

      // isFlipped: topological "lower" sensor is physically farther from the beam; determines which
      // sensor is physically inner/outer for stub direction.
      SOA_COLUMN(bool, isFlipped),

      // Detector location: true for barrel, false for endcap
      SOA_COLUMN(bool, isBarrel),

      // Module is flat (not tilted) - only meaningful for barrel
      // Flat means module axis is nearly radial (|cos(tiltAngle)| > cos(0.1) ~= 0.995)
      // Corresponds to tilt angles near 0 (inner->outer nearly radial)
      SOA_COLUMN(bool, isFlat),

      // Forward endcap: true if z > 0 (forward), false if z < 0 (backward)
      // Only meaningful for endcap modules (isBarrel == false)
      SOA_COLUMN(bool, isFwdEndcap),

      // Layer number within OT (0-5, corresponding to TBPS/TBSS layers)
      SOA_COLUMN(uint8_t, layer),

      // Maximum allowed bend for pT threshold (radians)
      // Calculated as: maxBend = 0.3 * B * separation / minPt
      // Stubs with |bend| > maxBend are rejected
      SOA_COLUMN(float, maxBend),

      // Normalized physical inner->outer sensor vector in global coords (consistent with tiltAngle,
      // for sign-independent dPhiDr).
      SOA_COLUMN(float, globalLowUpNormX),
      SOA_COLUMN(float, globalLowUpNormY),
      SOA_COLUMN(float, globalLowUpNormZ),

      // Local x-axis direction in global coords (parallax correction).
      SOA_COLUMN(float, localXInGlobalX),
      SOA_COLUMN(float, localXInGlobalY),
      SOA_COLUMN(float, localXInGlobalZ),

      // Per-sensor surface frames (position + rotation). Used by the stub bend-error formula (sigma_phi
      // propagated through each sensor's own frame) and as the source for CAModulesSoA::innerSensorFrame.
      SOA_COLUMN(SOAFrame<float>, lowerSensorFrame),
      SOA_COLUMN(SOAFrame<float>, upperSensorFrame))

  using StackedModuleGeometrySoA = StackedModuleGeometryLayout<>;
  using StackedModuleGeometryView = StackedModuleGeometrySoA::View;
  using StackedModuleGeometryConstView = StackedModuleGeometrySoA::ConstView;

}  // namespace reco

#endif  // RecoTracker_PixelSeeding_interface_StackedModuleGeometrySoA_h
