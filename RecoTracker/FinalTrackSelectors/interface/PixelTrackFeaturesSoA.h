#ifndef RecoTracker_FinalTrackSelectors_PixelTrackFeaturesSoA_h
#define RecoTracker_FinalTrackSelectors_PixelTrackFeaturesSoA_h

#include "DataFormats/SoATemplate/interface/SoALayout.h"

// Column ABI: cols 1-17 are the fit/covariance features, bound positionally to the TorchScript
// model. Cols 18-27 append the 10 CA hit/stub features; 28-31 the hit-topology and radial-extent
// features; 32-35 the merged-collection provenance; 36-42 the pixel-cluster charge/shape block.
// Cols 18-42 are written only when useHitFeatures=true. Every model consumes a PREFIX of this
// layout, so appending a column never changes an existing model's score.
//
// Full 42-column ABI (1-based; the trainer feeds features in this order):
//    1 chi2            2 dzError        3 dxyError       4 eta            5 nHits
//    6 phi             7 phiError       8 pt             9 qOverPtError  10 dzBS
//   11 dxyBS          12 nLayers       13 cotThetaError  14 covCotThetaDz 15 covDxyQOverPt
//   16 covPhiDxy      17 covPhiQOverPt 18 caFitChi2      19 psFrac        20 r0
//   21 nPS            22 spanZ         23 nStubs         24 logChi2Stub   25 kErr
//   26 dcaEst         27 nBarrel       28 rzChi2         29 meanStubKappa 30 leverArm
//   31 rMax           32 nAttached     33 nOTExtras      34 iterationId   35 ndof
//   36 minCharge      37 meanCharge    38 minChargeNorm  39 maxSizeY      40 meanSizeY
//   41 maxSizeX       42 nLowCharge
GENERATE_SOA_LAYOUT(PixelTrackFeaturesSoALayout,
                    SOA_COLUMN(float, chi2),
                    SOA_COLUMN(float, dzError),
                    SOA_COLUMN(float, dxyError),
                    SOA_COLUMN(float, eta),
                    SOA_COLUMN(float, nHits),
                    SOA_COLUMN(float, phi),
                    SOA_COLUMN(float, phiError),
                    SOA_COLUMN(float, pt),
                    SOA_COLUMN(float, qOverPtError),
                    SOA_COLUMN(float, dzBS),
                    SOA_COLUMN(float, dxyBS),
                    SOA_COLUMN(float, nLayers),
                    SOA_COLUMN(float, cotThetaError),
                    SOA_COLUMN(float, covCotThetaDz),
                    SOA_COLUMN(float, covDxyQOverPt),
                    SOA_COLUMN(float, covPhiDxy),
                    SOA_COLUMN(float, covPhiQOverPt),
                    SOA_COLUMN(float, caFitChi2),
                    SOA_COLUMN(float, psFrac),
                    SOA_COLUMN(float, r0),
                    SOA_COLUMN(float, nPS),
                    SOA_COLUMN(float, spanZ),
                    SOA_COLUMN(float, nStubs),
                    SOA_COLUMN(float, logChi2Stub),
                    SOA_COLUMN(float, kErr),
                    SOA_COLUMN(float, dcaEst),
                    SOA_COLUMN(float, nBarrel),
                    // cols 28-29: hit-topology features from caTrackFeatures::fill (rzKappaOut).
                    SOA_COLUMN(float, rzChi2),
                    SOA_COLUMN(float, meanStubKappa),
                    // cols 30-31: radial-extent features (leverArm=rMax-r0, rMax=max hit radius),
                    // from the same fill() walk. Appended to keep cols 1-29 a contiguous prefix.
                    SOA_COLUMN(float, leverArm),
                    SOA_COLUMN(float, rMax),
                    // cols 32-35: merged-collection provenance (useHitFeatures-only). nAttached/
                    // nOTExtras count the per-track CSR hit span; iterationId is constant on a
                    // single-arm collection. ABI order: nAttached, nOTExtras, iterationId, ndof.
                    SOA_COLUMN(float, nAttached),
                    SOA_COLUMN(float, nOTExtras),
                    SOA_COLUMN(float, iterationId),
                    SOA_COLUMN(float, ndof),
                    // cols 36-42: pixel-cluster charge/shape block, aggregated over the track's
                    // pixel hits only. Charge is in electrons; minChargeNorm is path-length
                    // normalised: Q*|sin(theta)| in the pixel barrel (detectorIndex < 864),
                    // Q*|cos(theta)| in the endcap. nLowCharge counts hits with normalised charge
                    // < 7000 e. Cluster sizes are the raw signed 1/8-pixel values (clipped at 127,
                    // negated at sensor edges), used as stored. Sentinel: -1 on all seven columns
                    // for a track with no usable pixel hit (and for padding rows), matching the
                    // nano producer. ABI order: minCharge, meanCharge, minChargeNorm, maxSizeY,
                    // meanSizeY, maxSizeX, nLowCharge.
                    SOA_COLUMN(float, minCharge),
                    SOA_COLUMN(float, meanCharge),
                    SOA_COLUMN(float, minChargeNorm),
                    SOA_COLUMN(float, maxSizeY),
                    SOA_COLUMN(float, meanSizeY),
                    SOA_COLUMN(float, maxSizeX),
                    SOA_COLUMN(float, nLowCharge));

using PixelTrackFeaturesSoA = PixelTrackFeaturesSoALayout<>;

// Width of the feature vector the compact-forest scorer assembles positionally (== number of
// columns in the layout above). The model loader validates every int8 split-feature index
// against it; it is the declared bound of the f[] initializer in the scoring kernels.
inline constexpr int kNPixelTrackFeatures = 42;

// Define the SoA layout for track scores (output)
GENERATE_SOA_LAYOUT(PixelTrackScoresSoALayout, SOA_COLUMN(float, score))

using PixelTrackScoresSoA = PixelTrackScoresSoALayout<>;

#endif
