#include "RecoTracker/PixelSeeding/interface/IntermediateHitTriplets.h"
#include "DataFormats/Common/interface/Wrapper.h"

// Per-built-triplet training-dataset host SoA (CA_TRIPLET_DUMP product): lets the device product be
// transcribed to a host collection read by TripletFeaturesTableProducer. Type-info only, harmless
// when CA_TRIPLET_DUMP is off.
#include "RecoTracker/PixelSeeding/interface/TripletDumpHost.h"

// CA-ordered module geometry host collection (EventSetup product of CAGeometryESProducer, read
// through the CAGeometryRecord): needed for the serial backend product and the device product's
// host transcription. Type-info only.
#include "RecoTracker/PixelSeeding/interface/CAGeometrySoA.h"
#include "RecoTracker/PixelSeeding/interface/CAGeometryHost.h"

#include <vector>

namespace RecoPixelVertexing_PixelTriplets {
  struct dictionary {
    IntermediateHitTriplets iht;
    edm::Wrapper<IntermediateHitTriplets> wiht;
  };
}  // namespace RecoPixelVertexing_PixelTriplets
