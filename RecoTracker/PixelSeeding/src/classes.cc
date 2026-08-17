#include "DataFormats/Portable/interface/PortableHostCollectionReadRules.h"
#include "RecoTracker/PixelSeeding/interface/TripletDumpHost.h"

// ROOT read rules for the per-built-triplet training-dataset host SoA (CA_TRIPLET_DUMP product),
// as for any single-layout PortableHostCollection. Type-info only, harmless when the dump is off.
SET_PORTABLEHOSTCOLLECTION_READ_RULES(PortableHostCollection<caStructures::TripletDumpSoA>);
