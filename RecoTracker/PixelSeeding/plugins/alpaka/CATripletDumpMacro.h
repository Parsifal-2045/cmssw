#ifndef RecoTracker_PixelSeeding_plugins_alpaka_CATripletDumpMacro_h
#define RecoTracker_PixelSeeding_plugins_alpaka_CATripletDumpMacro_h

// Toggle for the built-triplet training-dataset dump (one row per built triplet carrying
// the in-kernel triplet feature vector, used to train the inline triplet DNN). OFF in
// production: every dump-related branch is #ifdef'd on this macro. Lives in its own
// minimal header so the producer/generator/kernels can see the toggle without pulling in
// the heavy CATripletCuts.h. Recipe: RecoTracker/PixelSeeding/test/models/README.md.
// #define CA_TRIPLET_DUMP

#endif  // RecoTracker_PixelSeeding_plugins_alpaka_CATripletDumpMacro_h
