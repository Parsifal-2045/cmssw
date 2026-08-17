#ifndef RecoTracker_PixelSeeding_plugins_alpaka_CASizingDumpMacro_h
#define RecoTracker_PixelSeeding_plugins_alpaka_CASizingDumpMacro_h

// Toggle for per-event per-container sizing dump (re-fit capacity curves). OFF by
// default: every dump block is #ifdef'd on this macro, so a default build compiles
// none of it. On: each stage prints [CA Sizing] iter=<stage> <key>=<value> per event.
// Lives in its own minimal header so the kernels, extension, and merger share the toggle.
// #define CA_SIZING_DUMP

#endif  // RecoTracker_PixelSeeding_plugins_alpaka_CASizingDumpMacro_h
