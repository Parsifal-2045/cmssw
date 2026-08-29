// Split-build TU: fused CA main fit, fast phase, Phase2OTStubs (N=3..10). One phase per TU: a fused
// kernel pulls all N of its phase, so separating phases bounds each nvcc partition to the per-N size
// (cf. BrokenLineFit_stubsMainBL). The per-bin reference kernels stay in stubsMainBL.
#include "BrokenLineFitKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  BLFIT_MAIN_FUSED_FAST_SIG(Phase2OTStubs);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
