// Split-build TU: fused extended-N refit, Outlier phase. One phase per TU: a fused kernel pulls all ten N,
// so separating phases bounds each nvcc partition (cf. per-N split in BrokenLineFit_refitLo/Hi).
#include "BrokenLineFitKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  BLFIT_REFIT_FUSED_OUTLIER_SIG();

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
