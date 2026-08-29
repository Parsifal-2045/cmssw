// Split-build TU: fused CA main fit, 4 non-stubs traits (N=3..6), both phases. Phase2OTStubs in
// BrokenLineFit_mainFusedStubs{Fast,Fit}.dev.cc. Each (phase,traits) in exactly one TU (ODR).
#include "BrokenLineFitKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  BLFIT_MAIN_FUSED_FAST_SIG(Phase1);
  BLFIT_MAIN_FUSED_FAST_SIG(Phase2);
  BLFIT_MAIN_FUSED_FAST_SIG(Phase2OT);
  BLFIT_MAIN_FUSED_FAST_SIG(HIonPhase1);

  BLFIT_MAIN_FUSED_FIT_SIG(Phase1);
  BLFIT_MAIN_FUSED_FIT_SIG(Phase2);
  BLFIT_MAIN_FUSED_FIT_SIG(Phase2OT);
  BLFIT_MAIN_FUSED_FIT_SIG(HIonPhase1);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
