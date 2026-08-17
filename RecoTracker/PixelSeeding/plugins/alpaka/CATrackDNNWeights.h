#ifndef RecoTracker_PixelSeeding_plugins_alpaka_CATrackDNNWeights_h
#define RecoTracker_PixelSeeding_plugins_alpaka_CATrackDNNWeights_h

// Dispatcher for the per-iteration loose/tight track-DNN weight banks:
//   - CATrackDNNWeights_displaced.h : namespace caTrackDNN_displaced (default bank)
//   - CATrackDNNWeights_prompt.h    : namespace caTrackDNN_prompt
// The bank is chosen at runtime (DnnBank, CADnnBank.h); caTrackDNN_eval::score<BANK>() is specialised per bank.
// The `caTrackDNN` alias is the displaced bank, for code reading kNFeat/kDefaultThreshold without a bank.

#include "CATrackDNNWeights_displaced.h"
#include "CATrackDNNWeights_prompt.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace caTrackDNN = caTrackDNN_displaced;  // default bank alias

  static_assert(caTrackDNN_prompt::kNFeat == caTrackDNN_displaced::kNFeat &&
                    caTrackDNN_prompt::kNH1 == caTrackDNN_displaced::kNH1 &&
                    caTrackDNN_prompt::kNH2 == caTrackDNN_displaced::kNH2,
                "prompt and displaced loose/tight-DNN banks must share the network shape");
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // RecoTracker_PixelSeeding_plugins_alpaka_CATrackDNNWeights_h
