/**
 * PixelTrackForestHighPuritySelector
 *
 * GPU/Accelerator module performing displaced HighPurity pixel-track selection with a
 * gradient-boosted decision tree ("forest"), composed of:
 *   1. CA-based quality preselection
 *   2. Feature extraction
 *   3. Compact gradient-boosted tree inference (custom alpaka traversal kernel)
 *   4. Score-based filtering
 *   5. Track/hit compaction and output production
 *
 * This is the tree counterpart of PixelTrackTorchHighPuritySelector (which runs a Torch MLP via
 * a per-stream batched forward). A Torch forest would be replicated per stream, which is not
 * affordable at O(100) streams per GPU; this module instead shares one read-only copy per device.
 *
 * Pipeline: TracksSoA -> CA preselection -> feature extraction -> compact tree inference ->
 * score filtering -> output TrackSoA compaction.
 *
 * Tree inference: the model is a gradient-boosted tree (XGBoost) exported to a compact binary
 * (int8 split-feature, fp32 threshold/leaf value, int32 children, per-tree root offsets, base
 * logit). It is loaded once per process into a shared GlobalCache (DispTreeCache), lifted to the
 * device once as a read-only buffer and traversed by a hand-written alpaka kernel, so the model
 * (O(10 MB)) lives once per device, not once per stream.
 */

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include <cstdint>
#include <fstream>
#include <optional>
#include <utility>

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaCore/interface/CopyToDeviceCache.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackSoA/interface/TrackDefinitions.h"
#include "DataFormats/TrackSoA/interface/TracksDevice.h"
#include "DataFormats/TrackSoA/interface/TracksHost.h"
#include "DataFormats/TrackSoA/interface/alpaka/TracksSoACollection.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/OTRecHitsSoACollection.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"

#include "RecoTracker/FinalTrackSelectors/interface/PixelTrackFeaturesSoA.h"
#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackFeaturesDeviceCollection.h"
#include "RecoTracker/FinalTrackSelectors/plugins/alpaka/PixelTrackTorchHighPuritySelectorKernels.h"

// #define PIXEL_TRACK_HP_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // Per-process shared compact gradient-boosted tree. Loaded once per process (from the FileInPath
  // binary) into host buffers; cms::alpakatools::CopyToDeviceCache then lifts a read-only copy to
  // every visible device at construction, and get(queue) returns the buffer resident on that queue's
  // device. Read-only, so concurrency-safe with no mutex. Per-device is required: the framework
  // spreads streams across all visible GPUs, so a single-device buffer would be an illegal access
  // from streams running on the other device(s).
  // Binary format: int32 nNodes, int32 nTrees, float baseLogit, then int8 feat[N] (-1=leaf),
  // float val[N] (threshold / leaf value), int32 left[N], int32 right[N], int32 roots[nTrees].
  // feat[] is a position in the feature vector the scoring kernel assembles, so load() range-checks
  // it against kNPixelTrackFeatures: the format has no version field, and a stale/foreign ABI would
  // otherwise show up only as an out-of-bounds device read and garbage scores.
  template <typename T>
  using TreeArrayCache = cms::alpakatools::CopyToDeviceCache<Device, cms::alpakatools::host_buffer<T[]>>;

  struct DispTreeCache {
    explicit DispTreeCache(std::string const& path) : DispTreeCache(load(path)) {}

    int nNodes = 0, nTrees = 0;
    float baseLogit = 0.f;
    TreeArrayCache<int8_t> feat;
    TreeArrayCache<float> val;
    TreeArrayCache<int32_t> left, right, roots;

  private:
    // Host-side staging read from the binary, consumed by the delegating constructor below.
    struct HostTree {
      int nNodes, nTrees;
      float baseLogit;
      cms::alpakatools::host_buffer<int8_t[]> feat;
      cms::alpakatools::host_buffer<float[]> val;
      cms::alpakatools::host_buffer<int32_t[]> left, right, roots;
    };
    static HostTree load(std::string const& path) {
      std::ifstream in(path, std::ios::binary);
      if (!in)
        throw cms::Exception("PixelTrackConfiguration") << "cannot open compact tree binary: " << path;
      int32_t nN = 0, nt = 0;
      float bl = 0.f;
      in.read(reinterpret_cast<char*>(&nN), 4);
      in.read(reinterpret_cast<char*>(&nt), 4);
      in.read(reinterpret_cast<char*>(&bl), 4);
      auto feat = cms::alpakatools::make_host_buffer<int8_t[]>(nN);
      auto val = cms::alpakatools::make_host_buffer<float[]>(nN);
      auto left = cms::alpakatools::make_host_buffer<int32_t[]>(nN);
      auto right = cms::alpakatools::make_host_buffer<int32_t[]>(nN);
      auto roots = cms::alpakatools::make_host_buffer<int32_t[]>(nt);
      in.read(reinterpret_cast<char*>(feat.data()), std::streamsize(nN) * sizeof(int8_t));
      in.read(reinterpret_cast<char*>(val.data()), std::streamsize(nN) * sizeof(float));
      in.read(reinterpret_cast<char*>(left.data()), std::streamsize(nN) * sizeof(int32_t));
      in.read(reinterpret_cast<char*>(right.data()), std::streamsize(nN) * sizeof(int32_t));
      in.read(reinterpret_cast<char*>(roots.data()), std::streamsize(nt) * sizeof(int32_t));
      if (!in)
        throw cms::Exception("PixelTrackConfiguration") << "compact tree binary truncated/corrupt: " << path;
      // Feature-index range check. The format has no version field, so the split indices are the
      // only thing that can betray a stale ABI. The scorer reads f[0..kNPixelTrackFeatures-1]
      // with no bound check of its own, so an index >= width is an out-of-bounds device read;
      // a value below -1 is silently mistaken for a leaf (traversal tests `>= 0`). Reject both
      // here. Models narrower than the current width stay valid (they never index the tail).
      for (int n = 0; n < nN; ++n) {
        const int fi = feat[n];
        if (fi >= kNPixelTrackFeatures || fi < -1)
          throw cms::Exception("PixelTrackConfiguration")
              << "compact tree binary " << path << ": node " << n << " splits on feature index " << fi
              << ", outside the [0, " << kNPixelTrackFeatures - 1 << "] range the selector provides (-1 = leaf). "
              << "The model was exported against a different feature ABI than PixelTrackFeaturesSoA's "
              << kNPixelTrackFeatures << " columns.";
      }
      return HostTree{nN, nt, bl, std::move(feat), std::move(val), std::move(left), std::move(right), std::move(roots)};
    }
    explicit DispTreeCache(HostTree h)
        : nNodes(h.nNodes),
          nTrees(h.nTrees),
          baseLogit(h.baseLogit),
          feat(h.feat),
          val(h.val),
          left(h.left),
          right(h.right),
          roots(h.roots) {}
  };

  // SynchronizingEDProducer (two-phase acquire/produce) rather than FixedQueueEDProducer: the
  // latter allocates one Queue per stream and refuses ExternalWork, so a two-phase module
  // cannot keep the fixed-queue property. All of one event's work still runs on one queue (the
  // EDMetadata/queue created in acquire is carried into produce), the shared tree is uploaded
  // once per device at construction and read-only afterwards, and per-event memory is
  // queue-aware via the caching allocators.
  class PixelTrackForestHighPuritySelector : public stream::SynchronizingEDProducer<edm::GlobalCache<DispTreeCache>> {
    using TkSoADevice = reco::TracksSoACollection;
    using HitsOnDevice = reco::TrackingRecHitsSoACollection;
    using OTHitsOnDevice = reco::OTRecHitsSoACollection;
    using TrackHitSoA = ::reco::TrackHitSoA;

  public:
    explicit PixelTrackForestHighPuritySelector(const edm::ParameterSet&, const DispTreeCache*);
    static std::unique_ptr<DispTreeCache> initializeGlobalCache(const edm::ParameterSet& iConfig) {
      return std::make_unique<DispTreeCache>(iConfig.getParameter<edm::FileInPath>("model").fullPath());
    }
    static void globalEndJob(DispTreeCache*) {}
    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    // acquire() runs the whole selection chain asynchronously and issues one async D2H of the
    // two counts the compaction kernel commits. produce() reads the landed counts (plain host
    // memory by then -- no wait), emplaces the collection and publishes the counts as host products.
    void acquire(device::Event const&, device::EventSetup const&) override;
    void produce(device::Event&, device::EventSetup const&) override;

    const device::EDGetToken<TkSoADevice> pixelTrackToken_;
    const int maxNumberOfTracks_;
    const int maxPreselectedTracks_;
    const int minNumberOfHits_;
    const int avgHitsPerTrack_;
    const pixelTrack::Quality minimumTrackQuality_;
    const double scoreThreshold_;
    // dxy-aware threshold ramp: threshold goes from scoreThresholdLowDxy at |dxyBS|=0 down to
    // scoreThreshold at |dxyBS|>=dxyRampKnee. scoreThresholdLowDxy < 0 disables it (flat threshold).
    const double scoreThresholdLowDxy_;
    const double dxyRampKnee_;
    // Hit/stub features: when true, the 10 CA hit/stub features (caTrackFeatures::fill) and the
    // hit-topology, radial-extent, provenance and pixel-cluster columns are appended to the 17
    // fit/cov features. This module requires true (the constructor rejects false); the merged-hits
    // product is consumed only when enabled.
    const bool useHitFeatures_;
    device::EDGetToken<HitsOnDevice> mergedHitsToken_;
    // Raw OT-rechit SoA, so the hit-feature walk can resolve the bit30-tagged raw-OT ids a track
    // may carry. Consumed only when useHitFeatures_ (same guard as mergedHitsToken_).
    device::EDGetToken<OTHitsOnDevice> otRecHitsSoAToken_;
    const device::EDPutToken<TkSoADevice> tokenTrackOut_;
    // Host products carrying the counts this module actually selected, so a consumer can size
    // against the counts instead of against this collection's capacity (both exact by
    // construction -- see PixelTrackFilterKernel).
    const edm::EDPutTokenT<uint32_t> tokenNTracksOut_;
    const edm::EDPutTokenT<uint32_t> tokenNKeptHitsOut_;

    // Seam state (per stream), cleared at acquire entry so a skipped produce() can never leak:
    //  - pendingTracks_: compacted output collection, built in acquire, emplaced in produce.
    //  - countsHost_: 2-word host mirror ([0] tracks, [1] hits); async copy issued at end of
    //    acquire, landed by the time produce runs.
    std::optional<TkSoADevice> pendingTracks_;
    std::optional<cms::alpakatools::host_buffer<uint32_t[]>> countsHost_;
  };

  PixelTrackForestHighPuritySelector::PixelTrackForestHighPuritySelector(const edm::ParameterSet& iConfig,
                                                                         const DispTreeCache*)
      : SynchronizingEDProducer(iConfig),
        pixelTrackToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelTrackSrc"))),
        maxNumberOfTracks_(iConfig.getParameter<int>("maxNumberOfTracks")),
        maxPreselectedTracks_(iConfig.getParameter<int>("maxPreselectedTracks")),
        minNumberOfHits_(iConfig.getParameter<int>("minNumberOfHits")),
        avgHitsPerTrack_(iConfig.getParameter<int>("avgHitsPerTrack")),
        minimumTrackQuality_(pixelTrack::qualityByName(iConfig.getParameter<std::string>("minimumTrackQuality"))),
        scoreThreshold_(iConfig.getParameter<double>("scoreThreshold")),
        scoreThresholdLowDxy_(iConfig.getParameter<double>("scoreThresholdLowDxy")),
        dxyRampKnee_(iConfig.getParameter<double>("dxyRampKnee")),
        useHitFeatures_(iConfig.getParameter<bool>("useHitFeatures")),
        tokenTrackOut_(produces()),
        tokenNTracksOut_(produces("nTracks")),
        tokenNKeptHitsOut_(produces("nKeptHits")) {
    if (useHitFeatures_) {
      mergedHitsToken_ = consumes(iConfig.getParameter<edm::InputTag>("mergedHitsSrc"));
      otRecHitsSoAToken_ = consumes(iConfig.getParameter<edm::InputTag>("otRecHitsSoASrc"));
    }
    if (minimumTrackQuality_ == pixelTrack::Quality::notQuality) {
      throw cms::Exception("PixelTrackConfiguration")
          << iConfig.getParameter<std::string>("minimumTrackQuality") + " is not a pixelTrack::Quality";
    }
    if (minimumTrackQuality_ < pixelTrack::Quality::dup) {
      throw cms::Exception("PixelTrackConfiguration")
          << iConfig.getParameter<std::string>("minimumTrackQuality") + " not supported";
    }
    if (maxPreselectedTracks_ > maxNumberOfTracks_) {
      throw cms::Exception("PixelTrackConfiguration") << "maxPreselectedTracks must be <= maxNumberOfTracks";
    }
    // The tree scorer reads all columns unconditionally; columns 18-42 are written only under
    // useHitFeatures, so refuse useHitFeatures=false (would score uninitialised memory).
    if (!useHitFeatures_) {
      throw cms::Exception("PixelTrackConfiguration")
          << "PixelTrackForestHighPuritySelector requires useHitFeatures = True (the tree scorer "
             "reads the full feature vector, whose columns 18 and above are hit-derived)";
    }
  }

  void PixelTrackForestHighPuritySelector::acquire(device::Event const& iEvent, device::EventSetup const&) {
    /*
    Processing steps:
      1. CA-based preselection of tracks
      2. Feature extraction (track + hit SoA)
      3. Compact tree inference
      4. Score-based filtering
      5. Track compaction and output production
      6. Async readback of the two counts the compaction committed
*/
    // Reset seam state so a skipped produce() can never leak into the next event.
    pendingTracks_.reset();
    countsHost_.reset();

    // Retrieve tokens
    auto& queue = iEvent.queue();
    const auto& tracks = iEvent.get(pixelTrackToken_).view();

    // Instantiate the necessary objects in memory
    //  - Temporary storage for filtering
    auto d_nPreselectedTracks = cms::alpakatools::make_device_buffer<int>(queue);
    auto d_nSelectedTracks = cms::alpakatools::make_device_buffer<int>(queue);
    auto d_preselectedTrackIndices = cms::alpakatools::make_device_buffer<int[]>(queue, maxNumberOfTracks_);
    auto d_selectedTrackIndices = cms::alpakatools::make_device_buffer<int[]>(queue, maxPreselectedTracks_);
    auto d_trackHitCounts = cms::alpakatools::make_device_buffer<int[]>(queue, maxPreselectedTracks_);
    auto d_selectedTrackHitOffsets = cms::alpakatools::make_device_buffer<int[]>(queue, maxPreselectedTracks_);
    auto d_preselectionOffsets = cms::alpakatools::make_device_buffer<int[]>(queue, maxNumberOfTracks_);

    // NO pre-fill of the seven buffers above: each is written before read over a range containing
    // every index any consumer reads. Buffers whose consumer CAN read an unwritten element keep
    // their fill (preselectionMask, selectionMask, selectedTrackHitCounts: the scans sweep the whole
    // capacity while the producing kernels stop at the track count).

    //  - Features and scores containers
    PixelTrackFeaturesOnDevice trackFeatures(queue, maxPreselectedTracks_);
    PixelTrackScoresOnDevice trackScoresOnDevice(queue, maxPreselectedTracks_);

    // Optional debug definitions
#ifdef PIXEL_TRACK_HP_DEBUG
    auto h_nPreselectedTracks = cms::alpakatools::make_host_buffer<int>(queue);
    auto h_nSelectedTracks = cms::alpakatools::make_host_buffer<int>(queue);
    auto nPreselectedTracks = 0;
    auto nSelectedTracks = 0;
    // Helper to copy the number of kept tracks back to host (debug only)
    auto fetchnPreselectedTracks = [&]() {
      alpaka::memcpy(queue, h_nPreselectedTracks, d_nPreselectedTracks);
      alpaka::wait(queue);
      return *h_nPreselectedTracks;
    };
    auto fetchnSelectedTracks = [&]() {
      alpaka::memcpy(queue, h_nSelectedTracks, d_nSelectedTracks);
      alpaka::wait(queue);
      return *h_nSelectedTracks;
    };
#endif

    // 1. CA-based preselection of tracks
    //  Launch first kernel to look which tracks need to be filtered out
    //  based on quality criteria from the CA
    launchCAPreselection(queue,
                         maxNumberOfTracks_,
                         minNumberOfHits_,
                         minimumTrackQuality_,
                         tracks.tracks(),
                         alpaka::getPtrNative(d_preselectedTrackIndices),
                         alpaka::getPtrNative(d_preselectionOffsets),
                         alpaka::getPtrNative(d_nPreselectedTracks));

#ifdef PIXEL_TRACK_HP_DEBUG
    nPreselectedTracks = fetchnPreselectedTracks();
    std::cout << "PixelTrackForestHighPuritySelector::Prefiltered tracks=" << nPreselectedTracks << "\n";
#endif

    // 2. Feature extraction (track + hit SoA)
    //  The merged TrackingRecHitsSoA (same product the CA indexed: its trackHits().id()
    //  point into it) is needed ONLY for the hit/stub features. Fetch its view when enabled;
    //  otherwise pass a null view + nHitsTot=0 (the kernel reads `hits` only under useHitFeatures).
    ::reco::TrackingRecHitConstView mergedHitsView{};
    int nHitsTot = 0;
    // OT-rechit view for resolving tagged OT extras (empty view + 0 when not needed).
    ::reco::OTRecHitsConstView otHitsView{};
    uint32_t nOTHits = 0;
    if (useHitFeatures_) {
      const auto& mergedHits = iEvent.get(mergedHitsToken_);
      mergedHitsView = mergedHits.view().trackingHits();
      nHitsTot = mergedHitsView.metadata().size();
      const auto& otHits = iEvent.get(otRecHitsSoAToken_);
      otHitsView = otHits.const_view().otRecHits();
      nOTHits = otHitsView.metadata().size();
    }
    launchFeaturesExtractor(queue,
                            maxPreselectedTracks_,
                            tracks.tracks(),
                            tracks.trackHits(),
                            mergedHitsView,
                            nHitsTot,
                            otHitsView,
                            nOTHits,
                            useHitFeatures_,
                            alpaka::getPtrNative(d_preselectedTrackIndices),
                            alpaka::getPtrNative(d_nPreselectedTracks),
                            trackFeatures.view(),
                            alpaka::getPtrNative(d_trackHitCounts));

    // 3. Tree inference: a custom traversal of the compact gradient-boosted tree, reading the buffers
    // resident on THIS queue's device from the per-device shared cache -> no per-stream model copy.
    auto const& tc = *globalCache();
    launchTreeScore(queue,
                    maxPreselectedTracks_,
                    alpaka::getPtrNative(tc.feat.get(queue)),
                    alpaka::getPtrNative(tc.val.get(queue)),
                    alpaka::getPtrNative(tc.left.get(queue)),
                    alpaka::getPtrNative(tc.right.get(queue)),
                    alpaka::getPtrNative(tc.roots.get(queue)),
                    tc.nTrees,
                    tc.baseLogit,
                    trackFeatures.const_view(),
                    alpaka::getPtrNative(d_nPreselectedTracks),
                    trackScoresOnDevice.view());

    // 4. Score-based filtering
    launchScoreFilter(queue,
                      maxPreselectedTracks_,
                      scoreThreshold_,
                      scoreThresholdLowDxy_,
                      dxyRampKnee_,
                      trackFeatures.const_view(),
                      trackScoresOnDevice.view(),
                      alpaka::getPtrNative(d_preselectedTrackIndices),
                      alpaka::getPtrNative(d_nPreselectedTracks),
                      alpaka::getPtrNative(d_trackHitCounts),
                      alpaka::getPtrNative(d_selectedTrackIndices),
                      alpaka::getPtrNative(d_nSelectedTracks),
                      alpaka::getPtrNative(d_selectedTrackHitOffsets));

#ifdef PIXEL_TRACK_HP_DEBUG
    nSelectedTracks = fetchnSelectedTracks();
    std::cout << "PixelTrackForestHighPuritySelector::Filtered tracks=" << nSelectedTracks << "\n";
#endif
    // 2-word device scratch the compaction kernel fills with the counts it commits. No pre-fill:
    // PixelTrackFilterKernel assigns both words unconditionally under once_per_grid before the copy.
    // Function scope is safe (caching allocator re-hands a freed block only after the queue's free event completes).
    auto d_selectedCounts = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, 2u);
    auto tracks_out = launchProduceOutputTracks(queue,
                                                maxPreselectedTracks_,
                                                avgHitsPerTrack_,
                                                tracks.tracks(),
                                                tracks.trackHits(),
                                                alpaka::getPtrNative(d_selectedTrackIndices),
                                                alpaka::getPtrNative(d_nSelectedTracks),
                                                alpaka::getPtrNative(d_selectedTrackHitOffsets),
                                                alpaka::getPtrNative(d_selectedCounts));
    // Async 8-byte D2H into a pinned buffer; no wait here (the framework's between-phase drain delivers it).
    countsHost_.emplace(cms::alpakatools::make_host_buffer<uint32_t[]>(queue, 2u));
    alpaka::memcpy(queue, *countsHost_, d_selectedCounts);
    pendingTracks_.emplace(std::move(tracks_out));
  }

  void PixelTrackForestHighPuritySelector::produce(device::Event& iEvent, device::EventSetup const&) {
    // The counts landed while the framework waited on this event's queue: plain host memory here,
    // read with no synchronisation of any kind.
    uint32_t nTracksSel = 0;
    uint32_t nKeptHitsSel = 0;
    if (countsHost_) {
      nTracksSel = countsHost_->data()[0];
      nKeptHitsSel = countsHost_->data()[1];
    }

    if (pendingTracks_) {
      iEvent.emplace(tokenTrackOut_, std::move(*pendingTracks_));
    } else {
      // acquire() has no early-out today, so this cannot trigger; it keeps a future early-out from
      // dereferencing an empty optional and guarantees all three products always exist.
      auto& queue = iEvent.queue();
      TkSoADevice empty(queue, 0, 0);
      auto nTracks_d = cms::alpakatools::make_device_view(queue, empty.view().tracks().nTracks());
      alpaka::memset(queue, nTracks_d, 0);
      iEvent.emplace(tokenTrackOut_, std::move(empty));
    }
    iEvent.emplace(tokenNTracksOut_, nTracksSel);
    iEvent.emplace(tokenNKeptHitsOut_, nKeptHitsSel);

    pendingTracks_.reset();
    countsHost_.reset();
  }

  void PixelTrackForestHighPuritySelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("pixelTrackSrc", {"hltPhase2PixelTracksSoA"});
    desc.add<int>("maxNumberOfTracks", 100000);
    desc.add<int>("maxPreselectedTracks", 10000);
    desc.add<int>("minNumberOfHits", 0);
    desc.add<int>("avgHitsPerTrack", 8);
    desc.add<std::string>("minimumTrackQuality", "tight");
    // Compact gradient-boosted tree binary (NOT a TorchScript .pt): int32 nNodes, int32 nTrees,
    // float baseLogit, then int8 split-feature[N] (-1=leaf), float threshold/leaf[N],
    // int32 left[N]/right[N], int32 roots[nTrees]. Split-feature indices are positions in
    // the PixelTrackFeaturesSoA column order; the loader throws if any is >= kNPixelTrackFeatures.
    desc.add<edm::FileInPath>("model");
    desc.add<double>("scoreThreshold", 0.5);
    // dxy-aware threshold ramp: scoreThresholdLowDxy<0 disables it (flat scoreThreshold). When
    // >=0, the cut ramps from scoreThresholdLowDxy at |dxyBS|=0 to scoreThreshold at |dxyBS|>=dxyRampKnee[cm]
    // -- a more aggressive cut at low reco displacement (where the OT iteration mostly makes fakes).
    desc.add<double>("scoreThresholdLowDxy", -1.0);
    desc.add<double>("dxyRampKnee", 2.0);
    // Hit/stub features. This module requires useHitFeatures = True (the tree scorer reads the full
    // feature vector); mergedHitsSrc is the merged TrackingRecHitsSoA the CA indexed (its
    // trackHits().id() point into it) and is only consumed when useHitFeatures=true.
    desc.add<bool>("useHitFeatures", true);
    desc.add<edm::InputTag>("mergedHitsSrc", {"hltPhase2PixelRecHitsStubsMerger"});
    // Raw OT-rechit SoA (only consumed when useHitFeatures=true): resolves the bit30-tagged
    // raw-OT hit ids so OT-extended tracks are scored on their full hit content.
    desc.add<edm::InputTag>("otRecHitsSoASrc", {"hltPixelSeedingOTRecHitsSoA"});
    descriptions.addWithDefaultLabel(desc);
  }
};  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(PixelTrackForestHighPuritySelector);
