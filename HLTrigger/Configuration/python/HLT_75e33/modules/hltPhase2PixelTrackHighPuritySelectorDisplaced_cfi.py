import FWCore.ParameterSet.Config as cms

# Stage-2 high-purity selector of the displaced iteration: a gradient-boosted decision
# forest, as opposed to the Torch MLP used by the prompt arm.
hltPhase2PixelTrackHighPuritySelectorDisplaced = cms.EDProducer('PixelTrackForestHighPuritySelector@alpaka',
    pixelTrackSrc = cms.InputTag('hltPhase2PixelTracksSoADisplaced'),
    maxNumberOfTracks = cms.int32(10_000),
    # Hard cap on the number of tracks that get scored, sized for PU200 with headroom.
    maxPreselectedTracks = cms.int32(3_072),
    minNumberOfHits = cms.int32(0),
    avgHitsPerTrack = cms.int32(8),
    minimumTrackQuality = cms.string('tight'),
    # Build the full 31-feature vector (17 fit/covariance + 10 hit/stub CA features +
    # rzChi2, meanStubKappa, leverArm, rMax); the deployed forest reads features 0..26.
    useHitFeatures = cms.bool(True),
    mergedHitsSrc = cms.InputTag('hltPhase2PixelRecHitsStubsMerger'),
    # Compact custom binary (not TorchScript), loaded once per process into a shared GlobalCache.
    model = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/PixelTrackTorchHighPuritySelector/disp_tree31_wp_20260821.bin'),
    # dxy-aware working point: the cut ramps from scoreThresholdLowDxy at |dxyBS| = 0
    # to scoreThreshold at |dxyBS| >= dxyRampKnee.
    scoreThreshold = cms.double(0.1824),        # high-dxy arm
    scoreThresholdLowDxy = cms.double(-1.0),    # low-dxy (near-prompt) arm
    dxyRampKnee = cms.double(3.0)               # cm
)
