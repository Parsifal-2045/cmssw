import FWCore.ParameterSet.Config as cms

hltPhase2PixelTrackTorchHighPuritySelector = cms.EDProducer('PixelTrackTorchHighPuritySelector@alpaka',
    pixelTrackSrc = cms.InputTag('hltPhase2PixelTracksSoA'),
    maxNumberOfTracks = cms.int32(2*60*1024),
    maxPreselectedTracks = cms.int32(9_984),
    minNumberOfHits = cms.int32(0),
    avgHitsPerTrack = cms.int32(8),
    minimumTrackQuality = cms.string('tight'),
    model = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/PixelTrackTorchHighPuritySelector/pixel_track_classifier_FP16.pt'),
    scoreThreshold = cms.double(0.4),
    batchSize = cms.int32(4_992)
)

# Stub chain (phase2CAStubs): replaces the Torch MLP with a 31-feature gradient-boosted
# forest (17 fit/covariance + 10 hit/stub CA features + 4 hit-walk features). toReplaceWith keeps the module label,
# so all downstream wiring and the SerialSync clone stay untouched.
_hltPhase2PixelTrackForestHighPuritySelector = cms.EDProducer('PixelTrackForestHighPuritySelector@alpaka',
    pixelTrackSrc = cms.InputTag('hltPhase2PixelTracksSoA'),
    maxNumberOfTracks = cms.int32(2*60*1024),
    maxPreselectedTracks = cms.int32(9_984),     # sized for PU200 with headroom, so nothing is truncated
    minNumberOfHits = cms.int32(0),
    avgHitsPerTrack = cms.int32(8),
    minimumTrackQuality = cms.string('tight'),
    useHitFeatures = cms.bool(True),             # build the full 31-feature vector
    mergedHitsSrc = cms.InputTag('hltPhase2PixelRecHitsStubsMerger'),
    # Compact gradient-boosted forest binary, shared per device.
    model = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/PixelTrackTorchHighPuritySelector/prompt_tree31_ew4m0_20260829.bin'),
    # Flat working point on this arm: the dxy-dependent threshold ramp is disabled, so the
    # same score cut applies at every impact parameter.
    scoreThreshold = cms.double(0.185),
    scoreThresholdLowDxy = cms.double(-1.0),
    dxyRampKnee = cms.double(3.0)
)

from Configuration.ProcessModifiers.phase2CAStubs_cff import phase2CAStubs
phase2CAStubs.toReplaceWith(hltPhase2PixelTrackTorchHighPuritySelector,
                            _hltPhase2PixelTrackForestHighPuritySelector)
