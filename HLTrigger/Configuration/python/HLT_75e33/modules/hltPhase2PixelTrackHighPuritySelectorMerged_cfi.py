import FWCore.ParameterSet.Config as cms

# Final high-purity selection of the two-iteration stub chain (phase2CAStubs & pixelTrackMask):
# it scores the merger output, i.e. a merged, OT-extended, refit and de-duplicated collection
# that neither of the per-iteration selectors has seen.
hltPhase2PixelTrackHighPuritySelectorMerged = cms.EDProducer('PixelTrackForestHighPuritySelector@alpaka',
    pixelTrackSrc = cms.InputTag('hltPhase2PixelTracksSoAMerger'),

    # Caps on the tracks kept and scored, sized for PU200 with ample headroom.
    maxNumberOfTracks = cms.int32(2*60*1024),
    maxPreselectedTracks = cms.int32(9_984),
    minNumberOfHits = cms.int32(0),
    # The merger's attach walk adds OT hits on top of the CA hit content.
    avgHitsPerTrack = cms.int32(16),
    minimumTrackQuality = cms.string('tight'),

    # Merged hit collection plus the raw OT-rechit SoA, so that bit30-tagged raw-OT ids resolve.
    useHitFeatures = cms.bool(True),
    mergedHitsSrc = cms.InputTag('hltPhase2PixelRecHitsStubsMerger'),
    otRecHitsSoASrc = cms.InputTag('hltPixelSeedingOTRecHitsSoA'),

    # Forest trained on the merged collection; the per-iteration models do not apply here.
    model = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/PixelTrackTorchHighPuritySelector/'
                           'merged_tree42_finalhpfrw_20260823.bin'),
    # dxy-aware working point: the cut ramps from scoreThresholdLowDxy at |dxyBS| = 0 to
    # scoreThreshold at |dxyBS| >= dxyRampKnee, and is tuned for high recall on primary matches.
    scoreThreshold = cms.double(0.02),          # displaced arm: loose, so no vertex-position bin loses efficiency
    scoreThresholdLowDxy = cms.double(0.0422),  # prompt-like arm
    dxyRampKnee = cms.double(1.0)               # cm
)
