import FWCore.ParameterSet.Config as cms

hltPhase2PixelTrackSoATableProducer = cms.EDProducer("HLTPixelTrackSoATableProducer",
    trackSrc = cms.InputTag("hltPhase2PixelTrackTorchHighPuritySelector"),
)

# Two-iteration stub chain (pixelTrackMask): the table source is the merger's
# (merged + OT-extended + refit + de-duplicated) SoA, not the prompt arm alone.
from Configuration.ProcessModifiers.pixelTrackMask_cff import pixelTrackMask
pixelTrackMask.toModify(hltPhase2PixelTrackSoATableProducer, trackSrc = "hltPhase2PixelTracksSoAMerger")
