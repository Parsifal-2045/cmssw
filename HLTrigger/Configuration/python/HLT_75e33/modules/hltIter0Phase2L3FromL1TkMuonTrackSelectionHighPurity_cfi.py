import FWCore.ParameterSet.Config as cms

hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks")
)

_hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity = cms.EDProducer("MuonIOTracksDNNSelector",
    tracks = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks"),
    decisionThreshold = cms.double(0.9456754326820374),
    nFeatures = cms.int32(36),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/general_track_selector.onnx")
)

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
from Configuration.ProcessModifiers.trackingLST_cff import trackingLST
(phase2CAExtension & trackingLST).toReplaceWith(
    hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity,
    _hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity
)
