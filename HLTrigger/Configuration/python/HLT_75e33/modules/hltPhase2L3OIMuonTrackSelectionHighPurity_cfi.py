import FWCore.ParameterSet.Config as cms

hltPhase2L3OIMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltPhase2L3OIMuCtfWithMaterialTracks")
)

_hltPhase2L3OIMuonTrackSelectionHighPurityPixelSelector = cms.EDProducer("MuonOITracksDNNSelector",
    decisionThreshold = cms.double(0.9775444269180298),
    nFeatures = cms.int32(26),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/OI_pixel_model.onnx")
)

_hltPhase2L3OIMuonTrackSelectionHighPurityGeneralSelector = cms.EDProducer("MuonOITracksDNNSelector",
    decisionThreshold = cms.double(0.9548693895339966),
    nFeatures = cms.int32(26),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/OI_general_model.onnx")
)

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
phase2CAExtension.toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _hltPhase2L3OIMuonTrackSelectionHighPurityPixelSelector
)

from Configuration.ProcessModifiers.trackingLST_cff import trackingLST
(phase2CAExtension & trackingLST).toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _hltPhase2L3OIMuonTrackSelectionHighPurityGeneralSelector
)
