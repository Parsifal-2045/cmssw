import FWCore.ParameterSet.Config as cms

hltPhase2L3OIMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltPhase2L3OIMuCtfWithMaterialTracks")
)

from RecoMuon.L3TrackFinder.MuonOITracksDNNSelector import MuonOITracksDNNSelector as _MuonOITracksDNNSelector

_hltPhase2L3OIMuonTrackSelectionHighPurityPixelSelector = _MuonOITracksDNNSelector(
    decisionThreshold = cms.double(0.9775444269180298),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/OI_pixel_model.onnx")
)

_hltPhase2L3OIMuonTrackSelectionHighPurityGeneralSelector = _MuonOITracksDNNSelector(
    decisionThreshold = cms.double(0.9548693895339966),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/OI_general_model.onnx")
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
(phase2MuonPixelTracksSelector & phase2CAExtension).toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _hltPhase2L3OIMuonTrackSelectionHighPurityPixelSelector
)

from Configuration.ProcessModifiers.phase2MuonGeneralTracksSelector_cff import phase2MuonGeneralTracksSelector
(phase2MuonGeneralTracksSelector).toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _hltPhase2L3OIMuonTrackSelectionHighPurityGeneralSelector
)
