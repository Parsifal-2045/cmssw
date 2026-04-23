import FWCore.ParameterSet.Config as cms

hltPhase2L3OIMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltPhase2L3OIMuCtfWithMaterialTracks")
)

_pixelDNNSelector = cms.EDProducer('MuonOITracksDNNSelector',
    tracks = cms.InputTag('hltPhase2L3OIMuCtfWithMaterialTracks'),
    standaloneMuons = cms.InputTag('hltL2MuonsFromL1TkMuon', 'UpdatedAtVtx'),
    modelPath = cms.string('RecoMuon/L3TrackFinder/data/OI_pixel_model.onnx'),
    decisionThreshold = cms.double(0.9775444269180298),
    useStandaloneMuonFeatures = cms.bool(True),
    nFeatures = cms.int32(26),
)

_seedsDNNSelector = _pixelDNNSelector.clone(
    decisionThreshold = cms.double(0.9548693895339966),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/OI_general_model.onnx")
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.ngtScouting_cff import ngtScouting
(phase2MuonPixelTracksSelector | ngtScouting).toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _pixelDNNSelector
)

from Configuration.ProcessModifiers.phase2MuonSeedsSelector_cff import phase2MuonSeedsSelector
(phase2MuonSeedsSelector).toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _seedsDNNSelector
)
