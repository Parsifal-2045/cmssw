import FWCore.ParameterSet.Config as cms

hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks")
)

_hltSingleIterL3FromL1TkMuonPixelTrackSelection = cms.EDProducer("TrackProducer",
    AlgorithmName = cms.string('hltIter0'),
    Fitter = cms.string('FlexibleKFFittingSmoother'),
    GeometricInnerState = cms.bool(True),
    MeasurementTracker = cms.string(''),
    MeasurementTrackerEvent = cms.InputTag("hltMeasurementTrackerEvent"),
    NavigationSchool = cms.string(''),
    Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
    SimpleMagneticField = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TrajectoryInEvent = cms.bool(False),
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    clusterRemovalInfo = cms.InputTag(""),
    src = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates"),
    useHitsSplitting = cms.bool(False),
    useSimpleMF = cms.bool(False)
)

_hltSingleIterL3FromL1TkMuonGeneralTrackSelectionHighPurity = cms.EDProducer("MuonIOTracksDNNSelector",
    tracks = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks"),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    decisionThreshold = cms.double(0.9456754326820374),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(36),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/general_track_selector.onnx")
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
(phase2MuonPixelTracksSelector & phase2CAExtension).toReplaceWith(
    hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity,
    _hltSingleIterL3FromL1TkMuonPixelTrackSelection
)

from Configuration.ProcessModifiers.phase2MuonGeneralTracksSelector_cff import phase2MuonGeneralTracksSelector
phase2MuonGeneralTracksSelector.toReplaceWith(
    hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity,
    _hltSingleIterL3FromL1TkMuonGeneralTrackSelectionHighPurity
)
