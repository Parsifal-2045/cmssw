import FWCore.ParameterSet.Config as cms

hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
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

from RecoMuon.L3TrackFinder.MuonTracksSelectorFromL1TkMuon import MuonTracksSelectorFromL1TkMuon as _MuonTracksSelectorFromL1TkMuon
_hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks = _MuonTracksSelectorFromL1TkMuon(
    TrackInputCollection = "hltGeneralTracks",
    trackMaxEta = 3.0,
    maxDz = 1,
    maxDr = 0.1,
    maxChi2 = 9
)

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
from Configuration.ProcessModifiers.trackingLST_cff import trackingLST
(phase2CAExtension & trackingLST).toReplaceWith(
    hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks,
    _hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks
)
