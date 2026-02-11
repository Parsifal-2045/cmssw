import FWCore.ParameterSet.Config as cms

hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracks"),
    InputVertexCollection = cms.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices"),
    SeedCreatorPSet = cms.PSet(
        refToPSet_ = cms.string('hltPhase2SeedFromProtoTracks')
    ),
    TTRHBuilder = cms.string('WithTrackAngle'),
    originHalfLength = cms.double(0.3),
    originRadius = cms.double(0.1),
    useEventsWithNoVertex = cms.bool(True),
    usePV = cms.bool(False),
    useProtoTrackKinematics = cms.bool(False)
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
(phase2MuonPixelTracksSelector & phase2CAExtension).toModify(
    hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks,
    InputCollection=cms.InputTag("hltPhase2MuonPixelTracksDNNSelector"),
    InputVertexCollection=cms.InputTag(""),
)
