import FWCore.ParameterSet.Config as cms

hltPhase2L3FromL1TkMuonPixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltPhase2PixelTrackFilterByKinematics"),
    Fitter = cms.InputTag("hltPhase2PixelFitterByHelixProjections"),
    SeedingHitSets = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets"),
    passLabel = cms.string('')
)

_hltPhase2L3FromL1TkMuonPixelTracks = cms.EDProducer("TrackSelectorByRegion",
 produceMask = cms.bool(False),
 produceTrackCollection = cms.bool(True),
 regions = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions"),
 tracks = cms.InputTag("hltPhase2PixelTracks")
)

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
phase2CAExtension.toReplaceWith(
    hltPhase2L3FromL1TkMuonPixelTracks,
    _hltPhase2L3FromL1TkMuonPixelTracks
)

