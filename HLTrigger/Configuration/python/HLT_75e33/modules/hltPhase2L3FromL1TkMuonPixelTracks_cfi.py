import FWCore.ParameterSet.Config as cms

hltPhase2L3FromL1TkMuonPixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltPhase2PixelTrackFilterByKinematics"),
    Fitter = cms.InputTag("hltPhase2PixelFitterByHelixProjections"),
    SeedingHitSets = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets"),
    passLabel = cms.string('')
)

#_hltPhase2L3FromL1TkMuonPixelTracks = cms.EDProducer("TrackSelectorByRegion",
#    produceMask = cms.bool(False),
#    produceTrackCollection = cms.bool(True),
#    regions = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions"),
#    tracks = cms.InputTag("hltPhase2PixelTracksCAExtension")
#)

from RecoMuon.L3TrackFinder.MuonTracksSelectorFromL1TkMuon import MuonTracksSelectorFromL1TkMuon as _MuonTracksSelectorFromL1TkMuon
_hltPhase2L3FromL1TkMuonPixelTracks = _MuonTracksSelectorFromL1TkMuon(
    TrackInputCollection = "hltPhase2PixelTracksCAExtension",
    trackMaxEta = 3.0,
    maxDz = 1,
    maxDr = 0.1,
    maxChi2 = 9
)

from ..modules.hltPhase2PixelTracksCAExtension_cfi import hltPhase2PixelTracksCAExtension as _hltPhase2PixelTracksCAExtension

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
phase2CAExtension.toReplaceWith(
    hltPhase2L3FromL1TkMuonPixelTracks,
    #_hltPhase2PixelTracksCAExtension
    _hltPhase2L3FromL1TkMuonPixelTracks
)

