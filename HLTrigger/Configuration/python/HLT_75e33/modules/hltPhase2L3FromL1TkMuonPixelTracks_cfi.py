import FWCore.ParameterSet.Config as cms

hltPhase2L3FromL1TkMuonPixelTracks = cms.EDProducer("PixelTrackProducer",
    Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
    Filter = cms.InputTag("hltPhase2PixelTrackFilterByKinematics"),
    Fitter = cms.InputTag("hltPhase2PixelFitterByHelixProjections"),
    SeedingHitSets = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets"),
    passLabel = cms.string('')
)

from RecoMuon.L3TrackFinder.MuonTracksSelectorFromL1TkMuon import MuonTracksSelectorFromL1TkMuon as _MuonTracksSelectorFromL1TkMuon
_hltPhase2L3FromL1TkMuonPixelTracks = _MuonTracksSelectorFromL1TkMuon(
    TrackInputCollection = "hltPhase2PixelTracksCAExtension",
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
(phase2MuonPixelTracksSelector & phase2CAExtension).toReplaceWith(
    hltPhase2L3FromL1TkMuonPixelTracks,
    _hltPhase2L3FromL1TkMuonPixelTracks
)
