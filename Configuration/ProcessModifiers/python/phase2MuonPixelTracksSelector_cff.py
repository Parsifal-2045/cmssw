import FWCore.ParameterSet.Config as cms

# This modifier, together with phase2CAExtension,
# enables single iteration muon Inside-Out tracking
# by selecting extended PixelTracks matching L1TkMuons,
# fitering them with a DNN and finally extending them with CKF
phase2MuonPixelTracksSelector = cms.Modifier()
