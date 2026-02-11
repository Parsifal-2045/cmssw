import FWCore.ParameterSet.Config as cms

# This modifier enables single iteration muon
# Inside-Out tracking by selecting general tracks
# matching L1TkMuons and fitering them with a DNN.
# The resulting tracks are then used as Inside-Out
# tracker muon tracks
phase2MuonGeneralTracksSelector = cms.Modifier()
