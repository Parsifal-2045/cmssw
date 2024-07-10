import FWCore.ParameterSet.Config as cms

L3IOFIRST = True

phase2L2MuonTracksToReuse = cms.EDProducer('phase2HLTMuonSelectorForL3',
 l2Muons = cms.InputTag("hltL2MuonsFromL1TkMuon"),
 l2MuonsUpdVtx = cms.InputTag("hltL2MuonsFromL1TkMuon", "UpdatedAtVtx"),
 TrajToL2AssociationMap = cms.InputTag("hltL2MuonsFromL1TkMuon"),
 l3Tracks = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMerged") if L3IOFIRST else cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
 #l3IOTracks = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMerged"),
 #l3OITracks = cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
 IOFirst = cms.bool(L3IOFIRST),
 matchingDr = cms.double(0.02),
 applyL3Filters = cms.bool(True),
 MinNhits = cms.int32(1),
 MaxNormalizedChi2 = cms.double(20.0),
 MinNmuonHits = cms.int32(0),
 MaxPtDifference = cms.double(999.0)
)


# hltPhase2L3OIL3MuonCandidates L3 muonCandidates OI

# Use tracks for both OI / IO for now
# Match in dR to see which seeds to reuse 

# OI tracks (high purity selection)
# hltPhase2L3OIMuonTrackSelectionHighPurity 

# IO tracks (high purity selection)
# hltIter2Phase2L3FromL1TkMuonMerged 