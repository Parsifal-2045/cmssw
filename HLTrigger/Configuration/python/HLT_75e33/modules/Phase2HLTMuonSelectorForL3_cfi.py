import FWCore.ParameterSet.Config as cms

L3IOFIRST = True

phase2L3FilteredObjects = cms.EDProducer('phase2HLTMuonSelectorForL3',
 l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
 l2MuonsUpdVtx = cms.InputTag("hltL2MuonsFromL1TkMuon", "UpdatedAtVtx"),
 l3Tracks = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMerged") if L3IOFIRST else cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
 IOFirst = cms.bool(L3IOFIRST),
 matchingDr = cms.double(0.02),
 applyL3Filters = cms.bool(True),
 MinNhits = cms.int32(1),
 MaxNormalizedChi2 = cms.double(5.0),
 MinNhitsMuons = cms.int32(0),
 MinNhitsPixel = cms.int32(1),
 MinNhitsTracker = cms.int32(6),
 MaxPtDifference = cms.double(999.0)
)
