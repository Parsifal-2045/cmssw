import FWCore.ParameterSet.Config as cms

process = cms.Process("FORESTSMOKE")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

process.load("HLTrigger.Configuration.HLT_75e33.modules.hltPhase2MuonPixelTracksHighPurityForest_cfi")
process.load("HLTrigger.Configuration.HLT_75e33.modules.hltPhase2MuonIOTrackSelectionHighPurityForest_cfi")
process.load("HLTrigger.Configuration.HLT_75e33.modules.hltPhase2L3OIMuonTrackSelectionHighPurityForest_cfi")

from HLTrigger.Configuration.HLT_75e33.modules.hltPhase2L3OIMuonTrackSelectionHighPurityForest_cfi import _pixelOIForestSelector, _seedsOIForestSelector

process.pixelOIForestSelector = _pixelOIForestSelector.clone()
process.generalOIForestSelector = _seedsOIForestSelector.clone()

process.p = cms.Path(
    process.hltPhase2MuonPixelTracksHighPurityForest
    + process.hltPhase2MuonIOTrackSelectionHighPurityForest
    + process.pixelOIForestSelector
    + process.generalOIForestSelector
)
