import FWCore.ParameterSet.Config as cms

# Smoke test of the four muon HP forest selectors: constructs every module
# (parameter-set validation, plugin registration) and loads + validates its
# compact forest against the feature extractor (muonhp::CompactForest::load).
# No events are processed: the selectors need HLT tracking products. The
# per-track behaviour is validated on real events by the training-side
# cross-check (muonHighPurityTrackSelection/production/features_validation).

process = cms.Process("FORESTSMOKE")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(0))
process.MessageLogger.cerr.MuonIOTracksForestSelector = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.MessageLogger.cerr.MuonOITracksForestSelector = cms.untracked.PSet(limit=cms.untracked.int32(-1))

process.load("HLTrigger.Configuration.HLT_75e33.modules.hltPhase2MuonPixelTracksHighPurityForest_cfi")
process.load("HLTrigger.Configuration.HLT_75e33.modules.hltPhase2MuonIOTrackSelectionHighPurityForest_cfi")
from HLTrigger.Configuration.HLT_75e33.modules.hltPhase2L3OIMuonTrackSelectionHighPurityForest_cfi import (
    _pixelOIForestSelector,
    _seedsOIForestSelector,
)

process.pixelOIForestSelector = _pixelOIForestSelector.clone()
process.seedsOIForestSelector = _seedsOIForestSelector.clone()

process.p = cms.Path(
    process.hltPhase2MuonPixelTracksHighPurityForest
    + process.hltPhase2MuonIOTrackSelectionHighPurityForest
    + process.pixelOIForestSelector
    + process.seedsOIForestSelector
)
