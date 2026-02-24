import FWCore.ParameterSet.Config as cms

hltDisplacedMuonsSiClusters = cms.EDProducer("TrackClusterRemoverPhase2",
    TrackQuality = cms.string(''),
    maxChi2 = cms.double(9.0),
    mightGet = cms.optional.untracked.vstring,
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    phase2OTClusters = cms.InputTag("hltSiPhase2Clusters"),
    phase2pixelClusters = cms.InputTag("hltSiPixelClusters"),
    trackClassifier = cms.InputTag(""),
    trajectories = cms.InputTag("hltPhase2L3MuonMerged")
)
