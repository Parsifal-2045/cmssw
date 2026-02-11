import FWCore.ParameterSet.Config as cms

hltIter2Phase2L3FromL1TkMuonClustersRefRemoval = cms.EDProducer("TrackClusterRemoverPhase2",
    TrackQuality = cms.string('highPurity'),
    maxChi2 = cms.double(16.0),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    oldClusterRemovalInfo = cms.InputTag(""),
    overrideTrkQuals = cms.InputTag(""),
    phase2OTClusters = cms.InputTag("hltSiPhase2Clusters"),
    phase2pixelClusters = cms.InputTag("hltSiPixelClusters"),
    trackClassifier = cms.InputTag("","QualityMasks"),
    trajectories = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity")
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
(phase2MuonPixelTracksSelector & phase2CAExtension).toModify(
    hltIter2Phase2L3FromL1TkMuonClustersRefRemoval,
    TrackQuality = 'undefQuality'
)
