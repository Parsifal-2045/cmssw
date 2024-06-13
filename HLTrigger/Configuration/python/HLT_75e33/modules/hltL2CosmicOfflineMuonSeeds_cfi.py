import FWCore.ParameterSet.Config as cms

hltL2CosmicOfflineMuonSeeds = cms.EDProducer( "CosmicMuonSeedGenerator",
    MaxSeeds = cms.int32( 1000 ),
    CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
    EnableDTMeasurement = cms.bool( True ),
    MaxCSCChi2 = cms.double( 300.0 ),
    MaxDTChi2 = cms.double( 300.0 ),
    DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
    EnableCSCMeasurement = cms.bool( True ),
    ForcePointDown = cms.bool( False )
)