import FWCore.ParameterSet.Config as cms

HLTPSetDisplacedMuonCkfTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(10.0),  # permissive for displaced
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(10),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(10),
    maxLostHitsFraction = cms.double(0.25),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(2.5),
    minimumNumberOfHits = cms.int32(4),
    nSigmaMinPt = cms.double(4.0),
    pixelSeedExtension = cms.bool(True),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)
