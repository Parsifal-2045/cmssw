import FWCore.ParameterSet.Config as cms

Phase2hltL2MuonSeedsFromL1TkMuon = cms.EDProducer('Phase2L2MuonSeedCreator',
  InputObjects = cms.InputTag('l1tTkMuonsGmt'),
  CSCRecSegmentLabel = cms.InputTag('hltCscSegments'),
  DTRecSegmentLabel = cms.InputTag('hltDt4DSegments'),
  MinPL1Tk = cms.double(3.5),
  MaxPL1Tk = cms.double(200),
  stubMatchDPhi = cms.double(0.05),
  stubMatchDTheta = cms.double(0.1),
  extrapolationWindowClose = cms.double(0.1),
  extrapolationWindowFar = cms.double(0.05),
  maximumEtaBarrel = cms.double(0.7),
  maximumEtaOverlap = cms.double(1.3),
  Propagator = cms.string('SteppingHelixPropagatorAny'),
  ServiceParameters = cms.PSet(
      Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny'),
      RPCLayers = cms.bool(True),
      UseMuonNavigation = cms.untracked.bool(True)
  ),
  mightGet = cms.optional.untracked.vstring
)
