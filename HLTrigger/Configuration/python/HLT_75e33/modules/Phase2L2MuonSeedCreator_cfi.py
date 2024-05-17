import FWCore.ParameterSet.Config as cms

Phase2hltL2MuonSeedsFromL1TkMuon = cms.EDProducer('Phase2L2MuonSeedCreator',
  InputObjects = cms.InputTag('l1tTkMuonsGmt'),
  CSCRecSegmentLabel = cms.InputTag('hltCscSegments'),
  DTRecSegmentLabel = cms.InputTag('hltDt4DSegments'),
  MinPL1Tk = cms.double(3.5),
  MaxPL1Tk = cms.double(200),
  StubMatchDR = cms.double(0.25),
  maximumEtaBarrel = cms.double(0.7),
  maximumEtaOverlap = cms.double(1.3),
  mightGet = cms.optional.untracked.vstring
)
