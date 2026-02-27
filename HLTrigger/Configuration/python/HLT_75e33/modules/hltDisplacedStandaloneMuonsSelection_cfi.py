import FWCore.ParameterSet.Config as cms

hltDisplacedStandaloneMuonsSelection = cms.EDFilter("TrackSelector",
    src = cms.InputTag("hltDisplacedStandaloneMuons"),
    cut = cms.string("numberOfValidHits >= 4 && normalizedChi2 < 25 && pt > 2.5")
)
