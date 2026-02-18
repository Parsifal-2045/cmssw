import FWCore.ParameterSet.Config as cms

from RecoMuon.MuonSeedGenerator.MuonSegmentRemover import MuonSegmentRemover as _MuonSegmentRemover
hltMaskedMuonHits = _MuonSegmentRemover(
    standaloneMuons = cms.InputTag("hltL2MuonsFromL1TkMuon"),
    cscSegments = cms.InputTag("hltCscSegments"),
    dtSegments = cms.InputTag("hltDt4DSegments"),
)

from ..modules.hltDisplacedStandaloneMuonSeeds_cfi import *
from ..modules.hltDisplacedStandaloneMuons_cfi import *

HLTDisplacedMuonsSequence = cms.Sequence(
    hltMaskedMuonHits
    + hltDisplacedStandaloneMuonSeeds
    + hltDisplacedStandaloneMuons
)
