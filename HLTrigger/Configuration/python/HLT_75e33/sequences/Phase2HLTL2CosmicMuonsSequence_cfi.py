import FWCore.ParameterSet.Config as cms

from ..modules.hltL2CosmicMuonsFromL1Muon_cfi import *
from ..modules.hltL2CosmicMuonSeedsFromL1Muon_cfi import *
from ..modules.hltL2CosmicOfflineMuonSeeds_cfi import *

Phase2HLTL2CosmicMuonsSequence = cms.Sequence(hltL2CosmicOfflineMuonSeeds+hltL2CosmicMuonSeedsFromL1Muon+hltL2CosmicMuonsFromL1Muon) 