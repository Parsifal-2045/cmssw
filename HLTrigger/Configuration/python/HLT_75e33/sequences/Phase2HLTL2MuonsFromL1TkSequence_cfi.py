import FWCore.ParameterSet.Config as cms

from ..sequences.muonlocalrecoSequence_cfi import *
from ..modules.Phase2L2MuonSeedCreator_cfi import *


Phase2HLTL2MuonsFromL1TkSequence = cms.Sequence(muonlocalrecoSequence+Phase2L2MuonSeedCreator+hltL2MuonsFromL1TkMuon)