import FWCore.ParameterSet.Config as cms

from ..sequences.muonlocalrecoSequence_cfi import *
from ..modules.hltL2MuonSeedsFromL1TkMuon_cfi import * # includes PHASE2_TAG
from ..modules.hltL2MuonsFromL1TkMuon_cfi import *

Phase2HLTL2MuonsFromL1TkSequence = cms.Sequence(muonlocalrecoSequence+hltL2MuonSeedsFromL1TkMuon+hltL2MuonsFromL1TkMuon)