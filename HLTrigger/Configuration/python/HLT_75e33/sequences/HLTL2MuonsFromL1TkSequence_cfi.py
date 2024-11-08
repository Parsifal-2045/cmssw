import FWCore.ParameterSet.Config as cms

from ..modules.hltL2MuonSeedsFromL1TkMuon_cfi import *
from ..modules.hltL2MuonsFromL1TkMuon_cfi import *
from ..modules.hltL2OfflineMuonSeeds_cfi import *
from ..sequences.HLTMuonlocalrecoSequence_cfi import *

HLTL2MuonsFromL1TkSequence = cms.Sequence(HLTMuonlocalrecoSequence
                                         +hltL2OfflineMuonSeeds
                                         +hltL2MuonSeedsFromL1TkMuon
                                         +hltL2MuonsFromL1TkMuon)


Phase2HLTL2MuonsFromL1TkSequence = cms.Sequence(HLTMuonlocalrecoSequence
                                               +hltL2MuonSeedsFromL1TkMuon
                                               +hltL2MuonsFromL1TkMuon)

from Configuration.ProcessModifiers.phase2Muon_cff import phase2Muon
phase2Muon.toReplaceWith(HLTL2MuonsFromL1TkSequence, Phase2HLTL2MuonsFromL1TkSequence)
