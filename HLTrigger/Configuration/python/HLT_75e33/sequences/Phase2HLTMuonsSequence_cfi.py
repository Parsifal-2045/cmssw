import FWCore.ParameterSet.Config as cms

from ..sequences.HLTIter0Phase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTIter2Phase2L3FromL1TkSequence_cfi import *
from ..sequences.Phase2HLTL2MuonsFromL1TkSequence_cfi import * # includes PHASE2_TAG from HLTrigger/Configuration/python/HLT_75e33/modules/hltL2MuonSeedsFromL1TkMuon_cfi.py
from ..sequences.HLTPhase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTPhase2L3MuonsSequence_cfi import *
from ..sequences.HLTPhase2L3OISequence_cfi import *
from ..modules.Phase2HLTMuonSelectorForL3_cfi import * # includes L3IOFIRST

if L3IOFIRST :
    Phase2HLTMuonsSequence = cms.Sequence(Phase2HLTL2MuonsFromL1TkSequence
                                          +HLTPhase2L3FromL1TkSequence
                                          +HLTIter0Phase2L3FromL1TkSequence
                                          +HLTIter2Phase2L3FromL1TkSequence
                                          +phase2L3FilteredObjects
                                          +HLTPhase2L3OISequence
                                          +HLTPhase2L3MuonsSequence)

else :
    Phase2HLTMuonsSequence = cms.Sequence(Phase2HLTL2MuonsFromL1TkSequence
                                          +HLTPhase2L3OISequence
                                          +phase2L3FilteredObjects
                                          +HLTPhase2L3FromL1TkSequence
                                          +HLTIter0Phase2L3FromL1TkSequence
                                          +HLTIter2Phase2L3FromL1TkSequence
                                          +HLTPhase2L3MuonsSequence)


# Original sequence
#Phase2HLTMuonsSequence = cms.Sequence(Phase2HLTL2MuonsFromL1TkSequence
#                                      +HLTPhase2L3OISequence
#                                      +HLTPhase2L3FromL1TkSequence
#                                      +HLTIter0Phase2L3FromL1TkSequence
#                                      +HLTIter2Phase2L3FromL1TkSequence
#                                      +HLTPhase2L3MuonsSequence)