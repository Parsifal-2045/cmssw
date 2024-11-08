import FWCore.ParameterSet.Config as cms

from ..sequences.HLTIter0Phase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTIter2Phase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTL2MuonsFromL1TkSequence_cfi import *
from ..sequences.HLTPhase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTPhase2L3MuonsSequence_cfi import *
from ..sequences.HLTPhase2L3OISequence_cfi import *

HLTMuonsSequence = cms.Sequence(HLTL2MuonsFromL1TkSequence
                               +HLTPhase2L3OISequence
                               +HLTPhase2L3FromL1TkSequence
                               +HLTIter0Phase2L3FromL1TkSequence
                               +HLTIter2Phase2L3FromL1TkSequence
                               +HLTPhase2L3MuonsSequence)

from ..modules.Phase2HLTMuonSelectorForL3_cfi import *

# The IO first HLT Muons sequence
Phase2HLTMuonsSequenceIOFirst = cms.Sequence(Phase2HLTL2MuonsFromL1TkSequence
                                     +HLTPhase2L3FromL1TkSequence
                                     +HLTIter0Phase2L3FromL1TkSequence
                                     +HLTIter2Phase2L3FromL1TkSequence
                                     +phase2L3FilteredObjects
                                     +HLTPhase2L3OISequence
                                     +HLTPhase2L3MuonsSequence)
# The OI first HLT Muons sequence
Phase2HLTMuonsSequenceOIFirst = cms.Sequence(Phase2HLTL2MuonsFromL1TkSequence
                                     +HLTPhase2L3OISequence
                                     +phase2L3FilteredObjects
                                     +HLTPhase2L3FromL1TkSequence
                                     +HLTIter0Phase2L3FromL1TkSequence
                                     +HLTIter2Phase2L3FromL1TkSequence
                                     +HLTPhase2L3MuonsSequence)

from Configuration.ProcessModifiers.phase2Muon_cff import phase2Muon, L3IOFIRST
if L3IOFIRST :
    phase2Muon.toReplaceWith(HLTMuonsSequence, Phase2HLTMuonsSequenceIOFirst)
else:
    phase2Muon.toReplaceWith(HLTMuonsSequence, Phase2HLTMuonsSequenceOIFirst)
