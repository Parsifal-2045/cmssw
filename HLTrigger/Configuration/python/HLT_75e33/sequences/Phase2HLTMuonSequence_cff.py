import FWCore.ParameterSet.Config as cms

from ..sequences.HLTIter0Phase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTIter2Phase2L3FromL1TkSequence_cfi import *
from ..sequences.Phase2HLTL2MuonsFromL1TkSequence_cfi import *
from ..sequences.HLTPhase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTPhase2L3MuonsSequence_cfi import *
from ..sequences.HLTPhase2L3OISequence_cfi import *

Phase2HLTMuonSequence = cms.Sequence(Phase2HLTL2MuonsFromL1TkSequence+HLTPhase2L3OISequence+HLTPhase2L3FromL1TkSequence+HLTIter0Phase2L3FromL1TkSequence+HLTIter2Phase2L3FromL1TkSequence+HLTPhase2L3MuonsSequence)