import FWCore.ParameterSet.Config as cms

from ..sequences.HLTPhase2MuonPixelTracksFromL1TkSequence_cfi import *
from ..sequences.HLTIter0Phase2L3FromL1TkSequence_cfi import *
from ..sequences.HLTIter2Phase2L3FromL1TkSequence_cfi import *

HLTPhase2L3MuonsIOSequence = cms.Sequence(
    HLTPhase2MuonPixelTracksFromL1TkSequence
    + HLTIter0Phase2L3FromL1TkSequence
    + HLTIter2Phase2L3FromL1TkSequence
)

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
phase2CAExtension.toReplaceWith(
    HLTPhase2L3MuonsIOSequence,
    HLTPhase2L3MuonsIOSequence.copyAndExclude([HLTIter2Phase2L3FromL1TkSequence])
)
