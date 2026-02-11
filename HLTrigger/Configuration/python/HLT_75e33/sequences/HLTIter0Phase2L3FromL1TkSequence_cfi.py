import FWCore.ParameterSet.Config as cms

from ..modules.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier_cfi import *
from ..modules.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity_cfi import *

HLTIter0Phase2L3FromL1TkSequence = cms.Sequence(
    hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks
    + hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates
    + hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks
    + hltIter0Phase2L3FromL1TkMuonTrackCutClassifier
    + hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
(phase2MuonPixelTracksSelector & phase2CAExtension).toReplaceWith(
    HLTIter0Phase2L3FromL1TkSequence,
    HLTIter0Phase2L3FromL1TkSequence.copyAndExclude(
            [hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks,
            hltIter0Phase2L3FromL1TkMuonTrackCutClassifier]
    )
)

from ..sequences.HLTTrackingSequence_cfi import *
_HLTIter0Phase2L3FromL1TkSequenceGeneralTracks = cms.Sequence(
    HLTTrackingSequence
    + hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks
    + hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity
)

from Configuration.ProcessModifiers.phase2MuonGeneralTracksSelector_cff import phase2MuonGeneralTracksSelector
phase2MuonGeneralTracksSelector.toReplaceWith(
    HLTIter0Phase2L3FromL1TkSequence,
    _HLTIter0Phase2L3FromL1TkSequenceGeneralTracks
)
