import FWCore.ParameterSet.Config as cms

from ..modules.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracks_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelVertices_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonTrimmedPixelVertices_cfi import *

HLTPhase2MuonPixelTracksFromL1TkSequence = cms.Sequence(
    hltPhase2L3FromL1TkMuonPixelLayerQuadruplets
    + hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions
    + hltPhase2L3FromL1TkMuonPixelTracksHitDoublets
    + hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets
    + hltPhase2L3FromL1TkMuonPixelTracks
    + hltPhase2L3FromL1TkMuonPixelVertices
    + hltPhase2L3FromL1TkMuonTrimmedPixelVertices
)

from ..sequences.HLTPhase2PixelTracksAndVerticesSequence_cfi import *
from ..modules.hltPhase2MuonPixelTracksDNNSelector_cfi import *
_HLTPhase2MuonPixelTracksFromL1TkSequence = cms.Sequence(
    HLTPhase2PixelTracksAndVerticesSequence.copyAndExclude(
        [hltPhase2PixelTracksAndHighPtStepTrackingRegions]
    )
    + hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions 
    + hltPhase2L3FromL1TkMuonPixelTracks
    + hltPhase2MuonPixelTracksDNNSelector
)

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
phase2CAExtension.toReplaceWith(
    HLTPhase2MuonPixelTracksFromL1TkSequence,
    _HLTPhase2MuonPixelTracksFromL1TkSequence
)

