import FWCore.ParameterSet.Config as cms

from ..modules.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracks_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonPixelVertices_cfi import *
from ..modules.hltPhase2L3FromL1TkMuonTrimmedPixelVertices_cfi import *

from ..modules.hltPhase2PixelFitterByHelixProjections_cfi import hltPhase2PixelFitterByHelixProjections
from ..modules.hltPhase2PixelTrackFilterByKinematics_cfi import hltPhase2PixelTrackFilterByKinematics

HLTPhase2MuonPixelTracksFromL1TkSequence = cms.Sequence(
    hltPhase2L3FromL1TkMuonPixelLayerQuadruplets
    + hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions
    + hltPhase2L3FromL1TkMuonPixelTracksHitDoublets
    + hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets
    + hltPhase2PixelTrackFilterByKinematics
    + hltPhase2PixelFitterByHelixProjections
    + hltPhase2L3FromL1TkMuonPixelTracks
    + hltPhase2L3FromL1TkMuonPixelVertices
    + hltPhase2L3FromL1TkMuonTrimmedPixelVertices
)

from ..sequences.HLTPhase2PixelTracksAndVerticesSequence_cfi import *
from ..modules.hltPhase2MuonPixelTracks_cfi import *
from ..modules.hltPhase2MuonPixelTracksHighPurity_cfi import *
_HLTPhase2MuonPixelTracksSelectorSequence = cms.Sequence(
    HLTBeamSpotSequence
    + hltPhase2OtRecHitsSoA
    + hltPhase2PixelRecHitsExtendedSoA
    + hltPhase2PixelTracksSoA
    + hltPhase2PixelTracksCAExtension
    + hltPhase2MuonPixelTracks
    + hltPhase2MuonPixelTracksHighPurity
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.ngtScouting_cff import ngtScouting
(phase2MuonPixelTracksSelector | ngtScouting).toReplaceWith(
    HLTPhase2MuonPixelTracksFromL1TkSequence,
    _HLTPhase2MuonPixelTracksSelectorSequence
)
