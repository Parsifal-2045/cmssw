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

from Configuration.ProcessModifiers.trackingLST_cff import trackingLST

from ..modules.hltInputLST_cfi import *
from ..modules.hltInitialStepMkFitSeeds_cfi import *
from ..modules.hltInitialStepSeedTracksLST_cfi import *
from ..modules.hltInitialStepSeeds_cfi import *
from ..modules.hltInitialStepTrackCandidates_cfi import *
from ..modules.hltInitialStepTrackCandidatesMkFit_cfi import *
from ..modules.hltInitialStepTrackCutClassifier_cfi import *
from ..modules.hltInitialStepTrackSelectionHighPurity_cfi import *
from ..modules.hltInitialStepTracks_cfi import *
from ..modules.hltInitialStepTrajectorySeedsLST_cfi import *
from ..modules.hltLST_cfi import *
from ..modules.hltSiPhase2RecHits_cfi import *
from ..sequences.HLTMkFitInputSequence_cfi import *

from RecoMuon.L3TrackFinder.MuonSeedsSelectorFromL1TkMuon import MuonSeedsSelectorFromL1TkMuon as _MuonSeedsSelectorFromL1TkMuon 
hltIOMuonTracksSeeds = _MuonSeedsSelectorFromL1TkMuon()

#hltInitialStepSeedTracksLST = cms.EDProducer(
#    "TrackFromSeedProducer",
#    src = cms.InputTag("hltInitialStepTrajectorySeedsLST"),
#    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
#    TTRHBuilder = cms.string("hltESPTTRHBuilderWithoutRefit")
#)

hltIOMuonsMkFitSeeds = hltInitialStepMkFitSeeds.clone(
    seeds = "hltIOMuonTracksSeeds",
)

hltIOMuonsTrackCandidatesMkFit = hltInitialStepTrackCandidatesMkFit.clone(
    seeds = "hltIOMuonsMkFitSeeds",
)

hltIOMuonsTrackCandidates = hltInitialStepTrackCandidates.clone()
trackingLST.toModify(hltIOMuonsTrackCandidates,
                     mkFitSeeds = "hltIOMuonsMkFitSeeds",
                     seeds = "hltIOMuonTracksSeeds",
                     tracks = "hltIOMuonsTrackCandidatesMkFit",
)

_HLTIter0Phase2L3FromL1TkSequenceGeneralTracksLST = cms.Sequence(
    HLTItLocalRecoSequence
    + HLTOtLocalRecoSequence
    + hltTrackerClusterCheck
    + HLTPhase2PixelTracksAndVerticesSequence
    + hltInitialStepSeeds
    + hltInitialStepSeedTracksLST
    + hltSiPhase2RecHits
    + hltInputLST
    + hltLST
    + hltInitialStepTrajectorySeedsLST
    + hltIOMuonTracksSeeds
    + HLTMkFitInputSequence
    + hltIOMuonsMkFitSeeds
    + hltIOMuonsTrackCandidatesMkFit
    + hltIOMuonsTrackCandidates
    + hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks
    + hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity 
)

(phase2MuonGeneralTracksSelector & trackingLST).toReplaceWith(
    HLTIter0Phase2L3FromL1TkSequence,
    _HLTIter0Phase2L3FromL1TkSequenceGeneralTracksLST
)
