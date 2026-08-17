import FWCore.ParameterSet.Config as cms

from Validation.RecoTrack.HLTmultiTrackValidator_cfi import *
from SimGeneral.TrackingAnalysis.trackingParticleNumberOfLayersProducer_cff import *
from Validation.RecoTrack.cutsRecoTracks_cfi import cutsRecoTracks as _cutsRecoTracks

hltTrackValidator = hltMultiTrackValidator.clone(
    label = [
        "hltPixelTracks",
        "hltIter0PFlowTrackSelectionHighPurity",
        "hltIter1PFlowTrackSelectionHighPurity",
        "hltIter1Merged",
        "hltIter2PFlowTrackSelectionHighPurity",
        "hltIter2Merged",
        "hltMergedTracks"
    ]
)

# Pixel-less track selector
hltPixelLessTracks = _cutsRecoTracks.clone(
    throwOnMissing = cms.bool(False), # HLT collection might be missing
    src = "hltMergedTracks",
    beamSpot = "hltOnlineBeamSpot",
    minLayer = 3,
    maxPixelHit = 0
)

# Tracks with at least one pixel hit
hltWithPixelTracks = _cutsRecoTracks.clone(
    throwOnMissing = cms.bool(False), # HLT collection might be missing
    src = "hltMergedTracks",
    beamSpot = "hltOnlineBeamSpot",
    minLayer = 3,
    minPixelHit = 1
)

hltMultiTrackValidationTask = cms.Task(
    hltTPClusterProducer
    , trackingParticleNumberOfLayersProducer
    , hltTrackAssociatorByHits
)
hltMultiTrackValidation = cms.Sequence(
    hltTrackValidator,
    hltMultiTrackValidationTask
)

def _modifyForRun3(trackvalidator):
    trackvalidator.label = ["hltPixelTracks", "hltIter0PFlowCtfWithMaterialTracks", "hltIter0PFlowTrackSelectionHighPurity", "hltDoubletRecoveryPFlowCtfWithMaterialTracks", "hltDoubletRecoveryPFlowTrackSelectionHighPurity", "hltMergedTracks"]

from Configuration.Eras.Modifier_run3_common_cff import run3_common
run3_common.toModify(hltTrackValidator, _modifyForRun3)

from Configuration.ProcessModifiers.phase2CAStubs_cff import phase2CAStubs
from Configuration.ProcessModifiers.pixelTrackMask_cff import pixelTrackMask

# Legacy converter template for per-iteration monitoring: the CA iterations write only SoA
# collections, so the validator needs its own conversion. This is the default flavour, without
# phase2CAStubs, and mirrors the baseline hltPhase2PixelTracks converter.
_hltPhase2PixelTracksIterConverter = cms.EDProducer("PixelTrackProducerFromSoAAlpaka",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    trackSrc = cms.InputTag("hltPhase2PixelTrackTorchHighPuritySelector"),
    pixelRecHitLegacySrc = cms.InputTag("hltSiPixelRecHits"),
    outerTrackerRecHitSrc = cms.InputTag("hltSiPhase2RecHits"),
    outerTrackerRecHitSoAConverterSrc = cms.InputTag("hltPhase2OtRecHitsSoA"),
    minNumberOfHits = cms.int32(0),
    minQuality = cms.string('tight'),
    useOTExtension = cms.bool(True),
    requireQuadsFromConsecutiveLayers = cms.bool(False)
)
# Stub-chain flavour: hits from the seeding OT rechit SoA + stub SoA, stubs expanded
# back to their two legacy rechits. Applied to the template before the clones below.
_hltPhase2PixelTracksIterConverterWithStubs = cms.EDProducer("PixelTrackProducerFromSoAAlpaka",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    trackSrc = cms.InputTag("hltPhase2PixelTrackTorchHighPuritySelector"),
    pixelRecHitLegacySrc = cms.InputTag("hltSiPixelRecHits"),
    outerTrackerRecHitSrc = cms.InputTag("hltSiPhase2RecHits"),
    otRecHitsSoASrc = cms.InputTag("hltPixelSeedingOTRecHitsSoA"),
    stubsSoASrc = cms.InputTag("hltOTStubProducer"),
    minNumberOfHits = cms.int32(0),
    minQuality = cms.string('tight'),
    useOTExtension = cms.bool(True),
    expandStubs = cms.bool(True),
    requireQuadsFromConsecutiveLayers = cms.bool(False)
)
phase2CAStubs.toReplaceWith(_hltPhase2PixelTracksIterConverter, _hltPhase2PixelTracksIterConverterWithStubs)

# Four per-iteration states: prompt and displaced, each seen before and after its HP selector.
# Each clones the template (inherits the phase2CAStubs-selected converter flavour).
hltPhase2PixelTracksPromptOnly = _hltPhase2PixelTracksIterConverter.clone(
    trackSrc = "hltPhase2PixelTrackTorchHighPuritySelector",
)

hltPhase2PixelTracksPromptOnlyNoHP = _hltPhase2PixelTracksIterConverter.clone(
    trackSrc = "hltPhase2PixelTracksSoA",
)

hltPhase2PixelTracksDisplacedOnly = _hltPhase2PixelTracksIterConverter.clone(
    trackSrc = "hltPhase2PixelTrackHighPuritySelectorDisplaced",
)

hltPhase2PixelTracksDisplacedOnlyNoHP = _hltPhase2PixelTracksIterConverter.clone(
    trackSrc = "hltPhase2PixelTracksSoADisplaced",
)

# Fifth state: the merged collection before the chain's final high-purity selection. The
# collection after that selection is hltPhase2PixelTracks itself (pixelTrackMask).
hltPhase2PixelTracksMergedNoHP = _hltPhase2PixelTracksIterConverter.clone(
    trackSrc = "hltPhase2PixelTracksSoAMerger",
)

def _modifyForPhase2(trackvalidator):
    trackvalidator.label = ["hltGeneralTracks", "hltPhase2PixelTracks", "hltInitialStepTrackSelectionHighPurity", "hltPixelLessTracks", "hltWithPixelTracks"]

# The per-iteration labels the two-iteration chain monitors on top of the merged hltPhase2PixelTracks.
# Shared by every two-iteration label list below (plain Phase-2 and ngt scouting alike).
_stubsCATwoIterPixelLabels = [
    "hltPhase2PixelTracksPromptOnly",          # prompt iteration alone, after its HighPurity selector
    "hltPhase2PixelTracksPromptOnlyNoHP",      # prompt iteration alone, before its HighPurity selector
    "hltPhase2PixelTracksDisplacedOnly",       # displaced iteration alone, after its HighPurity selector
    "hltPhase2PixelTracksDisplacedOnlyNoHP",   # displaced iteration alone, before its HighPurity selector
    "hltPhase2PixelTracksMergedNoHP",          # both iterations merged, before the final HighPurity selector
]                                              # (after it = hltPhase2PixelTracks, already in every list)

def _modifyForStubsCATwoIter(trackvalidator):
    trackvalidator.label = ["hltGeneralTracks", "hltPhase2PixelTracks"] + _stubsCATwoIterPixelLabels + [
        "hltInitialStepTrackSelectionHighPurity",
        "hltPixelLessTracks",
        "hltWithPixelTracks",
    ]

from Configuration.Eras.Modifier_phase2_tracker_cff import phase2_tracker
phase2_tracker.toModify(hltTrackValidator, _modifyForPhase2)
phase2_tracker.toModify(hltPixelLessTracks, src = "hltGeneralTracks")
phase2_tracker.toModify(hltWithPixelTracks, src = "hltGeneralTracks")

## for Phase 2 only (no pixelless tracks in Run3) run the track selectors
phase2_tracker.toReplaceWith(
    hltMultiTrackValidation,
    cms.Sequence(
        hltPixelLessTracks +
        hltWithPixelTracks +
        hltTrackValidator,
        hltMultiTrackValidationTask
    )
)

from Configuration.ProcessModifiers.ngtScouting_cff import ngtScoutingBase
from Configuration.ProcessModifiers.hltPhase2LegacyTracking_cff import hltPhase2LegacyTracking
from Configuration.ProcessModifiers.trackingLST_cff import trackingLST

def _modifyForHLTPhase2LegacyTracking(trackvalidator):
    trackvalidator.label = ["hltGeneralTracks", "hltPhase2PixelTracks", "hltInitialStepTrackSelectionHighPurity", "hltHighPtTripletStepTrackSelectionHighPurity", "hltPixelLessTracks", "hltWithPixelTracks"]
hltPhase2LegacyTracking.toModify(hltTrackValidator, _modifyForHLTPhase2LegacyTracking)

def _modifyForNGTScouting(trackvalidator):
    trackvalidator.label = ["hltGeneralTracks", "hltPhase2PixelTracks"]
(ngtScoutingBase & ~trackingLST).toModify(hltTrackValidator, _modifyForNGTScouting)

def _modifyForNGTScoutingLST(trackvalidator):
    trackvalidator.label = ["hltGeneralTracks", "hltPhase2PixelTracks", "hltInitialStepTracksT4T5TCLST", "hltPixelLessTracks", "hltWithPixelTracks"]
(ngtScoutingBase & trackingLST).toModify(hltTrackValidator, _modifyForNGTScoutingLST)
# The per-iteration converters exist only for the SoA collections of the two-iteration chain and
# are scheduled and monitored only there (phase2CAStubs & pixelTrackMask). ngtScouting is a
# ModifierChain over phase2CAStubs and pixelTrackMask, so the plain two-iteration list must not
# fire under it: the scouting arms get their own combined lists below.
def _modifyForNGTScoutingStubsCATwoIter(trackvalidator):
    trackvalidator.label = ["hltGeneralTracks", "hltPhase2PixelTracks"] + _stubsCATwoIterPixelLabels

def _modifyForNGTScoutingLSTStubsCATwoIter(trackvalidator):
    trackvalidator.label = ["hltGeneralTracks", "hltPhase2PixelTracks"] + _stubsCATwoIterPixelLabels + [
        "hltInitialStepTracksT4T5TCLST",
        "hltPixelLessTracks",
        "hltWithPixelTracks",
    ]

(phase2CAStubs & pixelTrackMask & ~ngtScoutingBase).toModify(hltTrackValidator, _modifyForStubsCATwoIter)
(phase2CAStubs & pixelTrackMask & ngtScoutingBase & ~trackingLST).toModify(hltTrackValidator, _modifyForNGTScoutingStubsCATwoIter)
(phase2CAStubs & pixelTrackMask & ngtScoutingBase & trackingLST).toModify(hltTrackValidator, _modifyForNGTScoutingLSTStubsCATwoIter)

(phase2CAStubs & pixelTrackMask).toReplaceWith(
    hltMultiTrackValidation,
    cms.Sequence(
        hltPixelLessTracks +
        hltWithPixelTracks +
        hltPhase2PixelTracksPromptOnly +
        hltPhase2PixelTracksPromptOnlyNoHP +
        hltPhase2PixelTracksDisplacedOnly +
        hltPhase2PixelTracksDisplacedOnlyNoHP +
        hltPhase2PixelTracksMergedNoHP +
        hltTrackValidator,
        hltMultiTrackValidationTask
    )
)
