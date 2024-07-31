import FWCore.ParameterSet.Config as cms

from ..sequences.hgcalLocalRecoSequence_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTBtagDeepCSVSequencePFPuppi_cfi import *
from ..sequences.HLTBtagDeepFlavourSequencePFPuppi_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTTrackingV61Sequence_cfi import *
from ..sequences.localrecoSequence_cfi import *
from ..sequences.RawToDigiSequence_cfi import *

from ..modules.hltL2MuonSeedsFromL1TkMuon_cfi import PHASE2_TAG
from ..sequences.Phase2HLTMuonsSequence_cfi import *

if PHASE2_TAG:
    MC_BTV = cms.Path(HLTBeginSequence+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+Phase2HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+HLTBtagDeepCSVSequencePFPuppi+HLTBtagDeepFlavourSequencePFPuppi)
else:
    MC_BTV = cms.Path(HLTBeginSequence+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+HLTBtagDeepCSVSequencePFPuppi+HLTBtagDeepFlavourSequencePFPuppi)
