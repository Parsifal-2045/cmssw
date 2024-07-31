import FWCore.ParameterSet.Config as cms

from ..modules.hltPFPuppiHT_cfi import *
from ..modules.hltPFPuppiHT1070_cfi import *
from ..modules.l1tPFPuppiHT450off_cfi import *
from ..sequences.hgcalLocalRecoSequence_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTTrackingV61Sequence_cfi import *
from ..sequences.localrecoSequence_cfi import *
from ..sequences.RawToDigiSequence_cfi import *

from ..modules.hltL2MuonSeedsFromL1TkMuon_cfi import PHASE2_TAG
from ..sequences.Phase2HLTMuonsSequence_cfi import *

if PHASE2_TAG:
    HLT_PFPuppiHT1070 = cms.Path(HLTBeginSequence+l1tPFPuppiHT450off+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+Phase2HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+hltPFPuppiHT+hltPFPuppiHT1070+HLTEndSequence)
else:
    HLT_PFPuppiHT1070 = cms.Path(HLTBeginSequence+l1tPFPuppiHT450off+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+hltPFPuppiHT+hltPFPuppiHT1070+HLTEndSequence)
