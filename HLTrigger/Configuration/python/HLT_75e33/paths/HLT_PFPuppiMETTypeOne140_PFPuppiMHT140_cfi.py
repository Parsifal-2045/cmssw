import FWCore.ParameterSet.Config as cms

from ..modules.hltPFPuppiMETTypeOne_cfi import *
from ..modules.hltPFPuppiMETTypeOne140_cfi import *
from ..modules.hltPFPuppiMETTypeOneCorrector_cfi import *
from ..modules.hltPFPuppiMHT_cfi import *
from ..modules.hltPFPuppiMHT140_cfi import *
from ..modules.hltL1SeedsForPuppiMETFilter_cfi import *
from ..sequences.hgcalLocalRecoSequence_cfi import *
from ..sequences.HLTAK4PFPuppiJetsReconstruction_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.HLTParticleFlowSequence_cfi import *
from ..sequences.HLTPFPuppiMETReconstruction_cfi import *
from ..sequences.HLTTrackingV61Sequence_cfi import *
from ..sequences.localrecoSequence_cfi import *
from ..sequences.RawToDigiSequence_cfi import *

from ..modules.hltL2MuonSeedsFromL1TkMuon_cfi import PHASE2_TAG
from ..sequences.Phase2HLTMuonsSequence_cfi import *

if PHASE2_TAG:
    HLT_PFPuppiMETTypeOne140_PFPuppiMHT140 = cms.Path(HLTBeginSequence+hltL1SeedsForPuppiMETFilter+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+Phase2HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+HLTPFPuppiMETReconstruction+hltPFPuppiMETTypeOneCorrector+hltPFPuppiMETTypeOne+hltPFPuppiMETTypeOne140+hltPFPuppiMHT+hltPFPuppiMHT140+HLTEndSequence)
else:
    HLT_PFPuppiMETTypeOne140_PFPuppiMHT140 = cms.Path(HLTBeginSequence+hltL1SeedsForPuppiMETFilter+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+HLTPFPuppiMETReconstruction+hltPFPuppiMETTypeOneCorrector+hltPFPuppiMETTypeOne+hltPFPuppiMETTypeOne140+hltPFPuppiMHT+hltPFPuppiMHT140+HLTEndSequence)
