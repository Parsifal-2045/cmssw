import FWCore.ParameterSet.Config as cms

from ..modules.hltSingleAK4PFPuppiJet520_cfi import *
from ..modules.hltL1SeedsForPuppiJetFilter_cfi import *
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
    HLT_AK4PFPuppiJet520 = cms.Path(HLTBeginSequence+hltL1SeedsForPuppiJetFilter+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+Phase2HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+hltSingleAK4PFPuppiJet520+HLTEndSequence)
else :
    HLT_AK4PFPuppiJet520 = cms.Path(HLTBeginSequence+hltL1SeedsForPuppiJetFilter+RawToDigiSequence+hgcalLocalRecoSequence+localrecoSequence+HLTTrackingV61Sequence+HLTMuonsSequence+HLTParticleFlowSequence+HLTAK4PFPuppiJetsReconstruction+hltSingleAK4PFPuppiJet520+HLTEndSequence)
