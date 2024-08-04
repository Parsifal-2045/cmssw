import FWCore.ParameterSet.Config as cms

from ..modules.hltDoubleMuon7DZ1p0_cfi import *
from ..modules.hltL3fL1DoubleMu155fFiltered37_cfi import *
from ..modules.hltL3fL1DoubleMu155fPreFiltered27_cfi import *
from ..modules.hltPhase2L3MuonCandidates_cfi import *
from ..modules.hltPhase2PixelFitterByHelixProjections_cfi import *
from ..modules.hltPhase2PixelTrackFilterByKinematics_cfi import *
from ..sequences.HLTBeginSequence_cfi import *
from ..sequences.HLTEndSequence_cfi import *
from ..sequences.HLTMuonsSequence_cfi import *
from ..sequences.itLocalRecoSequence_cfi import *
from ..sequences.muonlocalrecoSequence_cfi import *
from ..sequences.otLocalRecoSequence_cfi import *

from ..sequences.Phase2HLTMuonsSequence_cfi import * # includes PHASE2_TAG from HLTrigger/Configuration/python/HLT_75e33/modules/hltL2MuonSeedsFromL1TkMuon_cfi.py

if PHASE2_TAG :
    HLT_Mu37_Mu27_FromL1TkMuon = cms.Path(HLTBeginSequence
    +hltDoubleMuon7DZ1p0
    +muonlocalrecoSequence
    +itLocalRecoSequence
    +otLocalRecoSequence
    +hltPhase2PixelFitterByHelixProjections
    +hltPhase2PixelTrackFilterByKinematics
    +Phase2HLTMuonsSequence
    +hltL3fL1DoubleMu155fPreFiltered27
    +hltL3fL1DoubleMu155fFiltered37
    +HLTEndSequence)

else : 
    HLT_Mu37_Mu27_FromL1TkMuon = cms.Path(HLTBeginSequence
    +hltDoubleMuon7DZ1p0
    +muonlocalrecoSequence
    +itLocalRecoSequence
    +otLocalRecoSequence
    +hltPhase2PixelFitterByHelixProjections
    +hltPhase2PixelTrackFilterByKinematics
    +HLTMuonsSequence
    +hltPhase2L3MuonCandidates
    +hltL3fL1DoubleMu155fPreFiltered27
    +hltL3fL1DoubleMu155fFiltered37
    +HLTEndSequence)
