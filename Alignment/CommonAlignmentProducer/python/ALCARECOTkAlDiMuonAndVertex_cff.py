import FWCore.ParameterSet.Config as cms

##################################################################
# Exact same configuration as TkAlZMuMu: extract mumu pairs
#################################################################
from Alignment.CommonAlignmentProducer.TkAlMuonSelectors_cfi import *
import Alignment.CommonAlignmentProducer.ALCARECOTkAlZMuMu_cff as confALCARECOTkAlZMuMu
ALCARECOTkAlDiMuonHLT = confALCARECOTkAlZMuMu.ALCARECOTkAlZMuMuHLT.clone()
ALCARECOTkAlDiMuonDCSFilter = confALCARECOTkAlZMuMu.ALCARECOTkAlZMuMuDCSFilter.clone()
ALCARECOTkAlDiMuon = confALCARECOTkAlZMuMu.ALCARECOTkAlZMuMu.clone()
ALCARECOTkAlDiMuon.GlobalSelector.muonSource = 'TkAlRelCombIsoMuonSelector'

##################################################################
# Tracks from the selected vertex
#################################################################
import Alignment.CommonAlignmentProducer.AlignmentTracksFromVertexSelector_cfi as TracksFromVertex
ALCARECOTkAlDiMuonVertexTracks = TracksFromVertex.AlignmentTracksFromVertexSelector.clone()

##################################################################
# for the GEN level information
##################################################################
TkAlDiMuonAndVertexGenMuonSelector = cms.EDFilter("GenParticleSelector",
                                                  src = cms.InputTag("genParticles"),
                                                  cut = cms.string("abs(pdgId) == 13"), # Select only muons
                                                  filter = cms.bool(False),
                                                  throwOnMissing = cms.untracked.bool(False))

##################################################################
# The sequence
#################################################################
seqALCARECOTkAlDiMuonAndVertex = cms.Sequence(ALCARECOTkAlDiMuonHLT+
                                              ALCARECOTkAlDiMuonDCSFilter+
                                              seqALCARECOTkAlRelCombIsoMuons+
                                              ALCARECOTkAlDiMuon+
                                              ALCARECOTkAlDiMuonVertexTracks+
                                              TkAlDiMuonAndVertexGenMuonSelector)

## customizations for the pp_on_AA eras
from Configuration.Eras.Modifier_pp_on_XeXe_2017_cff import pp_on_XeXe_2017
from Configuration.ProcessModifiers.pp_on_AA_cff import pp_on_AA
(pp_on_XeXe_2017 | pp_on_AA).toModify(ALCARECOTkAlDiMuonHLT,
                                      eventSetupPathsKey='TkAlZMuMuHI')
from Configuration.Eras.Modifier_phase2_tracker_cff import phase2_tracker
phase2_tracker.toModify(ALCARECOTkAlDiMuon, etaMin = -4, etaMax = 4)
