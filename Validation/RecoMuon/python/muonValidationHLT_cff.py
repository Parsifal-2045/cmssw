import FWCore.ParameterSet.Config as cms

from Validation.RecoMuon.selectors_cff import *
from Validation.RecoMuon.track_selectors_cff import *
from Validation.RecoMuon.associators_cff import *
from Validation.RecoMuon.histoParameters_cff import *

import Validation.RecoMuon.MuonTrackValidator_cfi
MTVhlt = Validation.RecoMuon.MuonTrackValidator_cfi.muonTrackValidator.clone(
# DEFAULTS ###################################
#    label_tp = "mix:MergedTrackTruth",
#    label_tp_refvector = False,
#    muonTPSelector = dict(muonTPSet),
##############################################
label_tp = ("TPmu"),
label_tp_refvector = True,
dirName = 'HLT/Muon/MuonTrack/',
#beamSpot = 'hltOfflineBeamSpot',
ignoremissingtrackcollection=True
)
MTVhlt.muonTPSelector.src = ("TPmu")
################################################

from SimMuon.MCTruth.SeedToTrackProducer_cfi import SeedToTrackProducer as _SeedToTrackProducer
muonL2SeedTracks = _SeedToTrackProducer.clone(L2seedsCollection = cms.InputTag("hltL2MuonSeedsFromL1TkMuon"))

# L2 standalone muons seeds
l2MuSeedV = MTVhlt.clone(
    associatormap = 'tpToL2SeedAssociation',
    label = ('muonL2SeedTracks',),
    muonHistoParameters = staSeedMuonHistoParameters
)

# L2 standalone muons
l2MuV = MTVhlt.clone(
    associatormap = 'tpToL2MuonAssociation',
    label = ('hltL2MuonsFromL1TkMuon:UpdatedAtVtx',),
    muonHistoParameters = staMuonHistoParameters
)

# L3 IO inner tracks
l3IOTkV = MTVhlt.clone(
    associatormap = 'tpToL3IOTkAssociation',
    label = ('hltIter2Phase2L3FromL1TkMuonMerged',),
    muonHistoParameters = trkMuonHistoParameters
)

# L3 OI inner tracks
l3OITkV = MTVhlt.clone(
    associatormap = 'tpToL3OITkAssociation',
    label = ('hltPhase2L3OIMuonTrackSelectionHighPurity',),
    muonHistoParameters = trkMuonHistoParameters
)

from HLTrigger.Configuration.HLT_75e33.modules.hltL2MuonSeedsFromL1TkMuon_cfi import PHASE2_TAG

from HLTrigger.Configuration.HLT_75e33.modules.Phase2HLTMuonSelectorForL3_cfi import L3IOFIRST
# Filtered objects (L2 Muons to reuse and L3 IO tracks filtered if IO first, L3 OI tracks filtered if OI first)
if PHASE2_TAG:    
    if L3IOFIRST:
        # L2 muons to reuse (in workflow where IO is done first and OI as a second pass)
        L2MuToReuseV = MTVhlt.clone(
            associatormap = 'tpToL2MuonToReuseAssociation',
            label = ('phase2L3FilteredObjects:L2MuToReuse',),
            muonHistoParameters = staMuonHistoParameters
        )

        # L3 IO inner tracks filtered
        l3IOTkFilteredV = MTVhlt.clone(
            associatormap = 'tpToL3IOTkFilteredAssociation',
            label = ('phase2L3FilteredObjects:L3IOTracksFiltered',),
            muonHistoParameters = trkMuonHistoParameters
        )

    else:
        # L3 OI inner tracks filtered
        l3OITkFilteredV = MTVhlt.clone(
            associatormap = 'tpToL3OITkFilteredAssociation',
            label = ('phase2L3FilteredObjects:L3OITracksFiltered',),
            muonHistoParameters = trkMuonHistoParameters
        )

# L3 inner tracks merged
l3TkMergedV = MTVhlt.clone(
    associatormap = 'tpToL3TkMergedAssociation',
    label = ('hltPhase2L3MuonMerged',),
    muonHistoParameters = trkMuonHistoParameters
)

# L3 global muons
l3GlbMuonV = MTVhlt.clone(
    associatormap = 'tpToL3GlbMuonMergedAssociation',
    label = ('hltPhase2L3GlbMuon',),
    muonHistoParameters = glbMuonHistoParameters
)

# L3 OI global muons
l3OIGlbMuonV = MTVhlt.clone(
    associatormap = 'tpToL3OIMuonAssociation',
    label = ('hltL3MuonsPhase2L3OI',),
    muonHistoParameters = glbMuonHistoParameters
)

# L3 Muons no ID not tracks nor standalone nor global cannot extract as they are
#l3NoIDV = MTVhlt.clone(
#    associatormap = 'tpToL3MuonNoIDAssociation',
#    label = ('hltPhase2L3MuonsNoID',),
#    muonHistoParameters = glbMuonHistoParameters
#)

# L3 Muons ID
l3MuIDTrackV = MTVhlt.clone(
    associatormap = 'tpToL3MuonIDAssociation',
    label = ('hltL3MuonIdTracks',),
    muonHistoParameters = glbMuonHistoParameters
)

#l2UpdMuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL2UpdMuonAssociation',
#    label = ('hltL2Muons:UpdatedAtVtx',),
#    muonHistoParameters = staUpdMuonHistoParameters
#)
#l3OITkMuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL3OITkMuonAssociation',
#    label = ('hltIterL3OIMuonTrackSelectionHighPurity:',),
#    muonHistoParameters = trkMuonHistoParameters
#)
#l3TkMuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL3TkMuonAssociation',
#    label = ('hltIterL3MuonMerged:',),
#    muonHistoParameters = trkMuonHistoParameters
#)
#l3IOFromL1TkMuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL3FromL1TkMuonAssociation',
#    label = ('hltIterL3MuonAndMuonFromL1Merged:',),
#    muonHistoParameters = trkMuonHistoParameters
#)
#l0l3FromL1TkMuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL0L3FromL1TkMuonAssociation',
#    label = ('hltIter0IterL3FromL1MuonTrackSelectionHighPurity:',),
#    muonHistoParameters = trkMuonHistoParameters
#)
#l3GlbMuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL3GlbMuonAssociation',
#    label = ('hltIterL3GlbMuon:',),
#    muonHistoParameters = glbMuonHistoParameters
#)
#l3NoIDMuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL3NoIDMuonAssociation',
#    label = ('hltIterL3MuonsNoIDTracks:',),
#    muonHistoParameters = glbMuonHistoParameters
#)
#l3MuonMuTrackV = MTVhlt.clone(
#    associatormap = 'tpToL3MuonAssociation',
#    label = ('hltIterL3MuonsTracks:',),
#    muonHistoParameters = glbMuonHistoParameters
#)
#
# The full Muon HLT validation sequence
#

if PHASE2_TAG: 
    if L3IOFIRST: 
        muonValidationHLT_seq = cms.Sequence(
            muonL2SeedTracks + tpToL2SeedAssociation + l2MuSeedV
            +tpToL2MuonAssociation + l2MuV
            +tpToL2MuonToReuseAssociation + L2MuToReuseV
            +tpToL3IOTkAssociation + l3IOTkV
            +tpToL3OITkAssociation + l3OITkV
            +tpToL3IOTkFilteredAssociation + l3IOTkFilteredV
            +tpToL3TkMergedAssociation + l3TkMergedV
            +tpToL3GlbMuonMergedAssociation + l3GlbMuonV
            +tpToL3OIMuonAssociation + l3OIGlbMuonV
            +hltL3MuonIdTracks_seq + tpToL3MuonIDAssociation + l3MuIDTrackV
            )
    else: 
        muonValidationHLT_seq = cms.Sequence(
            muonL2SeedTracks + tpToL2SeedAssociation + l2MuSeedV
            +tpToL2MuonAssociation + l2MuV
            +tpToL3IOTkAssociation + l3IOTkV
            +tpToL3OITkAssociation + l3OITkV
            +tpToL3OITkFilteredAssociation + l3OITkFilteredV
            +tpToL3TkMergedAssociation + l3TkMergedV
            +tpToL3GlbMuonMergedAssociation + l3GlbMuonV
            +tpToL3OIMuonAssociation + l3OIGlbMuonV
            +hltL3MuonIdTracks_seq + tpToL3MuonIDAssociation + l3MuIDTrackV
            )
else: 
    muonValidationHLT_seq = cms.Sequence(
        muonL2SeedTracks + tpToL2SeedAssociation + l2MuSeedV
        +tpToL2MuonAssociation + l2MuV
        +tpToL3IOTkAssociation + l3IOTkV
        +tpToL3OITkAssociation + l3OITkV
        +tpToL3TkMergedAssociation + l3TkMergedV
        +tpToL3GlbMuonMergedAssociation + l3GlbMuonV
        +tpToL3OIMuonAssociation + l3OIGlbMuonV
        +hltL3MuonIdTracks_seq + tpToL3MuonIDAssociation + l3MuIDTrackV
        )

recoMuonValidationHLT_seq = cms.Sequence(
    cms.SequencePlaceholder("TPmu") +
    muonValidationHLT_seq
    )
