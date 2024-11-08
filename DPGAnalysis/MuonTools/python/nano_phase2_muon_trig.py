import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *

from PhysicsTools.NanoAOD.l1trig_cff import *

from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

from Validation.RecoMuon.muonValidationHLT_cff import *

# L1Tk Muons
l1TkMuTable = cms.EDProducer("SimpleTriggerL1TkMuonFlatTableProducer",
    src = cms.InputTag("l1tTkMuonsGmt"),
    cut = cms.string(""), 
    name= cms.string("L1TkMu"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("phPt()", "float", doc="Physics pt"),
                         eta = Var("phEta()", "float", doc="#eta"),
                         phi = Var("phPhi()", "float", doc="#phi (rad)"),
                         dXY = Var("phD0()", "float", doc="dXY (cm)"),
                         dZ = Var("phZ0()", "float", doc="dZ (cm)")
                     )
)
    
# L2 offline seeds
l2SeedTable = cms.EDProducer("SimpleTrajectorySeedFlatTableProducer",
    src = cms.InputTag("hltL2OfflineMuonSeeds"),
    cut = cms.string(""), 
    name= cms.string("l2_seed_offline"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("startingState().pt()", "float", doc="p_T (GeV)"),
                     nHits = Var("nHits()", "int16", doc="")
                        )
)

# L2 seeds from L1Tk Muons
l2SeedFromL1TkMuonTable = cms.EDProducer("SimpleTrajectorySeedFlatTableProducer",
    src = cms.InputTag("hltL2MuonSeedsFromL1TkMuon"),
    cut = cms.string(""), 
    name= cms.string("l2_seed"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("startingState().pt()", "float", doc="p_T (GeV)"),
                         nHits = Var("nHits()", "int16", doc="")
                         )
)

# L2 standalone muons
l2MuTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltL2MuonsFromL1TkMuon"),
    cut = cms.string(""), 
    name= cms.string("l2_mu"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                         )
)

# L2 standalone muons updated at vertex
l2MuTableVtx = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltL2MuonsFromL1TkMuon","UpdatedAtVtx"),
    cut = cms.string(""), 
    name= cms.string("l2_mu_vtx"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                         )
)

# L3 IO inner tracks
l3TkIOTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMerged"),
    cut = cms.string(""), 
    name= cms.string("l3_tk_IO"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                         )
)

# L3 OI inner tracks
l3TkOITable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
    cut = cms.string(""), 
    name= cms.string("l3_tk_OI"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                        )
)

# L3 tracks merged
l3TkMergedTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltPhase2L3MuonMerged"),
    cut = cms.string(""), 
    name= cms.string("l3_tk_merged"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                        )
)

# L3 global muons
l3GlbMuTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltPhase2L3GlbMuon"),
    cut = cms.string(""), 
    name= cms.string("l3_mu_global"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                         )
)

# L3 Muons ID (tracks)
l3MuTkIDTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltL3MuonIdTracks"),
    cut = cms.string(""), 
    name= cms.string("l3_mu_ID"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                         )
)

# Phase2 L2 seeds from L1TkMu
phase2L2SeedFromL1TkMuonTable = cms.EDProducer("SimpleL2MuonTrajectorySeedFlatTableProducer",
    src = cms.InputTag("hltL2MuonSeedsFromL1TkMuon"),
    cut = cms.string(""), 
    name= cms.string("phase2_l2_seed"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("startingState().pt()", "float", doc="p_T (GeV)"),
                         nHits = Var("nHits()", "int16", doc="number of DT/CSC segments propagated to the seed"),
                         eta = Var("l1TkMu().phEta()", "float", doc="associated L1TkMu #eta"),
                         phi = Var("l1TkMu().phPhi()", "float", doc="associated L1TkMu #phi"),
                         localX = Var("startingState().parameters().position().x()", "float", doc="local x of the seed"),
                         localY = Var("startingState().parameters().position().y()", "float", doc="local y of the seed")
                         )
)

# L2 muons to reuse (IO first)
l2MuToReuseTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("phase2L3FilteredObjects", "L2MuToReuse"),
    cut = cms.string(""), 
    name= cms.string("l2_mu_to_reuse"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                         )
)

# L3 IO tracks filtered (IO first)
l3TkIOFilteredTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("phase2L3FilteredObjects:L3IOTracksFiltered"),
    cut = cms.string(""), 
    name= cms.string("l3_tk_IO_filtered"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                         )
)

# L1 Tracker Muons to reuse (OI first)
l1TkMuToReuseTable = cms.EDProducer("SimpleTriggerL1TkMuonFlatTableProducer",
src = cms.InputTag("phase2L3FilteredObjects", "L1TkMuToReuse"),
cut = cms.string(""), 
name= cms.string("L1TkMu_to_reuse"),
doc = cms.string(""),
extension = cms.bool(False), 
variables = cms.PSet(pt = Var("phPt()", "float", doc="Physics pt"),
                     eta = Var("phEta()", "float", doc="#eta"),
                     phi = Var("phPhi()", "float", doc="#phi (rad)"),
                     dXY = Var("phD0()", "float", doc="dXY (cm)"),
                     dZ = Var("phZ0()", "float", doc="dZ (cm)")
                 )
)

# L3 OI tracks filtered (OI first)
l3TkOIFilteredTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
src = cms.InputTag("phase2L3FilteredObjects:L3OITracksFiltered"),
cut = cms.string(""), 
name= cms.string("l3_tk_OI_filtered"),
doc = cms.string(""),
extension = cms.bool(False), 
variables = cms.PSet(pt = Var("pt()", "float", doc="p_T (GeV)"),
                         eta = Var("eta()", "float", doc="#eta"),
                         phi = Var("phi()", "float", doc="#phi (rad)"),
                         dXY = Var("dxy()", "float", doc="dXY (cm)"),
                         dZ = Var("dz()", "float", doc="dZ (cm)"),
                         t0 = Var("t0()", "float", doc="t0 (ns)"),
                         nPixelHits = Var("hitPattern().numberOfValidPixelHits()", "int16", doc=""),
                         nTrkLays = Var("hitPattern().trackerLayersWithMeasurement()", "int16", doc=""),
                         nMuHits = Var("hitPattern().numberOfValidMuonHits()", "int16", doc="")
                        )
)

# The muon trigger producers sequence
muTriggerProducers = cms.Sequence(recoMuonValidationHLT_seq
                              + l1TkMuTable
                              + l2SeedTable
                              + l2SeedFromL1TkMuonTable
                              + l2MuTable
                              + l2MuTableVtx
                              + l3TkIOTable
                              + l3TkOITable
                              + l3TkMergedTable
                              + l3GlbMuTable
                              + l3MuTkIDTable
)

# The Phase-2 IO first muon trigger producers sequence
muTriggerProducersIOFirst = cms.Sequence(recoMuonValidationHLT_seq
                              + l1TkMuTable
                              + phase2L2SeedFromL1TkMuonTable
                              + l2MuTable
                              + l2MuTableVtx
                              + l3TkIOTable
                              + l2MuToReuseTable
                              + l3TkIOFilteredTable
                              + l3TkOITable
                              + l3TkMergedTable
                              + l3GlbMuTable
                              + l3MuTkIDTable
)

# The Phase-2 OI first muon trigger producers sequence
muTriggerProducersOIFirst = cms.Sequence(recoMuonValidationHLT_seq
                              + l1TkMuTable
                              + phase2L2SeedFromL1TkMuonTable
                              + l2MuTable
                              + l2MuTableVtx
                              + l3TkOITable
                              + l1TkMuToReuseTable
                              + l3TkOIFilteredTable
                              + l3TkIOTable
                              + l3TkMergedTable
                              + l3GlbMuTable
                              + l3MuTkIDTable
)

from Configuration.ProcessModifiers.phase2Muon_cff import phase2Muon, L3IOFIRST
if L3IOFIRST:
    phase2Muon.toReplaceWith(muTriggerProducers, muTriggerProducersIOFirst)
else:
    phase2Muon.toReplaceWith(muTriggerProducers, muTriggerProducersOIFirst)
