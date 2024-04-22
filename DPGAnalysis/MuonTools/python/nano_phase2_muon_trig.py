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
    name= cms.string("l2_seed"),
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
    name= cms.string("l2_seed_froml1"),
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

# L2 cosmic offline seeds
l2CosmicSeedTable = cms.EDProducer("SimpleTrajectorySeedFlatTableProducer",
    src = cms.InputTag("hltL2CosmicOfflineMuonSeeds"),
    cut = cms.string(""), 
    name= cms.string("l2_cosmic_seed"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("startingState().pt()", "float", doc="p_T (GeV)"),
                         nHits = Var("nHits()", "int16", doc="")
                         )
)

# L2 cosmic seeds from L1 Muons
l2CosmicSeedFromL1MuonTable = cms.EDProducer("SimpleTrajectorySeedFlatTableProducer",
    src = cms.InputTag("hltL2CosmicMuonSeedsFromL1Muon"),
    cut = cms.string(""), 
    name= cms.string("l2_cosmic_seed_froml1"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("startingState().pt()", "float", doc="p_T (GeV)"),
                         nHits = Var("nHits()", "int16", doc="")
                         )
)

# L2 cosmic standalone muons
l2MuTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltL2CosmicMuonsFromL1Muon"),
    cut = cms.string(""), 
    name= cms.string("l2_cosmic_mu"),
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

# L3 OI global muons
l3OIMuTable = cms.EDProducer("SimpleTriggerTrackFlatTableProducer",
    src = cms.InputTag("hltL3MuonsPhase2L3OI"),
    cut = cms.string(""), 
    name= cms.string("l3_mu_OI"),
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
    src = cms.InputTag("hltPhase2L3MuonTracks"),
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

muTriggerProducers = cms.Sequence(l1TkMuTable
                                  + l2SeedTable
                                  + l2SeedFromL1TkMuonTable
                                  + l2CosmicSeedTable
                                  + l2CosmicSeedFromL1MuonTable
                                  + l3TkIOTable
                                  + l2MuTable
                                  + l2MuTableVtx
                                  + l3TkOITable
                                  + l3TkMergedTable
                                  + l3GlbMuTable
                                  + l3OIMuTable
                                  + l3MuTkIDTable
                                  + recoMuonValidationHLT_seq
                                )
