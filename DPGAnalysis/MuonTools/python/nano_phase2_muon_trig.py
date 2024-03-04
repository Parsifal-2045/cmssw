import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *


from PhysicsTools.NanoAOD.l1trig_cff import *

from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer


l1TkMuTable = cms.EDProducer("SimpleTriggerL1TkMuonFlatTableProducer",
    src = cms.InputTag("l1tTkMuonsGmt"),
    # minBX = cms.int32(-2),
    # maxBX = cms.int32(2),                           
    cut = cms.string(""), 
    name= cms.string("L1TkMu"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(l1ObjVars,
                         phPt = Var("phPt()", "float", doc="Physics pt")
                         # hwCharge = Var("hwCharge()","int16",doc="Charge (can be 0 if the charge measurement was not valid)"),
                         # hwChargeValid = Var("hwChargeValid()","int16",doc=""),
                         # tfMuonIndex = Var("tfMuonIndex()","uint16",doc="Index of muon at the uGMT input. 3 indices per link/sector/wedge. EMTF+ are 0-17, OMTF+ are 18-35, BMTF are 36-71, OMTF- are 72-89, EMTF- are 90-107"),
                         # hwTag = Var("hwTag()","int16",doc="not in L1 ntuples"),
                         # hwEtaAtVtx = Var("hwEtaAtVtx()","int16",doc="hardware eta estimated at the vertex"),
                         # hwPhiAtVtx = Var("hwPhiAtVtx()","int16",doc="hardware phi estimated at the vertex"),
                         # etaAtVtx = Var("etaAtVtx()",float,doc="eta estimated at the vertex"),
                         # phiAtVtx = Var("phiAtVtx()",float,doc="phi estimated at the vertex"),
                         # hwIsoSum = Var("hwIsoSum()","int16",doc="not in L1 ntuples"),
                         # hwDPhiExtra = Var("hwDPhiExtra()","int16",doc="Delta between Pseudo-rapidity at the muon system and the projected coordinate at the vertex in HW unit (for future l1t-integration-tag"),
                         # hwDEtaExtra = Var("hwDEtaExtra()","int16",doc="Delta between Azimuth at the muon system and the projected coordinate at the vertex in HW unit (for future l1t-integration-tag)"),
                         # hwRank = Var("hwRank()","int16",doc="not in L1Ntuples"),
                         # hwPtUnconstrained = Var("hwPtUnconstrained()","int16",doc=""),
                         # ptUnconstrained = Var("ptUnconstrained()",float,doc=""),
                         # hwDXY = Var("hwDXY()","uint16",doc=""),
                     )
)

l2MuTable = cms.EDProducer("SimpleTriggerL2TkMuonFlatTableProducer",
    src = cms.InputTag("hltL2MuonsFromL1TkMuon"),
    cut = cms.string(""), 
    name= cms.string("L2Mu"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(pt = Var("pt()", "float", doc="Physics pt")
                     )
)



muTriggerProducers = cms.Sequence(l1TkMuTable
                                  + l2MuTable
                                )
