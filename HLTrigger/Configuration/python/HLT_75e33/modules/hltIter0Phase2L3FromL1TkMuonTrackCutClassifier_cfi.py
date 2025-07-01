import FWCore.ParameterSet.Config as cms

hltIter0Phase2L3FromL1TkMuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
    beamspot = cms.InputTag("hltOnlineBeamSpot"),
    ignoreVertices = cms.bool(False),
    mva = cms.PSet(
        dr_par = cms.PSet(
            d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
            d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
            dr_exp = cms.vint32(4, 4, 2147483647),
            dr_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dr_par2 = cms.vdouble(0.3, 0.3, 3.40282346639e+38)
        ),
        dz_par = cms.PSet(
            dz_exp = cms.vint32(4, 4, 2147483647),
            dz_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
            dz_par2 = cms.vdouble(0.35, 0.35, 3.40282346639e+38)
        ),
        maxChi2 = cms.vdouble(3.40282346639e+38, 3.40282346639e+38, 3.40282346639e+38),
        maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
        maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
        maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
        maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
        maxLostLayers = cms.vint32(1, 1, 1),
        min3DLayers = cms.vint32(0, 3, 4),
        minLayers = cms.vint32(3, 3, 4),
        minNVtxTrk = cms.int32(3),
        minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
        minPixelHits = cms.vint32(0, 3, 4)
    ),
    qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
    src = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks"),
    vertices = cms.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices")
)

from Configuration.ProcessModifiers.phase2CAExtension_cff import phase2CAExtension
phase2CAExtension.toModify(
    hltIter0Phase2L3FromL1TkMuonTrackCutClassifier,
    ignoreVertices = True,
    mva = dict(
        dr_par = dict(
            d0err = [1.610380103264296148, 0.502722713333466853, 1.111607948895815268],
            d0err_par = [64.742912303519062789, 81.012436524573445240, 54.385119178719129707],
            dr_exp = [19, 10, 13],
            dr_par1 = [57.059481459985178731, 50.279776876120820361, 33.988307002383727706],
            dr_par2 = [74.068635810153196530, 34.471315687253465398, 45.226909827846640155]
        ),
        dz_par = dict(
            dz_exp = [10, 6, 9],
            dz_par1 = [45.512594312811913255, 53.117027806739699258, 82.032029513077333149],
            dz_par2 = [20.052051082517227343, 26.960171536570406658, 74.717977412227568834]
        ), 
    maxChi2 = [75.613179075702618093, 20.149173982721297449, 24.655685312520887464],
    maxChi2n = [75.613179075702618093, 20.149173982721297449, 24.655685312520887464],
    maxDr = [0.632564511692216502, 0.899389163528369062, 1.021550664333019887],
    maxDz = [55.531736492650743742, 4.605282832364866685, 35.255653778578881941],
    maxDzWrtBS = [65.655403829487170242, 62.902571518998129818, 71.481381340634541743],
    maxLostLayers = [6, 7, 4],
    min3DLayers = [4, 4, 4],
    minLayers = [3, 5, 8],
    minNVtxTrk = 5,
    minNdof = [0.789722702402609444, 1.496971638145484285, 1.666093571637450532],
    minPixelHits = [5, 1, 1]
    ),
    qualityCuts = [-0.023233126663624204, 0.132977980206018087, 0.490682422036849974],
    vertices = ""
)
