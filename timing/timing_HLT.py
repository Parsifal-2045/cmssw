# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: timing -s HLT:75e33 --processName=reHLT --conditions auto:phase2_realistic_T25 --geometry Extended2026D98 --era Phase2C17I13M9 --filein=input.root -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('reHLT',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HLTrigger.Configuration.HLT_75e33_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

import os
source_path = '/eos/cms/store/relval/CMSSW_14_0_0/RelValSingleMuPt100/GEN-SIM-RECO/PU_140X_mcRun4_realistic_v1_STD_2026D98_PU-v2/2580000/'
files = []
for path in os.listdir(source_path):
    if os.path.isfile(os.path.join(source_path, path)):
        files.append('file:' + os.path.join(source_path, path))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(files),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('timing nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string(''),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('timing_HLT.root'),
#     outputCommands = process.RECOSIMEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T25', '')

# Path and EndPath definitions
# process.HLT_AK4PFPuppiJet520 = cms.Path(process.HLTBeginSequence+process.l1tSinglePFPuppiJet230off+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltSingleAK4PFPuppiJet520+process.HLTEndSequence)
# process.HLT_Diphoton30_23_IsoCaloId_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTDiphoton3023IsoCaloIdL1SeededSequence+process.HLTEndSequence)
# process.HLT_Diphoton30_23_IsoCaloId_Unseeded = cms.Path(process.HLTBeginSequence+process.HLTDiphoton3023IsoCaloIdUnseededSequence+process.HLTEndSequence)
# process.HLT_DoubleEle23_12_Iso_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTDoubleEle2312IsoL1SeededSequence+process.HLTEndSequence)
# process.HLT_DoubleEle25_CaloIdL_PMS2_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTDoubleEle25CaloIdLPMS2L1SeededSequence+process.HLTEndSequence)
# process.HLT_DoubleEle25_CaloIdL_PMS2_Unseeded = cms.Path(process.HLTBeginSequence+process.HLTDoubleEle25CaloIdLPMS2UnseededSequence+process.HLTEndSequence)
# process.HLT_DoubleMediumChargedIsoPFTauHPS40_eta2p1 = cms.Path(process.HLTBeginSequence+process.hltPreDoublePFTauHPS+process.HLTParticleFlowSequence+process.HLTAK4PFJetsReconstruction+process.hltAK4PFJetsForTaus+process.HLTPFTauHPS+process.HLTHPSMediumChargedIsoPFTauSequence+process.hltHpsSelectedPFTausTrackPt1MediumChargedIsolation+process.hltHpsDoublePFTau40TrackPt1MediumChargedIsolation+process.HLTEndSequence)
# process.HLT_DoubleMediumDeepTauPFTauHPS35_eta2p1 = cms.Path(process.HLTBeginSequence+process.hltPreDoublePFTauHPS+process.HLTParticleFlowSequence+process.HLTAK4PFJetsReconstruction+process.hltAK4PFJetsForTaus+process.HLTPFTauHPS+process.HLTHPSDeepTauPFTauSequence+process.hltHpsSelectedPFTausMediumDitauWPDeepTau+process.hltHpsDoublePFTau35MediumDitauWPDeepTau+process.HLTEndSequence)
# process.HLT_DoublePFPuppiJets128_DoublePFPuppiBTagDeepCSV_2p4 = cms.Path(process.HLTBeginSequence+process.l1tDoublePFPuppiJet112offMaxEta2p4+process.l1tDoublePFPuppiJets112offMaxDeta1p6+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltDoublePFPuppiJets128MaxEta2p4+process.hltDoublePFPuppiJets128Eta2p4MaxDeta1p6+process.HLTBtagDeepCSVSequencePFPuppiModEta2p4+process.hltBTagPFPuppiDeepCSV0p865DoubleEta2p4+process.HLTEndSequence)
# process.HLT_DoublePFPuppiJets128_DoublePFPuppiBTagDeepFlavour_2p4 = cms.Path(process.HLTBeginSequence+process.l1tDoublePFPuppiJet112offMaxEta2p4+process.l1tDoublePFPuppiJets112offMaxDeta1p6+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltDoublePFPuppiJets128MaxEta2p4+process.hltDoublePFPuppiJets128Eta2p4MaxDeta1p6+process.HLTBtagDeepFlavourSequencePFPuppiModEta2p4+process.hltBTagPFPuppiDeepFlavour0p935DoubleEta2p4+process.HLTEndSequence)
# process.HLT_Ele115_NonIso_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTEle115NonIsoL1SeededSequence+process.HLTEndSequence)
# process.HLT_Ele26_WP70_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTEle26WP70L1SeededSequence+process.HLTEndSequence)
# process.HLT_Ele26_WP70_Unseeded = cms.Path(process.HLTBeginSequence+process.HLTEle26WP70UnseededSequence+process.HLTEndSequence)
# process.HLT_Ele32_WPTight_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTEle32WPTightL1SeededSequence+process.HLTEndSequence)
# process.HLT_Ele32_WPTight_Unseeded = cms.Path(process.HLTBeginSequence+process.HLTEle32WPTightUnseededSequence+process.HLTEndSequence)
# process.HLT_IsoMu24_FromL1TkMuon = cms.Path(process.HLTBeginSequence+process.hltL3fL1TkSingleMu22L3Filtered24Q+process.hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41+process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40+process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70+process.hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk+process.HLTEndSequence, cms.Task(process.HGCalRecHit, process.HGCalUncalibRecHit, process.MeasurementTrackerEvent, process.bunchSpacingProducer, process.hgcalDigis, process.hgcalLayerClustersEE, process.hgcalLayerClustersHSci, process.hgcalLayerClustersHSi, process.hgcalMergeLayerClusters, process.hltCsc2DRecHits, process.hltCscSegments, process.hltDt1DRecHits, process.hltDt4DSegments, process.hltEcalDetIdToBeRecovered, process.hltEcalDigis, process.hltEcalRecHit, process.hltEcalUncalibRecHit, process.hltFixedGridRhoFastjetAllCaloForEGamma, process.hltGemRecHits, process.hltGemSegments, process.hltHbhereco, process.hltHcalDigis, process.hltHfprereco, process.hltHfreco, process.hltHoreco, process.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks, process.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter2Phase2L3FromL1TkMuonClustersRefRemoval, process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent, process.hltIter2Phase2L3FromL1TkMuonMerged, process.hltIter2Phase2L3FromL1TkMuonPixelClusterCheck, process.hltIter2Phase2L3FromL1TkMuonPixelHitDoublets, process.hltIter2Phase2L3FromL1TkMuonPixelHitTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelSeeds, process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered, process.hltIter2Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltL2MuonFromL1TkMuonCandidates, process.hltL2MuonSeedsFromL1TkMuon, process.hltL2MuonsFromL1TkMuon, process.hltL2OfflineMuonSeeds, process.hltL3MuonsPhase2L3Links, process.hltL3MuonsPhase2L3OI, process.hltMe0RecHits, process.hltMe0Segments, process.hltParticleFlowClusterECALUncorrectedUnseeded, process.hltParticleFlowClusterECALUnseeded, process.hltParticleFlowClusterHBHE, process.hltParticleFlowClusterHCAL, process.hltParticleFlowRecHitECALUnseeded, process.hltParticleFlowRecHitHBHE, process.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracks, process.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets, process.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions, process.hltPhase2L3FromL1TkMuonPixelVertices, process.hltPhase2L3FromL1TkMuonTrimmedPixelVertices, process.hltPhase2L3GlbMuon, process.hltPhase2L3MuonCandidates, process.hltPhase2L3MuonGeneralTracks, process.hltPhase2L3MuonHighPtTripletStepClusters, process.hltPhase2L3MuonHighPtTripletStepHitDoublets, process.hltPhase2L3MuonHighPtTripletStepHitTriplets, process.hltPhase2L3MuonHighPtTripletStepSeedLayers, process.hltPhase2L3MuonHighPtTripletStepSeeds, process.hltPhase2L3MuonHighPtTripletStepTrackCandidates, process.hltPhase2L3MuonHighPtTripletStepTrackCutClassifier, process.hltPhase2L3MuonHighPtTripletStepTrackingRegions, process.hltPhase2L3MuonHighPtTripletStepTracks, process.hltPhase2L3MuonHighPtTripletStepTracksSelectionHighPurity, process.hltPhase2L3MuonInitialStepSeeds, process.hltPhase2L3MuonInitialStepTrackCandidates, process.hltPhase2L3MuonInitialStepTrackCutClassifier, process.hltPhase2L3MuonInitialStepTracks, process.hltPhase2L3MuonInitialStepTracksSelectionHighPurity, process.hltPhase2L3MuonMerged, process.hltPhase2L3MuonPixelTracks, process.hltPhase2L3MuonPixelTracksHitDoublets, process.hltPhase2L3MuonPixelTracksHitQuadruplets, process.hltPhase2L3MuonPixelTracksSeedLayers, process.hltPhase2L3MuonPixelTracksTrackingRegions, process.hltPhase2L3MuonPixelVertices, process.hltPhase2L3MuonTracks, process.hltPhase2L3Muons, process.hltPhase2L3MuonsEcalIsodR0p3dRVeto0p000, process.hltPhase2L3MuonsHcalIsodR0p3dRVeto0p000, process.hltPhase2L3MuonsHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00, process.hltPhase2L3MuonsNoID, process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07, process.hltPhase2L3OIL3MuonCandidates, process.hltPhase2L3OIL3Muons, process.hltPhase2L3OIL3MuonsLinksCombination, process.hltPhase2L3OIMuCtfWithMaterialTracks, process.hltPhase2L3OIMuonTrackCutClassifier, process.hltPhase2L3OIMuonTrackSelectionHighPurity, process.hltPhase2L3OISeedsFromL2Muons, process.hltPhase2L3OITrackCandidates, process.hltPhase2PixelFitterByHelixProjections, process.hltPhase2PixelTrackFilterByKinematics, process.hltPhase2TowerMakerForAll, process.hltRpcRecHits, process.siPhase2Clusters, process.siPixelClusterShapeCache, process.siPixelClusters, process.siPixelRecHits, process.trackerClusterCheck))
# process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon = cms.Path(process.HLTBeginSequence+process.hltL1TkDoubleMuFiltered7+process.hltL1TkSingleMuFiltered15+process.hltDoubleMuon7DZ1p0+process.hltL3fL1DoubleMu155fPreFiltered8+process.hltL3fL1DoubleMu155fFiltered17+process.hltDiMuon178RelTrkIsoFiltered0p4+process.hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2+process.HLTEndSequence, cms.Task(process.MeasurementTrackerEvent, process.hltCsc2DRecHits, process.hltCscSegments, process.hltDt1DRecHits, process.hltDt4DSegments, process.hltGemRecHits, process.hltGemSegments, process.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks, process.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter2Phase2L3FromL1TkMuonClustersRefRemoval, process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent, process.hltIter2Phase2L3FromL1TkMuonMerged, process.hltIter2Phase2L3FromL1TkMuonPixelClusterCheck, process.hltIter2Phase2L3FromL1TkMuonPixelHitDoublets, process.hltIter2Phase2L3FromL1TkMuonPixelHitTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelSeeds, process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered, process.hltIter2Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltL2MuonFromL1TkMuonCandidates, process.hltL2MuonSeedsFromL1TkMuon, process.hltL2MuonsFromL1TkMuon, process.hltL2OfflineMuonSeeds, process.hltL3MuonsPhase2L3Links, process.hltL3MuonsPhase2L3OI, process.hltMe0RecHits, process.hltMe0Segments, process.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracks, process.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets, process.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions, process.hltPhase2L3FromL1TkMuonPixelVertices, process.hltPhase2L3FromL1TkMuonTrimmedPixelVertices, process.hltPhase2L3GlbMuon, process.hltPhase2L3MuonCandidates, process.hltPhase2L3MuonGeneralTracks, process.hltPhase2L3MuonHighPtTripletStepClusters, process.hltPhase2L3MuonHighPtTripletStepHitDoublets, process.hltPhase2L3MuonHighPtTripletStepHitTriplets, process.hltPhase2L3MuonHighPtTripletStepSeedLayers, process.hltPhase2L3MuonHighPtTripletStepSeeds, process.hltPhase2L3MuonHighPtTripletStepTrackCandidates, process.hltPhase2L3MuonHighPtTripletStepTrackCutClassifier, process.hltPhase2L3MuonHighPtTripletStepTrackingRegions, process.hltPhase2L3MuonHighPtTripletStepTracks, process.hltPhase2L3MuonHighPtTripletStepTracksSelectionHighPurity, process.hltPhase2L3MuonInitialStepSeeds, process.hltPhase2L3MuonInitialStepTrackCandidates, process.hltPhase2L3MuonInitialStepTrackCutClassifier, process.hltPhase2L3MuonInitialStepTracks, process.hltPhase2L3MuonInitialStepTracksSelectionHighPurity, process.hltPhase2L3MuonMerged, process.hltPhase2L3MuonPixelTracks, process.hltPhase2L3MuonPixelTracksHitDoublets, process.hltPhase2L3MuonPixelTracksHitQuadruplets, process.hltPhase2L3MuonPixelTracksSeedLayers, process.hltPhase2L3MuonPixelTracksTrackingRegions, process.hltPhase2L3MuonPixelVertices, process.hltPhase2L3MuonTracks, process.hltPhase2L3Muons, process.hltPhase2L3MuonsNoID, process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p4, process.hltPhase2L3OIL3MuonCandidates, process.hltPhase2L3OIL3Muons, process.hltPhase2L3OIL3MuonsLinksCombination, process.hltPhase2L3OIMuCtfWithMaterialTracks, process.hltPhase2L3OIMuonTrackCutClassifier, process.hltPhase2L3OIMuonTrackSelectionHighPurity, process.hltPhase2L3OISeedsFromL2Muons, process.hltPhase2L3OITrackCandidates, process.hltPhase2PixelFitterByHelixProjections, process.hltPhase2PixelTrackFilterByKinematics, process.hltRpcRecHits, process.siPhase2Clusters, process.siPixelClusterShapeCache, process.siPixelClusters, process.siPixelRecHits, process.trackerClusterCheck))
# process.HLT_Mu37_Mu27_FromL1TkMuon = cms.Path(process.HLTBeginSequence+process.hltDoubleMuon7DZ1p0+process.hltL3fL1DoubleMu155fPreFiltered27+process.hltL3fL1DoubleMu155fFiltered37+process.HLTEndSequence, cms.Task(process.MeasurementTrackerEvent, process.hltCsc2DRecHits, process.hltCscSegments, process.hltDt1DRecHits, process.hltDt4DSegments, process.hltGemRecHits, process.hltGemSegments, process.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks, process.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter2Phase2L3FromL1TkMuonClustersRefRemoval, process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent, process.hltIter2Phase2L3FromL1TkMuonMerged, process.hltIter2Phase2L3FromL1TkMuonPixelClusterCheck, process.hltIter2Phase2L3FromL1TkMuonPixelHitDoublets, process.hltIter2Phase2L3FromL1TkMuonPixelHitTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelSeeds, process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered, process.hltIter2Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltL2MuonFromL1TkMuonCandidates, process.hltL2MuonSeedsFromL1TkMuon, process.hltL2MuonsFromL1TkMuon, process.hltL2OfflineMuonSeeds, process.hltL3MuonsPhase2L3Links, process.hltL3MuonsPhase2L3OI, process.hltMe0RecHits, process.hltMe0Segments, process.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracks, process.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets, process.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions, process.hltPhase2L3FromL1TkMuonPixelVertices, process.hltPhase2L3FromL1TkMuonTrimmedPixelVertices, process.hltPhase2L3GlbMuon, process.hltPhase2L3MuonCandidates, process.hltPhase2L3MuonMerged, process.hltPhase2L3MuonTracks, process.hltPhase2L3Muons, process.hltPhase2L3MuonsNoID, process.hltPhase2L3OIL3MuonCandidates, process.hltPhase2L3OIL3Muons, process.hltPhase2L3OIL3MuonsLinksCombination, process.hltPhase2L3OIMuCtfWithMaterialTracks, process.hltPhase2L3OIMuonTrackCutClassifier, process.hltPhase2L3OIMuonTrackSelectionHighPurity, process.hltPhase2L3OISeedsFromL2Muons, process.hltPhase2L3OITrackCandidates, process.hltPhase2PixelFitterByHelixProjections, process.hltPhase2PixelTrackFilterByKinematics, process.hltRpcRecHits, process.siPhase2Clusters, process.siPixelClusterShapeCache, process.siPixelClusters, process.siPixelRecHits))
# process.HLT_PFHT200PT30_QuadPFPuppiJet_70_40_30_30_TriplePFPuppiBTagDeepCSV_2p4 = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiHT400offMaxEta2p4+process.l1t1PFPuppiJet70offMaxEta2p4+process.l1t2PFPuppiJet55offMaxEta2p4+process.l1t4PFPuppiJet40offMaxEta2p4+process.l1t4PFPuppiJet25OnlineMaxEta2p4+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltPFPuppiCentralJetQuad30MaxEta2p4+process.hlt1PFPuppiCentralJet70MaxEta2p4+process.hlt2PFPuppiCentralJet40MaxEta2p4+process.hltHtMhtPFPuppiCentralJetsQuadC30MaxEta2p4+process.hltPFPuppiCentralJetsQuad30HT200MaxEta2p4+process.HLTBtagDeepCSVSequencePFPuppiModEta2p4+process.hltBTagPFPuppiDeepCSV0p38Eta2p4TripleEta2p4+process.HLTEndSequence)
# process.HLT_PFHT200PT30_QuadPFPuppiJet_70_40_30_30_TriplePFPuppiBTagDeepFlavour_2p4 = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiHT400offMaxEta2p4+process.l1t1PFPuppiJet70offMaxEta2p4+process.l1t2PFPuppiJet55offMaxEta2p4+process.l1t4PFPuppiJet40offMaxEta2p4+process.l1t4PFPuppiJet25OnlineMaxEta2p4+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltPFPuppiCentralJetQuad30MaxEta2p4+process.hlt1PFPuppiCentralJet70MaxEta2p4+process.hlt2PFPuppiCentralJet40MaxEta2p4+process.hltHtMhtPFPuppiCentralJetsQuadC30MaxEta2p4+process.hltPFPuppiCentralJetsQuad30HT200MaxEta2p4+process.HLTBtagDeepFlavourSequencePFPuppiModEta2p4+process.hltBTagPFPuppiDeepFlavour0p375Eta2p4TripleEta2p4+process.HLTEndSequence)
# process.HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepCSV_2p4 = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiHT400offMaxEta2p4+process.l1t1PFPuppiJet70offMaxEta2p4+process.l1t2PFPuppiJet55offMaxEta2p4+process.l1t4PFPuppiJet40offMaxEta2p4+process.l1t4PFPuppiJet25OnlineMaxEta2p4+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltPFPuppiCentralJetQuad30MaxEta2p4+process.hlt1PFPuppiCentralJet75MaxEta2p4+process.hlt2PFPuppiCentralJet60MaxEta2p4+process.hlt3PFPuppiCentralJet45MaxEta2p4+process.hlt4PFPuppiCentralJet40MaxEta2p4+process.hltHtMhtPFPuppiCentralJetsQuadC30MaxEta2p4+process.hltPFPuppiCentralJetsQuad30HT330MaxEta2p4+process.HLTBtagDeepCSVSequencePFPuppiModEta2p4+process.hltBTagPFPuppiDeepCSV0p31Eta2p4TripleEta2p4+process.HLTEndSequence)
# process.HLT_PFHT330PT30_QuadPFPuppiJet_75_60_45_40_TriplePFPuppiBTagDeepFlavour_2p4 = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiHT400offMaxEta2p4+process.l1t1PFPuppiJet70offMaxEta2p4+process.l1t2PFPuppiJet55offMaxEta2p4+process.l1t4PFPuppiJet40offMaxEta2p4+process.l1t4PFPuppiJet25OnlineMaxEta2p4+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltPFPuppiCentralJetQuad30MaxEta2p4+process.hlt1PFPuppiCentralJet75MaxEta2p4+process.hlt2PFPuppiCentralJet60MaxEta2p4+process.hlt3PFPuppiCentralJet45MaxEta2p4+process.hlt4PFPuppiCentralJet40MaxEta2p4+process.hltHtMhtPFPuppiCentralJetsQuadC30MaxEta2p4+process.hltPFPuppiCentralJetsQuad30HT330MaxEta2p4+process.HLTBtagDeepFlavourSequencePFPuppiModEta2p4+process.hltBTagPFPuppiDeepFlavour0p275Eta2p4TripleEta2p4+process.HLTEndSequence)
# process.HLT_PFPuppiHT1070 = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiHT450off+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.hltPFPuppiHT+process.hltPFPuppiHT1070+process.HLTEndSequence)
# process.HLT_PFPuppiMETTypeOne140_PFPuppiMHT140 = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiMET220off+process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.HLTPFPuppiMETReconstruction+process.hltPFPuppiMETTypeOneCorrector+process.hltPFPuppiMETTypeOne+process.hltPFPuppiMETTypeOne140+process.hltPFPuppiMHT+process.hltPFPuppiMHT140+process.HLTEndSequence)
# process.HLT_Photon108EB_TightID_TightIso_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTPhoton108EBTightIDTightIsoL1SeededSequence+process.HLTEndSequence)
# process.HLT_Photon108EB_TightID_TightIso_Unseeded = cms.Path(process.HLTBeginSequence+process.HLTPhoton108EBTightIDTightIsoUnseededSequence+process.HLTEndSequence)
# process.HLT_Photon187_L1Seeded = cms.Path(process.HLTBeginSequence+process.HLTPhoton187L1SeededSequence+process.HLTEndSequence)
# process.HLT_Photon187_Unseeded = cms.Path(process.HLTBeginSequence+process.HLTPhoton187UnseededSequence+process.HLTEndSequence)
# process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon = cms.Path(process.HLTBeginSequence+process.hltTripleMuon3DZ1p0+process.hltTripleMuon3DR0+process.hltL3fL1TkTripleMu533PreFiltered555+process.hltL3fL1TkTripleMu533L3Filtered1055+process.hltL3fL1TkTripleMu533L31055DZFiltered0p2+process.HLTEndSequence, cms.Task(process.MeasurementTrackerEvent, process.hltCsc2DRecHits, process.hltCscSegments, process.hltDt1DRecHits, process.hltDt4DSegments, process.hltGemRecHits, process.hltGemSegments, process.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks, process.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates, process.hltIter2Phase2L3FromL1TkMuonClustersRefRemoval, process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks, process.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent, process.hltIter2Phase2L3FromL1TkMuonMerged, process.hltIter2Phase2L3FromL1TkMuonPixelClusterCheck, process.hltIter2Phase2L3FromL1TkMuonPixelHitDoublets, process.hltIter2Phase2L3FromL1TkMuonPixelHitTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets, process.hltIter2Phase2L3FromL1TkMuonPixelSeeds, process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered, process.hltIter2Phase2L3FromL1TkMuonTrackCutClassifier, process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity, process.hltL2MuonFromL1TkMuonCandidates, process.hltL2MuonSeedsFromL1TkMuon, process.hltL2MuonsFromL1TkMuon, process.hltL2OfflineMuonSeeds, process.hltL3MuonsPhase2L3Links, process.hltL3MuonsPhase2L3OI, process.hltMe0RecHits, process.hltMe0Segments, process.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracks, process.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets, process.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets, process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions, process.hltPhase2L3FromL1TkMuonPixelVertices, process.hltPhase2L3FromL1TkMuonTrimmedPixelVertices, process.hltPhase2L3GlbMuon, process.hltPhase2L3MuonCandidates, process.hltPhase2L3MuonMerged, process.hltPhase2L3MuonTracks, process.hltPhase2L3Muons, process.hltPhase2L3MuonsNoID, process.hltPhase2L3OIL3MuonCandidates, process.hltPhase2L3OIL3Muons, process.hltPhase2L3OIL3MuonsLinksCombination, process.hltPhase2L3OIMuCtfWithMaterialTracks, process.hltPhase2L3OIMuonTrackCutClassifier, process.hltPhase2L3OIMuonTrackSelectionHighPurity, process.hltPhase2L3OISeedsFromL2Muons, process.hltPhase2L3OITrackCandidates, process.hltPhase2PixelFitterByHelixProjections, process.hltPhase2PixelTrackFilterByKinematics, process.hltRpcRecHits, process.siPhase2Clusters, process.siPixelClusterShapeCache, process.siPixelClusters, process.siPixelRecHits))
# process.L1T_DoubleNNTau52 = cms.Path(process.HLTL1Sequence+process.hltL1DoubleNNTau52)
# process.L1T_DoublePFPuppiJets112_2p4_DEta1p6 = cms.Path(process.HLTBeginSequence+process.l1tDoublePFPuppiJet112offMaxEta2p4+process.l1tDoublePFPuppiJets112offMaxDeta1p6+process.HLTEndSequence)
# process.L1T_PFHT400PT30_QuadPFPuppiJet_70_55_40_40_2p4 = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiHT400offMaxEta2p4+process.l1t1PFPuppiJet70offMaxEta2p4+process.l1t2PFPuppiJet55offMaxEta2p4+process.l1t4PFPuppiJet40offMaxEta2p4+process.l1t4PFPuppiJet25OnlineMaxEta2p4+process.HLTEndSequence)
# process.L1T_PFPuppiHT450off = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiHT450off+process.HLTEndSequence)
# process.L1T_PFPuppiMET220off = cms.Path(process.HLTBeginSequence+process.l1tPFPuppiMET220off+process.HLTEndSequence)
# process.L1T_SingleNNTau150 = cms.Path(process.HLTL1Sequence+process.hltL1SingleNNTau150)
# process.L1T_SinglePFPuppiJet230off = cms.Path(process.HLTBeginSequence+process.l1tSinglePFPuppiJet230off+process.HLTEndSequence)
# process.L1T_TkEle25TkEle12 = cms.Path(process.HLTBeginSequence+process.L1TTkEle25TkEle12Sequence+process.HLTEndSequence)
# process.L1T_TkEle36 = cms.Path(process.HLTBeginSequence+process.L1TTkEle36Sequence+process.HLTEndSequence)
# process.L1T_TkEm37TkEm24 = cms.Path(process.HLTBeginSequence+process.L1TTkEm37TkEm24Sequence+process.HLTEndSequence)
# process.L1T_TkEm51 = cms.Path(process.HLTBeginSequence+process.L1TTkEm51Sequence+process.HLTEndSequence)
# process.L1T_TkIsoEle22TkEm12 = cms.Path(process.HLTBeginSequence+process.L1TTkIsoEle22TkEm12Sequence+process.HLTEndSequence)
# process.L1T_TkIsoEle28 = cms.Path(process.HLTBeginSequence+process.L1TTkIsoEle28Sequence+process.HLTEndSequence)
# process.L1T_TkIsoEm22TkIsoEm12 = cms.Path(process.HLTBeginSequence+process.L1TTkIsoEm22TkIsoEm12Sequence+process.HLTEndSequence)
# process.L1T_TkIsoEm36 = cms.Path(process.HLTBeginSequence+process.L1TTkIsoEm36Sequence+process.HLTEndSequence)
# process.MC_BTV = cms.Path(process.HLTParticleFlowSequence+process.HLTAK4PFPuppiJetsReconstruction+process.HLTBtagDeepCSVSequencePFPuppi+process.HLTBtagDeepFlavourSequencePFPuppi)
# process.MC_Ele5_Open_L1Seeded = cms.Path(process.HLTBeginSequence+process.hltPreEle5OpenL1Seeded+process.HLTEle5OpenL1SeededSequence+process.HLTEndSequence)
# process.MC_Ele5_Open_Unseeded = cms.Path(process.HLTBeginSequence+process.hltPreEle5OpenUnseeded+process.HLTEle5OpenUnseededSequence+process.HLTEndSequence)
# process.MC_JME = cms.Path(process.HLTParticleFlowSequence+process.HLTJMESequence+process.hltPFPuppiHT+process.hltPFPuppiMHT)
process.endjob_step = cms.EndPath(process.endOfProcess)
# process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)


process.load("HLTrigger/Configuration/HLT_75e33/services/FastTimerService_cfi")
process.FastTimerService.enableDQM = False
process.FastTimerService.enableDQMbyModule = False
process.FastTimerService.enableDQMbyPath = False
process.FastTimerService.jsonFileName = 'full_timing.json'
process.FastTimerService.writeJSONSummary = cms.untracked.bool(True)

# Schedule definition
# process.schedule imported from cff in HLTrigger.Configuration
process.schedule.extend([process.endjob_step])#,process.RECOSIMoutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
