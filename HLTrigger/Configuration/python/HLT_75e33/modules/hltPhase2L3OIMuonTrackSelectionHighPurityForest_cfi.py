import FWCore.ParameterSet.Config as cms

# OI HighPurity selection. Base: cut-based highPurity cloner (the default when
# no selector process modifier is active). Under phase2MuonPixelTracksSelector
# / phase2MuonSeedsSelector the XGBoost forests take over.
hltPhase2L3OIMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
    copyExtras = cms.untracked.bool(True),
    copyTrajectories = cms.untracked.bool(False),
    minQuality = cms.string('highPurity'),
    originalMVAVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","MVAValues"),
    originalQualVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","QualityMasks"),
    originalSource = cms.InputTag("hltPhase2L3OIMuCtfWithMaterialTracks")
)

# OI XGBoost forests (compact .bin inference).
# MuonOITracksForestSelector extracts the 22-feature OI production set via
# the shared header
# RecoMuon/L3TrackFinder/interface/OITrackSelectorFeatures.h (named struct
# muonhp::OITrackFeatures).

_pixelOIForestSelector = cms.EDProducer('MuonOITracksForestSelector',
    tracks = cms.InputTag('hltPhase2L3OIMuCtfWithMaterialTracks'),
    standaloneMuons = cms.InputTag('hltL2MuonsFromL1TkMuon', 'UpdatedAtVtx'),
    # OI pixel forest v3: 300 trees, 22 features, trained 2026-09-23 (git ad9361e5e2); working points = validation per-pT-bin
    # F2 set points from production/oi/pixel/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/OI/muonHP_OI_pixelPath_forest_v3.bin'),
    decisionThreshold = cms.double(0.27376416325569153),
    ptBinEdges = cms.vdouble(0.0, 2.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.1006140485405922,  # pT [0, 2)
        0.27376416325569153,  # pT [2, 5)
        0.276641845703125,  # pT [5, 10)
        0.3339681029319763,  # pT [10, 50)
        0.3438858985900879,  # pT [50, 200)
        0.5035520195960999,  # pT [200, inf)
    ),
    nFeatures = cms.int32(22),
    dumpFeatures = cms.untracked.bool(False),
)

_seedsOIForestSelector = _pixelOIForestSelector.clone(
    # OI general forest v3: 200 trees, 22 features, trained 2026-09-23 (git ad9361e5e2); working points = validation per-pT-bin
    # F2 set points from production/oi/general/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/OI/muonHP_OI_seedsPath_forest_v3.bin'),
    decisionThreshold = cms.double(0.5891533493995667),
    decisionThresholds = cms.vdouble(
        0.47719821333885193,  # pT [0, 2)
        0.4623990058898926,  # pT [2, 5)
        0.32259517908096313,  # pT [5, 10)
        0.6471585631370544,  # pT [10, 50)
        0.6572598814964294,  # pT [50, 200)
        0.5691716074943542,  # pT [200, inf)
    ),
)

from Configuration.ProcessModifiers.phase2MuonPixelTracksSelector_cff import phase2MuonPixelTracksSelector
from Configuration.ProcessModifiers.ngtScouting_cff import ngtScouting
(phase2MuonPixelTracksSelector | ngtScouting).toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _pixelOIForestSelector
)

from Configuration.ProcessModifiers.phase2MuonSeedsSelector_cff import phase2MuonSeedsSelector
(phase2MuonSeedsSelector).toReplaceWith(
    hltPhase2L3OIMuonTrackSelectionHighPurity,
    _seedsOIForestSelector
)
