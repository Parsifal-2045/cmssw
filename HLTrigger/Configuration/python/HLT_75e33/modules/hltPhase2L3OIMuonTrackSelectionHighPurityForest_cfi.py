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
    # OI pixel forest v1: 200 trees, 22 features, trained 2026-10-10 (git e622949b3f)
    # working points = per-pT-bin recall-target set points of the validation split;
    # from production/oi/pixel/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/OI/muonHP_OI_pixelPath_forest_v1.bin'),
    decisionThreshold = cms.double(0.05863716080784798),
    ptBinEdges = cms.vdouble(0.0, 3.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.05613938346505165,  # pT [0, 3)
        0.045126959681510925,  # pT [3, 5)
        0.026730574667453766,  # pT [5, 10)
        0.05283714830875397,  # pT [10, 50)
        0.10757311433553696,  # pT [50, 200)
        0.6577504277229309,  # pT [200, inf)
    ),
    nFeatures = cms.int32(22),
    dumpFeatures = cms.untracked.bool(False),
)

_seedsOIForestSelector = _pixelOIForestSelector.clone(
    # OI general forest v1: 200 trees, 22 features, trained 2026-10-10 (git e622949b3f)
    # working points = per-pT-bin recall-target set points of the validation split;
    # from production/oi/general/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/OI/muonHP_OI_seedsPath_forest_v1.bin'),
    decisionThreshold = cms.double(0.07081226259469986),
    decisionThresholds = cms.vdouble(
        0.06444144994020462,  # pT [0, 3)
        0.03170956298708916,  # pT [3, 5)
        0.07118859142065048,  # pT [5, 10)
        0.03228765353560448,  # pT [10, 50)
        0.1375492364168167,  # pT [50, 200)
        0.8149058222770691,  # pT [200, inf)
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
