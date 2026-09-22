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
    modelPath = cms.string('RecoMuon/L3TrackFinder/data/OI_pixel_track_selector_forest.bin'),
    # OI-pixel forest (evt10 split, 5000-tree budget pruned to 1200 on
    # validation); per-pT-bin F2 working points from its training records.
    decisionThreshold = cms.double(0.4347943663597107),
    ptBinEdges = cms.vdouble(0.0, 2.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.36963194608688354,   # pT [0, 2)
        0.19431737065315247,   # pT [2, 5)
        0.4301176071166992,    # pT [5, 10)
        0.2698938846588135,    # pT [10, 50)
        0.47758954763412476,   # pT [50, 200)
        0.5737308859825134,    # pT [200, inf)
    ),
    useStandaloneMuonFeatures = cms.bool(True),
    nFeatures = cms.int32(22),
    dumpFeatures = cms.untracked.bool(False),
)

_seedsOIForestSelector = _pixelOIForestSelector.clone(
    modelPath = cms.string('RecoMuon/L3TrackFinder/data/OI_general_track_selector_forest.bin'),
    # OI-general forest (850 trees).
    decisionThreshold = cms.double(0.6497626304626465),
    decisionThresholds = cms.vdouble(
        0.42146751284599304,   # pT [0, 2)
        0.5860887765884399,    # pT [2, 5)
        0.45810478925704956,   # pT [5, 10)
        0.705085277557373,     # pT [10, 50)
        0.6939941644668579,    # pT [50, 200)
        0.6880988478660583,    # pT [200, inf)
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
