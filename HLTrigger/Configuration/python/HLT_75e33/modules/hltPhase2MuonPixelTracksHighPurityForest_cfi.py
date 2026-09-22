import FWCore.ParameterSet.Config as cms

# Pixel-track HighPurity selector using XGBoost forest (compact .bin inference).
# Identical 33-feature extraction as the DNN selector; only the inference
# backend changes (forest tree traversal vs ONNX Runtime).
# The .bin is 2.5x smaller than the ONNX model and ~2x faster at 3-12 tracks/event.
hltPhase2MuonPixelTracksHighPurityForest = cms.EDProducer('MuonIOTracksForestSelector',
    tracks = cms.InputTag('hltPhase2MuonPixelTracks'),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    modelPath = cms.string('RecoMuon/L3TrackFinder/data/pixel_track_selector_forest.bin'),
    # v2 model (evt10 split, 5000-tree budget pruned to 3250 on validation);
    # per-pT-bin F2 working points as set by the training campaign.
    # decisionThreshold (global F2) is the fallback used when ptBinEdges is empty.
    decisionThreshold = cms.double(0.5552404522895813),
    ptBinEdges = cms.vdouble(0.0, 2.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.34674736857414246,   # pT [0, 2)
        0.5552404522895813,    # pT [2, 5)
        0.7297627925872803,    # pT [5, 10)
        0.8456155061721802,    # pT [10, 50)
        0.9415441155433655,    # pT [50, 200)
        0.5552404522895813,    # pT [200, inf) - global fallback (no background in training)
    ),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(33),
    dumpFeatures = cms.untracked.bool(False),
)
