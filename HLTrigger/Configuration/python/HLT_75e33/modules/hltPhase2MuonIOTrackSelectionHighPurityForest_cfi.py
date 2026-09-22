import FWCore.ParameterSet.Config as cms

# IO-track (seeds) HighPurity selector using XGBoost forest (compact .bin inference).
# Identical 33-feature extraction as the DNN selector; only the inference
# backend changes (forest tree traversal vs ONNX Runtime).
# The .bin is 2.5x smaller than the ONNX model and ~2x faster at 3-12 tracks/event.
hltPhase2MuonIOTrackSelectionHighPurityForest = cms.EDProducer("MuonIOTracksForestSelector",
    tracks = cms.InputTag("hltPhase2MuonIOTracks"),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/seeds_track_selector_forest.bin"),
    # v2 model (evt10 split, 5000-tree budget pruned to 4050 on validation);
    # per-pT-bin F2 working points as set by the training campaign.
    # decisionThreshold (global F2) is the fallback used when ptBinEdges is empty.
    decisionThreshold = cms.double(0.631809413433075),
    ptBinEdges = cms.vdouble(0.0, 2.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.42168208956718445,   # pT [0, 2)
        0.6310316324234009,    # pT [2, 5)
        0.665044367313385,     # pT [5, 10)
        0.7129977941513062,    # pT [10, 50)
        0.6601066589355469,    # pT [50, 200)
        0.9321564435958862,    # pT [200, inf)
    ),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(33),
    dumpFeatures = cms.untracked.bool(False),
)
