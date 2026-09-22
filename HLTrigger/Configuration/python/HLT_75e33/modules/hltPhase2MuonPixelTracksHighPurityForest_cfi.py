import FWCore.ParameterSet.Config as cms

# Pixel-track HighPurity selector using XGBoost forest (compact .bin inference).
# Identical 33-feature extraction as the DNN selector; only the inference
# backend changes (forest tree traversal vs ONNX Runtime).
# The .bin is 2.5x smaller than the ONNX model and ~2x faster at 3-12 tracks/event.
hltPhase2MuonPixelTracksHighPurityForest = cms.EDProducer('MuonIOTracksForestSelector',
    tracks = cms.InputTag('hltPhase2MuonPixelTracks'),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    modelPath = cms.string('RecoMuon/L3TrackFinder/data/pixel_track_selector_forest.bin'),
    decisionThreshold = cms.double(0.627851665019989),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(33),
    dumpFeatures = cms.untracked.bool(False),
)
