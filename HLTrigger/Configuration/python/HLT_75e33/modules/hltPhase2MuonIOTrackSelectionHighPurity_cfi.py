import FWCore.ParameterSet.Config as cms

hltPhase2MuonIOTrackSelectionHighPurity = cms.EDProducer("MuonIOTracksDNNSelector",
    tracks = cms.InputTag("hltPhase2MuonIOTracks"),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    decisionThreshold = cms.double(0.990789532661438),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(44),
    dumpFeatures = cms.untracked.bool(False),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/seeds_selector_v6.onnx")
)
