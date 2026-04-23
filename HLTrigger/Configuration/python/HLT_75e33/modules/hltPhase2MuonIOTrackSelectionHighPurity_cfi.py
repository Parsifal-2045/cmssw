import FWCore.ParameterSet.Config as cms

hltPhase2MuonIOTrackSelectionHighPurity = cms.EDProducer("MuonIOTracksDNNSelector",
    tracks = cms.InputTag("hltPhase2MuonIOTracks"),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    decisionThreshold = cms.double(0.980633556842804),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(36),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/seeds_selector.onnx")
)
