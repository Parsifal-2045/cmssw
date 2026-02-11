import FWCore.ParameterSet.Config as cms

hltPhase2MuonPixelTracksDNNSelector = cms.EDProducer("MuonIOTracksDNNSelector",
    tracks = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracks"),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    decisionThreshold = cms.double(0.895016074180603),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(36),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/pixel_track_selector.onnx")
)
