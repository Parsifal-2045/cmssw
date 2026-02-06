import FWCore.ParameterSet.Config as cms

hltPhase2MuonPixelTracksDNNSelector = cms.EDProducer("MuonIOTracksDNNSelector",
    decisionThreshold = cms.double(0.895016074180603),
    nFeatures = cms.int32(36),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/pixel_track_selector.onnx")
)
