import FWCore.ParameterSet.Config as cms

hltPhase2MuonPixelTracksHighPurity = cms.EDProducer('MuonIOTracksDNNSelector',
    tracks = cms.InputTag('hltPhase2MuonPixelTracks'),
    l1TkMuons = cms.InputTag('l1tTkMuonsGmt'),
    modelPath = cms.string('RecoMuon/L3TrackFinder/data/pixel_track_selector.onnx'),
    decisionThreshold = cms.double(0.895016074180603),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(36)
)
