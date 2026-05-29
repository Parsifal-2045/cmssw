import FWCore.ParameterSet.Config as cms

hltPhase2MuonPixelTracksHighPurity = cms.EDProducer('MuonIOTracksDNNSelector',
    tracks = cms.InputTag('hltPhase2MuonPixelTracks'),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    decisionThreshold = cms.double(0.9868386387825012),
    useL1TkMuFeatures = cms.bool(True),
    useStubFeatures = cms.bool(True),
    nFeatures = cms.int32(44),
    dumpFeatures = cms.untracked.bool(False),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/tuned_pixel_track_selector_v6.onnx")
)
