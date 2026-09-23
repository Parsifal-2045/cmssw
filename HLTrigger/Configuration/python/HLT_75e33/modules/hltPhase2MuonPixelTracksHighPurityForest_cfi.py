import FWCore.ParameterSet.Config as cms

# Pixel-track (hltPhase2MuonPixelTracks, pixel-track chain) HighPurity selection with an XGBoost forest
# (MuonIOTracksForestSelector: 33 features from interface/IOTrackSelectorFeatures.h,
# compact forest in RecoMuon/L3TrackFinder/data/IO). Model and per-pT-bin
# working points come from the training records (muonHighPurityTrackSelection
# production/, deployed with deploy_cmssw.py).
hltPhase2MuonPixelTracksHighPurityForest = cms.EDProducer('MuonIOTracksForestSelector',
    tracks = cms.InputTag('hltPhase2MuonPixelTracks'),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    # IO pixel forest v3: 800 trees, 33 features, trained 2026-09-23 (git ad9361e5e2); working points = validation per-pT-bin
    # F2 set points from production/io/pixel/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/IO/muonHP_IO_pixelPath_forest_v3.bin'),
    decisionThreshold = cms.double(0.21345365047454834),
    ptBinEdges = cms.vdouble(0.0, 2.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.10667391866445541,  # pT [0, 2)
        0.21345365047454834,  # pT [2, 5)
        0.5016972422599792,  # pT [5, 10)
        0.5444226264953613,  # pT [10, 50)
        0.8626133799552917,  # pT [50, 200)
        0.21345365047454834,  # pT [200, inf) - global F2 fallback (no background in validation)
    ),
    nFeatures = cms.int32(33),
    dumpFeatures = cms.untracked.bool(False),
)
