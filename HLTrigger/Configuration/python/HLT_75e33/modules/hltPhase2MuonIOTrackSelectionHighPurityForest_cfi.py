import FWCore.ParameterSet.Config as cms

# IO-track (hltPhase2MuonIOTracks from LST seeds, seeds chain) HighPurity selection with an XGBoost forest
# (MuonIOTracksForestSelector: 33 features from interface/IOTrackSelectorFeatures.h,
# compact forest in RecoMuon/L3TrackFinder/data/IO). Model and per-pT-bin
# working points come from the training records (muonHighPurityTrackSelection
# production/, deployed with deploy_cmssw.py).
hltPhase2MuonIOTrackSelectionHighPurityForest = cms.EDProducer("MuonIOTracksForestSelector",
    tracks = cms.InputTag("hltPhase2MuonIOTracks"),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    # IO seeds forest v3: 1250 trees, 33 features, trained 2026-09-23 (git ad9361e5e2); working points = validation per-pT-bin
    # F2 set points from production/io/seeds/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/IO/muonHP_IO_seedsPath_forest_v3.bin'),
    decisionThreshold = cms.double(0.17353090643882751),
    ptBinEdges = cms.vdouble(0.0, 2.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.03589987754821777,  # pT [0, 2)
        0.16624124348163605,  # pT [2, 5)
        0.34693241119384766,  # pT [5, 10)
        0.47059762477874756,  # pT [10, 50)
        0.5600723624229431,  # pT [50, 200)
        0.8556238412857056,  # pT [200, inf)
    ),
    nFeatures = cms.int32(33),
    dumpFeatures = cms.untracked.bool(False),
)
