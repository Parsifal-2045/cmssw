import FWCore.ParameterSet.Config as cms

# Pixel-track (hltPhase2MuonPixelTracks, pixel-track chain) HighPurity selection with an XGBoost forest
# (MuonIOTracksForestSelector: track features, the stub count of the matched L1 tracker muon and the raw
# matching quantities of the candidates of every muon station, interface/IOTrackSelectorFeatures.h and
# MuonChamberMatch.h; compact forest in RecoMuon/L3TrackFinder/data/IO). Model and per-pT-bin
# working points come from the training records (muonHighPurityTrackSelection
# production/, deployed with deploy_cmssw.py).
hltPhase2MuonPixelTracksHighPurityForest = cms.EDProducer('MuonIOTracksForestSelector',
    tracks = cms.InputTag('hltPhase2MuonPixelTracks'),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    # IO pixel forest v1: 1100 trees, 46 features, trained 2026-10-10 (git e622949b3f)
    # working points = per-pT-bin recall-target set points, pinned on independent events where they populate the bin;
    # from production/io/pixel/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/IO/muonHP_IO_pixelPath_forest_v1.bin'),
    decisionThreshold = cms.double(0.7648054361343384),
    ptBinEdges = cms.vdouble(0.0, 3.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.26123544573783875,  # pT [0, 3)
        0.5858924984931946,  # pT [3, 5)
        0.9377449750900269,  # pT [5, 10)
        0.9989855289459229,  # pT [10, 50)
        0.9998182058334351,  # pT [50, 200)
        0.9997901320457458,  # pT [200, inf)
    ),
    nFeatures = cms.int32(46),
    dumpFeatures = cms.untracked.bool(False),
)
