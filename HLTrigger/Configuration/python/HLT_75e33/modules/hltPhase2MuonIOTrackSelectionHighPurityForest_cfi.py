import FWCore.ParameterSet.Config as cms

# IO-track (hltPhase2MuonIOTracks from LST seeds, seeds chain) HighPurity selection with an XGBoost forest
# (MuonIOTracksForestSelector: track features, the stub count of the matched L1 tracker muon and the raw
# matching quantities of the candidates of every muon station, interface/IOTrackSelectorFeatures.h and
# MuonChamberMatch.h; compact forest in RecoMuon/L3TrackFinder/data/IO). Model and per-pT-bin
# working points come from the training records (muonHighPurityTrackSelection
# production/, deployed with deploy_cmssw.py).
hltPhase2MuonIOTrackSelectionHighPurityForest = cms.EDProducer("MuonIOTracksForestSelector",
    tracks = cms.InputTag("hltPhase2MuonIOTracks"),
    l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
    # IO seeds forest v1: 400 trees, 46 features, trained 2026-10-10 (git e622949b3f)
    # working points = per-pT-bin recall-target set points, pinned on independent events where they populate the bin;
    # from production/io/seeds/thresholds.json (written by deploy_cmssw.py)
    modelPath = cms.FileInPath('RecoMuon/L3TrackFinder/data/IO/muonHP_IO_seedsPath_forest_v1.bin'),
    decisionThreshold = cms.double(0.591177761554718),
    ptBinEdges = cms.vdouble(0.0, 3.0, 5.0, 10.0, 50.0, 200.0),
    decisionThresholds = cms.vdouble(
        0.08881093561649323,  # pT [0, 3)
        0.4612923264503479,  # pT [3, 5)
        0.9276527166366577,  # pT [5, 10)
        0.9981892704963684,  # pT [10, 50)
        0.9987490773200989,  # pT [50, 200)
        0.999536395072937,  # pT [200, inf)
    ),
    nFeatures = cms.int32(46),
    dumpFeatures = cms.untracked.bool(False),
)
