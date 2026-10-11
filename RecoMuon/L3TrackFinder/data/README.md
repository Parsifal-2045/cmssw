# Muon HLT high-purity track selectors: forest models

Compact XGBoost forests (`.bin`, format in `interface/CompactForest.h`) of the
Phase-2 muon HLT high-purity selections, one per track family and HLT chain:

- `IO/`: inside-out tracks, read by `MuonIOTracksForestSelector` (46 features: track
  quantities and the raw matching quantities of the muon-station candidates, `interface/IOTrackSelectorFeatures.h`
  and `interface/MuonChamberMatch.h`);
- `OI/`: outside-in tracks, read by `MuonOITracksForestSelector` (22 features,
  `interface/OITrackSelectorFeatures.h`).

File names give family, HLT chain and model version (deployed: `v1`; every
deployed retraining gets a new version, history in the training repository
README). The selectors validate every file against the feature extractor when
they load it.

| File | Selector | HLT chain | Input tracks | Plugin / cfi module | Features | Trees | Test ROC-AUC / PR-AUC | Test precision / recall / fake rejection (per-bin WPs) | Trained | Training git | md5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `IO/muonHP_IO_pixelPath_forest_v1.bin` | IO pixel | pixel-track chain (phase2MuonPixelTracksSelector, ngtScouting) | `hltPhase2MuonPixelTracks` | `MuonIOTracksForestSelector` via `hltPhase2MuonPixelTracksHighPurityForest` | 46 | 1100 | 0.99498 / 0.99020 | 0.8781 / 0.9901 / 0.9287 | 2026-10-10 | `e622949b3f` | `042a283f9006688ab25a7fbe82ccc77b` |
| `IO/muonHP_IO_seedsPath_forest_v1.bin` | IO seeds | seeds chain (phase2MuonSeedsSelector) | `hltPhase2MuonIOTracks` | `MuonIOTracksForestSelector` via `hltPhase2MuonIOTrackSelectionHighPurityForest` | 46 | 400 | 0.99764 / 0.99211 | 0.8621 / 0.9889 / 0.9585 | 2026-10-10 | `e622949b3f` | `b418758c4ebc15e53c5827e7c7901c59` |
| `OI/muonHP_OI_pixelPath_forest_v1.bin` | OI pixel | pixel-track chain (phase2MuonPixelTracksSelector, ngtScouting) | `hltPhase2L3OIMuCtfWithMaterialTracks` | `MuonOITracksForestSelector` via `_pixelOIForestSelector` | 22 | 200 | 0.99151 / 0.86929 | 0.3433 / 0.9856 / 0.9200 | 2026-10-10 | `e622949b3f` | `e7fc5a9252172034f9df2707d108e24d` |
| `OI/muonHP_OI_seedsPath_forest_v1.bin` | OI general | seeds chain (phase2MuonSeedsSelector) | `hltPhase2L3OIMuCtfWithMaterialTracks` | `MuonOITracksForestSelector` via `_seedsOIForestSelector` | 22 | 200 | 0.99254 / 0.88157 | 0.3591 / 0.9840 / 0.9280 | 2026-10-10 | `e622949b3f` | `9bfee380b45bf307777ab0361fbd318a` |

Training, validation and working points: muonHighPurityTrackSelection
(`production/`, forest_pipeline.py). The cfi working points are written from
the trainings' `thresholds.json` by `production/deploy_cmssw.py`, and
`production/check_consistency.py --cmssw-src` checks that files, checksums and
cfi values agree with the training records.
