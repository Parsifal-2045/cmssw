# Muon HLT high-purity track selectors: forest models

Compact XGBoost forests (`.bin`, format in `interface/CompactForest.h`) of the
Phase-2 muon HLT high-purity selections, one per track family and HLT chain:

- `IO/`: inside-out tracks, read by `MuonIOTracksForestSelector` (33 features,
  `interface/IOTrackSelectorFeatures.h`);
- `OI/`: outside-in tracks, read by `MuonOITracksForestSelector` (22 features,
  `interface/OITrackSelectorFeatures.h`).

File names give family, HLT chain and model version (deployed: `v3`; every
deployed retraining gets a new version, history in the training repository
README). The selectors validate every file against the feature extractor when
they load it.

| File | Selector | HLT chain | Input tracks | Plugin / cfi module | Features | Trees | Test ROC-AUC / PR-AUC | Test precision / recall / fake rejection (per-bin WPs) | Trained | Training git | md5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `IO/muonHP_IO_pixelPath_forest_v3.bin` | IO pixel | pixel-track chain (phase2MuonPixelTracksSelector, ngtScouting) | `hltPhase2MuonPixelTracks` | `MuonIOTracksForestSelector` via `hltPhase2MuonPixelTracksHighPurityForest` | 33 | 800 | 0.99977 / 0.99963 | 0.9822 / 0.9959 / 0.9893 | 2026-09-23 | `ad9361e5e2` | `51a79055ff06afb9ed3e26ff35ff11c8` |
| `IO/muonHP_IO_seedsPath_forest_v3.bin` | IO seeds | seeds chain (phase2MuonSeedsSelector) | `hltPhase2MuonIOTracks` | `MuonIOTracksForestSelector` via `hltPhase2MuonIOTrackSelectionHighPurityForest` | 33 | 1250 | 0.99952 / 0.99880 | 0.9803 / 0.9930 / 0.9957 | 2026-09-23 | `ad9361e5e2` | `8c936b2ef077e5501f18c6676b5c0cc1` |
| `OI/muonHP_OI_pixelPath_forest_v3.bin` | OI pixel | pixel-track chain (phase2MuonPixelTracksSelector, ngtScouting) | `hltPhase2L3OIMuCtfWithMaterialTracks` | `MuonOITracksForestSelector` via `_pixelOIForestSelector` | 22 | 300 | 0.99971 / 0.99722 | 0.9710 / 0.9891 / 0.9987 | 2026-09-23 | `ad9361e5e2` | `87ab5578978d4fdc4fe0815cc223140b` |
| `OI/muonHP_OI_seedsPath_forest_v3.bin` | OI general | seeds chain (phase2MuonSeedsSelector) | `hltPhase2L3OIMuCtfWithMaterialTracks` | `MuonOITracksForestSelector` via `_seedsOIForestSelector` | 22 | 200 | 0.99983 / 0.99888 | 0.9828 / 0.9939 / 0.9982 | 2026-09-23 | `ad9361e5e2` | `78fa350f5572eecb31e3498cb72f51fa` |

Training, validation and working points: muonHighPurityTrackSelection
(`production/`, forest_pipeline.py). The cfi working points are written from
the trainings' `thresholds.json` by `production/deploy_cmssw.py`, and
`production/check_consistency.py --cmssw-src` checks that files, checksums and
cfi values agree with the training records.
