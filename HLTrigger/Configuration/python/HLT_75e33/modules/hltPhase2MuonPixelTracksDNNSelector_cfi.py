import FWCore.ParameterSet.Config as cms
#import os
#import json
#
## Load scaler parameters from JSON file
#def load_features():
#    json_path = os.path.join(
#        os.environ.get('CMSSW_BASE', ''),
#        'src/RecoMuon/L3TrackFinder/data/features.json'
#    )
#
#    if not os.path.exists(json_path):
#        raise RuntimeError(f"Features file not found: {json_path}")
#
#    with open(json_path, 'r') as f:
#        params = json.load(f)
#
#    return params
#
#features = load_features()
#
#hltPhase2MuonPixelTracksDNNSelector = cms.EDProducer("MuonPixelTracksDNNSelector",
#    decisionThreshold = cms.double(features['threshold']),
#    usePrunedFeatures = cms.bool(features['pruned_features_indices'] is not None),
#    prunedFeatureIndices = cms.vint32(features['pruned_features_indices'] if features['pruned_features_indices'] else []),
#    nFeatures = cms.int32(features['n_features']),
#    modelPath = cms.string("RecoMuon/L3TrackFinder/data/MuonPixelTracksSelector.onnx")
#)

hltPhase2MuonPixelTracksDNNSelector = cms.EDProducer("PixelTrackDNNSelector",
    decisionThreshold = cms.double(0.9084),
    nFeatures = cms.int32(36),
    modelPath = cms.string("RecoMuon/L3TrackFinder/data/pixel_track_selector.onnx")
)
