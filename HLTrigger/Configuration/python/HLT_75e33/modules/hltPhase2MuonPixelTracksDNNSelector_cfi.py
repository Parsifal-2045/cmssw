import FWCore.ParameterSet.Config as cms
import json
import os

# Load scaler parameters from JSON file
def load_scaler_params():
    json_path = os.path.join(
        os.environ.get('CMSSW_BASE', ''),
        'src/RecoMuon/L3TrackFinder/data/model_params.json'
    )

    if not os.path.exists(json_path):
        raise RuntimeError(f"Scaler params file not found: {json_path}")

    with open(json_path, 'r') as f:
        params = json.load(f)

    return params

scaler_params = load_scaler_params()

hltPhase2MuonPixelTracksDNNSelector = cms.EDProducer("MuonPixelTracksDNNSelector",
    # Scaler parameters from JSON
    scalerMean = cms.vdouble(scaler_params['mean']),
    scalerScale = cms.vdouble(scaler_params['scale']),
    # Model configuration
    decisionThreshold = cms.double(scaler_params['threshold']),
    usePrunedFeatures = cms.bool(scaler_params['pruned_features_indices'] is not None),
    prunedFeatureIndices = cms.vint32(scaler_params['pruned_features_indices'] if scaler_params['pruned_features_indices'] else []),
    nFeatures = cms.int32(scaler_params['n_features'])
)
