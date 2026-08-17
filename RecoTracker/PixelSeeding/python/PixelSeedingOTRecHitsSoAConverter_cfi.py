import FWCore.ParameterSet.Config as cms

# OT RecHits SoA Converter -- wrapper for @alpaka backend selection
pixelSeedingOTRecHitsSoAConverter = cms.EDProducer('PixelSeedingOTRecHitsSoAConverter@alpaka',
    otRecHitSource = cms.InputTag("siPhase2RecHits"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    psOnly = cms.bool(False)
)
# Alias matching the default otRecHitsSoA InputTag of the stub producer and the merger
otRecHitsSoAConverter = pixelSeedingOTRecHitsSoAConverter.clone()
