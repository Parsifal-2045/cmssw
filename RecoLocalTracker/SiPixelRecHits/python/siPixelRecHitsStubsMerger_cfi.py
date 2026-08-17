import FWCore.ParameterSet.Config as cms

siPixelRecHitsStubsMerger = cms.EDProducer('SiPixelRecHitsStubsMerger@alpaka',
    pixelRecHitsSoA = cms.InputTag('siPixelRecHitsPreSplittingAlpaka'),
    stubsSoA = cms.InputTag('otStubProducer'),
    otRecHitsSoA = cms.InputTag('otRecHitsSoAConverter'),
    produceHitMask = cms.bool(True),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')  # Auto-select backend
    )
)
