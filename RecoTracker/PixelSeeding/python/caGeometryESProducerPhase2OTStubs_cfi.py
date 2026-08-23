import FWCore.ParameterSet.Config as cms

# Per-topology alpaka ESProducer (Phase2OTStubs) for the GEOMETRY-only blocks of the CA geometry SoA
# (module surface frames + per-CA-layer classification). Shared by both OT-stub CA iterations and the merger.
# Carries nothing per iteration (graph/cuts/fishbone are in each iteration's own `geometry` PSet). nLayers
# MUST equal the number of layers every CA producer sharing this geometry is configured for.
caGeometryESProducerPhase2OTStubs = cms.ESProducer('CAGeometryESProducerPhase2OTStubs@alpaka',
    # 54 = pixelTopology::Phase2OTStubs::numberOfLayers (28 pixel + 26 outer-tracker CA layers).
    nLayers = cms.uint32(54),
    appendToDataLabel = cms.string(''),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')  # Empty string = use default backend
    )
)
