import FWCore.ParameterSet.Config as cms

# Per-topology alpaka ESProducer (Phase2OT) for the geometry-only blocks of the CA geometry SoA
# (module surface frames + per-CA-layer classification). Not used by the HLT menu, whose non-stub
# CA producers build their own geometry. nLayers MUST equal the number of layers every CA producer
# sharing this geometry is configured for (they throw on a mismatch).
caGeometryESProducerPhase2OT = cms.ESProducer('CAGeometryESProducerPhase2OT@alpaka',
    # 31 = pixelTopology::Phase2OT::numberOfLayers (28 pixel + 3 considered OT barrel layers,
    # phase2PixelTopology::nLayersPix + nLayersOT). NOT 54: that is the Phase2OTStubs layer table.
    nLayers = cms.uint32(31),
    appendToDataLabel = cms.string(''),
    alpaka = cms.untracked.PSet(
        backend = cms.untracked.string('')  # Empty string = use default backend
    )
)
