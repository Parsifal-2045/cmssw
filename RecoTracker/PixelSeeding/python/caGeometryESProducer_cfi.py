import FWCore.ParameterSet.Config as cms

# Generic entry point for the CA geometry ESProducer; aliases the Phase2OTStubs topology.
# nLayers MUST match the layer count every producer sharing this geometry is configured for.
from RecoTracker.PixelSeeding.caGeometryESProducerPhase2OTStubs_cfi import (
    caGeometryESProducerPhase2OTStubs as caGeometryESProducer,
)
