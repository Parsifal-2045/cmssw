import FWCore.ParameterSet.Config as cms

# EventSetup producer of the CA-ordered module geometry SoA for the Phase2OTStubs chain.
# Shared by both OT-stub CA iterations and the merger; carries nothing iteration-specific.
from RecoTracker.PixelSeeding.caGeometryESProducerPhase2OTStubs_cfi import (
    caGeometryESProducerPhase2OTStubs as _caGeometryESProducerPhase2OTStubs,
)

hltCAGeometryESProducerPhase2OTStubs = _caGeometryESProducerPhase2OTStubs.clone(
    # 54 CA layers (28 pixel + 26 OT): must equal the size of the `layers` table of every CA producer
    # sharing this geometry (hltPhase2PixelTracksSoAWithStubs_cfi and the displaced iteration).
    nLayers = 54,
)
