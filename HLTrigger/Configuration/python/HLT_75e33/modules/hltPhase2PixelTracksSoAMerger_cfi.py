import FWCore.ParameterSet.Config as cms

# Merges the high-purity tracks of both iterations, runs the OT-hit attach walk, the GBL refit
# and the final duplicate removal in one pass (merge -> attach -> refit -> dedup). The tables the
# attach walk is derived from are compiled-in constants (ExtDerivedTables.h).
hltPhase2PixelTracksSoAMerger = cms.EDProducer('PixelTracksSoAMerger@alpaka',
    # Both iterations' high-purity selector outputs, still un-extended: the OT-hit extension runs here.
    inputTkSoAs = cms.VInputTag("hltPhase2PixelTrackTorchHighPuritySelector", "hltPhase2PixelTrackHighPuritySelectorDisplaced"),
    minQuality = cms.string('tight'),
    matchFraction = cms.double(0.0),

    # Cross-arm twin merge: chi2 bound of the (phi, q/pT, cotTheta) compatibility gate; <= 0 disables it.
    twinMergeNSigma2 = cms.double(8.0),

    # Attach walk: host pre-gate and per-class windows.
    # Host-quality pre-gate on chi2/ndof; threshold on the fit's covariance convention.
    extHostMaxChi2Ndof = cms.double(1.69),         # reduced chi2/ndof host gate (<=0 = off)
    extMaxSharedOwners = cms.int32(2),             # N-way sharing of attached extras (top-N claim slots)
    extPixelGateChi2Cut = cms.double(3.0),         # pixel honest-calibration chi2 cut
    extRecallReachRelax = cms.double(2.5),         # pixel reachability-envelope slack (cm); 0 = identity
    extRecallPixelFirstBudget = cms.int32(3),      # pixel-first K-seat reserve on prefer-pixel hosts
    extMaxWalkLayers = cms.int32(6),               # runtime walk visit budget K (OT layers visited per host)

    # Ceiling on attach scratch candidate capacity (effective bound = min(merged track capacity, this)).
    extRefitMaxCandidates = cms.uint32(3072),

    # Far-first disc ordering on forward OT-less hosts; sharpens sigma(pT)/pT at 2.8<|eta|<3.5.
    extAttachFarMaxWin = cms.int32(1),

    # Attach pre-gate |eta| ceiling (top edge of the forward-eta pocket band).
    extMaxAbsEta = cms.double(4.5),

    # Derived selection knob: the probability mass of the correct-hit pull distribution that the
    # attach window has to contain. It ties the window, the gate and the ranking together.
    extDerivedEps = cms.double(0.55),
)
