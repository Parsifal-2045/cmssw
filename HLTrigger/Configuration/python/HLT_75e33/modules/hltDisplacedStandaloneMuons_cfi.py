import FWCore.ParameterSet.Config as cms

hltDisplacedStandaloneMuons = cms.EDProducer("StandAloneMuonProducer",
    InputObjects = cms.InputTag("hltDisplacedStandaloneMuonSeeds"),

    MuonTrajectoryBuilder = cms.string("StandAloneMuonTrajectoryBuilder"),

    ServiceParameters = cms.PSet(
        Propagators = cms.untracked.vstring(
            'hltESPFastSteppingHelixPropagatorAny',
            'hltESPFastSteppingHelixPropagatorOpposite'
        ),
        RPCLayers = cms.bool(True),
        UseMuonNavigation = cms.untracked.bool(True)
    ),

    STATrajBuilderParameters = cms.PSet(
        DoBackwardFilter = cms.bool(True),
        DoRefit = cms.bool(False),
        DoSeedRefit = cms.bool(False),

        FilterParameters = cms.PSet(
            NumberOfSigma = cms.double(3.0),
            FitDirection = cms.string('insideOut'),

            # Use filtered segments
            DTRecSegmentLabel = cms.InputTag("hltMaskedMuonHits"),
            CSCRecSegmentLabel = cms.InputTag("hltMaskedMuonHits"),
            RPCRecSegmentLabel = cms.InputTag(""),
            GEMRecSegmentLabel = cms.InputTag(""),
            ME0RecSegmentLabel = cms.InputTag(""),

            MaxChi2 = cms.double(1000.0),
            EnableRPCMeasurement = cms.bool(False),
            EnableDTMeasurement = cms.bool(True),
            EnableCSCMeasurement = cms.bool(True),
            EnableGEMMeasurement = cms.bool(False),
            EnableME0Measurement = cms.bool(False),

            Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),

            MuonTrajectoryUpdatorParameters = cms.PSet(
                MaxChi2 = cms.double(25.0),
                RescaleErrorFactor = cms.double(100.0),
                Granularity = cms.int32(0),
                ExcludeRPCFromFit = cms.bool(False),
                UseInvalidHits = cms.bool(True),
                RescaleError = cms.bool(False)
            )
        ),

        BWFilterParameters = cms.PSet(
            NumberOfSigma = cms.double(3.0),
            BWSeedType = cms.string('fromGenerator'),
            FitDirection = cms.string('outsideIn'),

            DTRecSegmentLabel = cms.InputTag("hltMaskedMuonHits"),
            CSCRecSegmentLabel = cms.InputTag("hltMaskedMuonHits"),
            RPCRecSegmentLabel = cms.InputTag(""),
            GEMRecSegmentLabel = cms.InputTag(""),
            ME0RecSegmentLabel = cms.InputTag(""),

            MaxChi2 = cms.double(100.0),
            EnableRPCMeasurement = cms.bool(False),
            EnableDTMeasurement = cms.bool(True),
            EnableCSCMeasurement = cms.bool(True),
            EnableGEMMeasurement = cms.bool(False),
            EnableME0Measurement = cms.bool(False),

            Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),

            MuonTrajectoryUpdatorParameters = cms.PSet(
                MaxChi2 = cms.double(25.0),
                RescaleErrorFactor = cms.double(100.0),
                Granularity = cms.int32(0),
                ExcludeRPCFromFit = cms.bool(False),
                UseInvalidHits = cms.bool(True),
                RescaleError = cms.bool(False)
            )
        ),

        NavigationType = cms.string('Standard'),
        SeedPropagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        SeedPosition = cms.string('in'),

        SeedTransformerParameters = cms.PSet(
            Fitter = cms.string('hltESPKFFittingSmootherForL2Muon'),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            NMinRecHits = cms.uint32(2),
            UseSubRecHits = cms.bool(False),
            Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            RescaleError = cms.double(100.0)
        )
    ),

    TrackLoaderParameters = cms.PSet(
        Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
        DoSmoothing = cms.bool(False),
        beamSpot = cms.InputTag("hltOnlineBeamSpot"),
        MuonUpdatorAtVertexParameters = cms.PSet(
            MaxChi2 = cms.double(1000000.0),
            BeamSpotPosition = cms.vdouble(0.0, 0.0, 0.0),
            Propagator = cms.string('hltESPFastSteppingHelixPropagatorOpposite'),
            BeamSpotPositionErrors = cms.vdouble(0.1, 0.1, 5.3)
        ),
        VertexConstraint = cms.bool(False),
        TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle')
    )
)
