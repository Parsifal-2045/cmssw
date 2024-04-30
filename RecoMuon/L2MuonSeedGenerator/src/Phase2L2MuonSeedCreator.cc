/**
 *  See header file for a description of this class.
 *
 */

#include "RecoMuon/L2MuonSeedGenerator/src/Phase2L2MuonSeedCreator.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"

#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/L1TMuonPhase2/interface/MuonStub.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "gsl/gsl_statistics.h"

/*
 * Constructor
 */
Phase2L2MuonSeedCreator::Phase2L2MuonSeedCreator(const edm::ParameterSet& pset)
    : l1TkMuCollToken_{consumes(pset.getParameter<edm::InputTag>("InputObjects"))},
      cscSegmentCollToken_{consumes(pset.getParameter<edm::InputTag>("CSCRecSegmentLabel"))},
      dtSegmentCollToken_{consumes(pset.getParameter<edm::InputTag>("DTRecSegmentLabel"))},
      cscGeometryToken_{esConsumes<CSCGeometry, MuonGeometryRecord>()},
      dtGeometryToken_{esConsumes<DTGeometry, MuonGeometryRecord>()},
      magneticFieldToken_{esConsumes<MagneticField, IdealMagneticFieldRecord>()},
      muonLayersToken_{esConsumes<MuonDetLayerGeometry, MuonRecoGeometryRecord>()},
      minMomentum_{pset.getParameter<double>("minimumSeedPt")},
      maxMomentum_{pset.getParameter<double>("maximumSeedPt")},
      defaultMomentum_{pset.getParameter<double>("defaultSeedPt")},
      debug_{pset.getParameter<bool>("DebugMuonSeed")},
      sysError_{pset.getParameter<double>("SeedPtSystematics")},
      DT12{pset.getParameter<std::vector<double>>("DT_12")},
      DT13{pset.getParameter<std::vector<double>>("DT_13")},
      DT14{pset.getParameter<std::vector<double>>("DT_14")},
      DT23{pset.getParameter<std::vector<double>>("DT_23")},
      DT24{pset.getParameter<std::vector<double>>("DT_24")},
      DT34{pset.getParameter<std::vector<double>>("DT_34")},
      DT12_1{pset.getParameter<std::vector<double>>("DT_12_1_scale")},
      DT12_2{pset.getParameter<std::vector<double>>("DT_12_2_scale")},
      DT13_1{pset.getParameter<std::vector<double>>("DT_13_1_scale")},
      DT13_2{pset.getParameter<std::vector<double>>("DT_13_2_scale")},
      DT14_1{pset.getParameter<std::vector<double>>("DT_14_1_scale")},
      DT14_2{pset.getParameter<std::vector<double>>("DT_14_2_scale")},
      DT23_1{pset.getParameter<std::vector<double>>("DT_23_1_scale")},
      DT23_2{pset.getParameter<std::vector<double>>("DT_23_2_scale")},
      DT24_1{pset.getParameter<std::vector<double>>("DT_24_1_scale")},
      DT24_2{pset.getParameter<std::vector<double>>("DT_24_2_scale")},
      DT34_1{pset.getParameter<std::vector<double>>("DT_34_1_scale")},
      DT34_2{pset.getParameter<std::vector<double>>("DT_34_2_scale")} {}

/*
 * Destructor
 */
// Phase2L2MuonSeedCreator::~Phase2L2MuonSeedCreator() {}

void Phase2L2MuonSeedCreator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("InputObjects", edm::InputTag("hltGmtStage2Digis"));
  desc.add<double>("MinPL1Tk", 3.5);
  desc.add<double>("MaxPL1Tk", 200);
  descriptions.add("Phase2L2MuonSeedCreator", desc);
}

void Phase2L2MuonSeedCreator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "Muon|RecoMuon|Phase2L2MuonSeedCreator";
  MuonPatternRecoDumper debug;

  auto output = std::make_unique<L2MuonTrajectorySeedCollection>();

  auto const l1TkMuColl = iEvent.getHandle(l1TkMuCollToken_);

  auto cscSegments = iEvent.getHandle(cscSegmentCollToken_);
  auto dtSegments = iEvent.getHandle(dtSegmentCollToken_);

  cscGeometry_ = iSetup.getHandle(cscGeometryToken_);
  dtGeometry_ = iSetup.getHandle(dtGeometryToken_);
  magneticField_ = iSetup.getHandle(magneticFieldToken_);

  double dRCone = 0.01;

  LogDebug(metname) << "Number of muons " << l1TkMuColl->size() << std::endl;

  //l1t::TrackerMuonRef::const_iterator it;
  l1t::TrackerMuonRef::key_type l1TkMuIndex = 0;

  for (; l1TkMuIndex != l1TkMuColl->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuColl, l1TkMuIndex);
    output->emplace_back(createSeed(*cscSegments, *dtSegments, l1TkMuRef, dRCone));
  }

  iEvent.put(std::move(output));
}

/*
 * createSeed
 *
 * Note type = 1 --> CSC
 *           = 2 --> Overlap
 *           = 3 --> DT
 */

L2MuonTrajectorySeed Phase2L2MuonSeedCreator::createSeed(CSCSegmentCollection cscSegments,
                                                         DTRecSegment4DCollection dtSegments,
                                                         l1t::TrackerMuonRef l1TkMuRef,
                                                         const double& dRCone) {
  double l1Pt = l1TkMuRef->phPt();
  double sptmean = minMomentum_;
  int charge = l1TkMuRef->phCharge();

  if (l1Pt < minMomentum_) {
    l1Pt = minMomentum_;
    sptmean = minMomentum_;
  } else if (l1Pt > maxMomentum_) {
    l1Pt = maxMomentum_;
    sptmean = maxMomentum_ * 0.25;
  }

  l1t::MuonStubRefVector stubRefs = l1TkMuRef->stubs();

  AlgebraicSymMatrix mat(5, 0);
  double p_err = 0;

  edm::OwnVector<TrackingRecHit> container;
  L2MuonTrajectorySeed theSeed;

  for (auto stub : stubRefs) {
    double bestDr2 = dRCone * dRCone;
    int bestSeg = -1;
    for (size_t s = 0; s != dtSegments.size(); ++s) {
      DTChamberId stubId = DTChamberId(stub->etaRegion(), stub->depthRegion(), stub->phiRegion());
      DetId segId = dtSegments[s].chamberId();
      if (!(stubId == segId)) {
        continue;
      }

      GlobalPoint segPos = dtGeometry_->idToDet(segId)->toGlobal(dtSegments[s].localPosition());
      // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
      double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
      if (dr2 < bestDr2) {
        bestDr2 = dr2;
        bestSeg = s;
      }
    }
    for (auto const& recHit : dtSegments[bestSeg].recHits()) {
      container.push_back(recHit);
    }
    DTChamberId segId = dtSegments[bestSeg].chamberId();
    const GeomDet* geomDet = dtGeometry_->idToDet(segId);

    // Fill the local trajectory parameters
    LocalPoint segPos = dtSegments[bestSeg].localPosition();
    LocalVector segDirFromPos = dtSegments[bestSeg].localDirection();
    LocalTrajectoryParameters param(segPos, segDirFromPos, charge);

    // Get the global direction
    GlobalVector mom = geomDet->toGlobal(dtSegments[bestSeg].localPosition()) - GlobalPoint();

    GlobalVector polar(
        GlobalVector::Spherical(mom.theta(), geomDet->toGlobal(dtSegments[bestSeg].localDirection()).phi(), 1.));
    // Scale magnitude of total momentum
    polar *= fabs(l1Pt) / polar.perp();
    // Get error matrix
    mat = dtSegments[bestSeg].parametersError().similarityT(dtSegments[bestSeg].projectionMatrix());
    p_err = (sptmean * sptmean) / (polar.mag() * polar.mag() * l1Pt * l1Pt);
    mat[0][0] = p_err;
    LocalTrajectoryError error(asSMatrix<5>(mat));

    // Create the TrajectoryStateOnSurface
    TrajectoryStateOnSurface tsos(param, error, geomDet->surface(), &*magneticField_);

    // Transform it in a TrajectoryStateOnSurface
    PTrajectoryStateOnDet seedTSOS = trajectoryStateTransform::persistentState(tsos, segId.rawId());
    theSeed = L2MuonTrajectorySeed(seedTSOS, container, alongMomentum, l1TkMuRef);
  }
  return theSeed;
}

DEFINE_FWK_MODULE(Phase2L2MuonSeedCreator);
