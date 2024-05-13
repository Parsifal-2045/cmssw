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
      minMomentum_{pset.getParameter<double>("MinPL1Tk")},
      maxMomentum_{pset.getParameter<double>("MaxPL1Tk")},
      dRCone_{pset.getParameter<double>("StubMatchDR")},
      maxEtaBarrel_{pset.getParameter<double>("maximumEtaBarrel")},
      maxEtaOverlap_{pset.getParameter<double>("maximumEtaOverlap")}{}

void Phase2L2MuonSeedCreator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("InputObjects", edm::InputTag("hltGmtStage2Digis"));
  desc.add<edm::InputTag>("CSCRecSegmentLabel", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("DTRecSegmentLabel", edm::InputTag("hltDt4DSegments"));
  desc.add<double>("MinPL1Tk", 3.5);
  desc.add<double>("MaxPL1Tk", 200);
  desc.add<double>("StubMatchDR", 0.25);
  desc.add<double>("maximumEtaBarrel", 0.7);
  desc.add<double>("maximumEtaOverlap", 1.3);
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

  LogDebug(metname) << "Number of muons " << l1TkMuColl->size() << std::endl;

  //l1t::TrackerMuonRef::const_iterator it;
  //l1t::TrackerMuonRef::key_type
  size_t l1TkMuIndex = 0;

  for (; l1TkMuIndex != l1TkMuColl->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuColl, l1TkMuIndex);
    Type muonType;
    if (std::abs(l1TkMuRef->phEta()) < maxEtaBarrel_) {
      muonType = barrel;
    } else if (std::abs(l1TkMuRef->phEta()) > maxEtaOverlap_) {
      muonType = endcap;
    } else {
      muonType = overlap;
    }
    output->emplace_back(createSeed(muonType, *cscSegments, *dtSegments, l1TkMuRef));
  }

  iEvent.put(std::move(output));
}

L2MuonTrajectorySeed Phase2L2MuonSeedCreator::createSeed(Type muonType,  // barrel (0), overlap (1), endcap (2)
                                                         CSCSegmentCollection cscSegments,
                                                         DTRecSegment4DCollection dtSegments,
                                                         l1t::TrackerMuonRef l1TkMuRef) {
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
    double bestDr2 = dRCone_ * dRCone_;
    int bestSegIndex = -1;
    bool bestInCsc = false;

    if (muonType == barrel or muonType == overlap) {
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
          bestSegIndex = s;
        }
      }
    }

    if (muonType == endcap or muonType == overlap) {
      for (size_t s = 0; s != cscSegments.size(); ++s) {
        int endcap = (l1TkMuRef->phEta() > 0) ? 1 : 2;  // endcap for CSC DetId (1 -> Forward, 2 -> Backward)
        CSCDetId stubId = CSCDetId(endcap, stub->depthRegion(), stub->etaRegion(), stub->phiRegion());
        CSCDetId segId = cscSegments[s].cscDetId();
        if (!(stubId == segId)) {
          continue;
        }

        GlobalPoint segPos = cscGeometry_->idToDet(segId)->toGlobal(cscSegments[s].localPosition());
        // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
        double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
        if (dr2 < bestDr2) {
          bestDr2 = dr2;
          bestSegIndex = s;
          bestInCsc = true;
        }
      }
    }

    // auto bestSegment = (bestInCsc) ? cscSegments[bestSegIndex] : dtSegments[bestSegIndex];  // CSCSegment vs DTRecSegment4D
    // auto segId = (bestInCsc) ? bestSegment->cscDetId() : bestSegment->chamberId();
    // const GeomDet* geomDet = (bestInCsc) ? cscGeometry_->idToDet(*segId) : dtGeometry_->idToDet(*segId);

    if (!bestInCsc) {
      auto bestSegment = dtSegments[bestSegIndex];
      auto segId = bestSegment.chamberId();
      const GeomDet* geomDet = dtGeometry_->idToDet(segId);

      for (auto const& recHit : bestSegment.recHits()) {
        container.push_back(recHit);
      }

      // Fill the local trajectory parameters
      LocalPoint segPos = bestSegment.localPosition();
      LocalVector segDirFromPos = bestSegment.localDirection();
      LocalTrajectoryParameters param(segPos, segDirFromPos, charge);

      // Get the global direction
      GlobalVector mom = geomDet->toGlobal(bestSegment.localPosition()) - GlobalPoint();

      GlobalVector polar(
          GlobalVector::Spherical(mom.theta(), geomDet->toGlobal(bestSegment.localDirection()).phi(), 1.));
      // Scale magnitude of total momentum
      polar *= fabs(l1Pt) / polar.perp();
      // Get error matrix
      mat = bestSegment.parametersError().similarityT(bestSegment.projectionMatrix());
      p_err = (sptmean * sptmean) / (polar.mag() * polar.mag() * l1Pt * l1Pt);
      mat[0][0] = p_err;
      LocalTrajectoryError error(asSMatrix<5>(mat));

      // Create the TrajectoryStateOnSurface
      TrajectoryStateOnSurface tsos(param, error, geomDet->surface(), &*magneticField_);

      // Transform it in a TrajectoryStateOnSurface
      PTrajectoryStateOnDet seedTSOS = trajectoryStateTransform::persistentState(tsos, segId.rawId());
      theSeed = L2MuonTrajectorySeed(seedTSOS, container, alongMomentum, l1TkMuRef);
    } else {
      auto bestSegment = cscSegments[bestSegIndex];
      auto segId = bestSegment.cscDetId();
      const GeomDet* geomDet = dtGeometry_->idToDet(segId);

      for (auto const& recHit : bestSegment.recHits()) {
        container.push_back(recHit);
      }

      // Fill the local trajectory parameters
      LocalPoint segPos = bestSegment.localPosition();
      LocalVector segDirFromPos = bestSegment.localDirection();
      LocalTrajectoryParameters param(segPos, segDirFromPos, charge);

      // Get the global direction
      GlobalVector mom = geomDet->toGlobal(bestSegment.localPosition()) - GlobalPoint();

      GlobalVector polar(
          GlobalVector::Spherical(mom.theta(), geomDet->toGlobal(bestSegment.localDirection()).phi(), 1.));
      // Scale magnitude of total momentum
      polar *= fabs(l1Pt) / polar.perp();
      // Get error matrix
      mat = bestSegment.parametersError().similarityT(bestSegment.projectionMatrix());
      p_err = (sptmean * sptmean) / (polar.mag() * polar.mag() * l1Pt * l1Pt);
      mat[0][0] = p_err;
      LocalTrajectoryError error(asSMatrix<5>(mat));

      // Create the TrajectoryStateOnSurface
      TrajectoryStateOnSurface tsos(param, error, geomDet->surface(), &*magneticField_);

      // Transform it in a TrajectoryStateOnSurface
      PTrajectoryStateOnDet seedTSOS = trajectoryStateTransform::persistentState(tsos, segId.rawId());
      theSeed = L2MuonTrajectorySeed(seedTSOS, container, alongMomentum, l1TkMuRef);
    }
  }
  return theSeed;
}

DEFINE_FWK_MODULE(Phase2L2MuonSeedCreator);
