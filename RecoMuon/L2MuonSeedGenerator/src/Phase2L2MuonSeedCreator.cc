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
      maxEtaOverlap_{pset.getParameter<double>("maximumEtaOverlap")} {
  produces<L2MuonTrajectorySeedCollection>();
}

void Phase2L2MuonSeedCreator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("InputObjects", edm::InputTag("l1tTkMuonsGmt"));
  desc.add<edm::InputTag>("CSCRecSegmentLabel", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("DTRecSegmentLabel", edm::InputTag("hltDt4DSegments"));
  desc.add<double>("MinPL1Tk", 3.5);
  desc.add<double>("MaxPL1Tk", 200);
  desc.add<double>("StubMatchDR", 0.1);
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

  std::cout << "Number of muons " << l1TkMuColl->size() << std::endl;

  //l1t::TrackerMuonRef::const_iterator it;
  //l1t::TrackerMuonRef::key_type
  size_t l1TkMuIndex = 0;

  for (; l1TkMuIndex != l1TkMuColl->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuColl, l1TkMuIndex);
    std::cout << "phEta: " << l1TkMuRef->phEta() << ", " << "phPhi: " << l1TkMuRef->phPhi() << std::endl;
    Type muonType = overlap;
    if (std::abs(l1TkMuRef->phEta()) < maxEtaBarrel_) {
      muonType = barrel;
    } else if (std::abs(l1TkMuRef->phEta()) > maxEtaOverlap_) {
      muonType = endcap;
    }
    output->emplace_back(createSeed(muonType, *cscSegments, *dtSegments, l1TkMuRef));
  }

  iEvent.put(std::move(output));
}

L2MuonTrajectorySeed Phase2L2MuonSeedCreator::createSeed(const Type muonType,  // barrel (0), overlap (1), endcap (2)
                                                         CSCSegmentCollection cscSegments,
                                                         DTRecSegment4DCollection dtSegments,
                                                         l1t::TrackerMuonRef l1TkMuRef) {
  double l1Pt = l1TkMuRef->phPt();
  //double sptmean = minMomentum_;
  //int charge = l1TkMuRef->phCharge();
  //
  //if (l1Pt < minMomentum_) {
  //  l1Pt = minMomentum_;
  //  sptmean = minMomentum_;
  //} else if (l1Pt > maxMomentum_) {
  //  l1Pt = maxMomentum_;
  //  sptmean = maxMomentum_ * 0.25;
  //}

  l1t::MuonStubRefVector stubRefs = l1TkMuRef->stubs();

  std::cout << "Number of stubs per L1TkMu: " << stubRefs.size() << '\n';
  std::cout << "Number of DT segments in event: " << dtSegments.size() << '\n';
  std::cout << "Number of CSC segments in event: " << cscSegments.size() << '\n';

  //AlgebraicSymMatrix mat(5, 0);
  //double p_err = 0;

  edm::OwnVector<TrackingRecHit> container;
  L2MuonTrajectorySeed theSeed;

  for (auto stub : stubRefs) {
    stub->print();
    double bestDr2 = dRCone_ * dRCone_;
    int bestSegIndex = -1;
    bool bestInCsc = false;

    switch (muonType) {
      case barrel:
        if (!stub->isBarrel()) {
          continue;
        }

        std::cout << "BARREL checking " << dtSegments.size() << " DT segments" << '\n';
        for (size_t s = 0; s != dtSegments.size(); ++s) {
          // wheel->eta station->depth sector->phi
          DTChamberId stubId = DTChamberId(stub->etaRegion(), stub->depthRegion(), stub->phiRegion() + 1); // converting online sectors to offline
          DTChamberId segId = dtSegments[s].chamberId();
          std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
                    << "segId: " << segId << ", RAW: " << segId.rawId() << '\n';

          if (!(stubId == segId)) {
            continue;
          }
          std::cout << "BARREL found match dtdetID" << '\n';

          GlobalPoint segPos = dtGeometry_->idToDet(segId)->toGlobal(dtSegments[s].localPosition());
          // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
          double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
          if (dr2 <= bestDr2) {
            std::cout << "BARREL found match dR" << '\n';
            bestDr2 = dr2;
            bestSegIndex = s;
          }
        }
        std::cout << "BARREL best segment: " << bestSegIndex << " found in csc? " << bestInCsc << '\n';
        break;

      case endcap:
        if (!stub->isEndcap()) {
          continue;
        }

        std::cout << "ENDCAP checking " << cscSegments.size() << " CSC segments" << '\n';
        for (size_t s = 0; s != cscSegments.size(); ++s) {
          int endcap = (l1TkMuRef->phEta() > 0) ? 1 : 2;  // endcap for CSC DetId (1 -> Forward, 2 -> Backward)
          CSCDetId stubId = CSCDetId(endcap, stub->depthRegion(), 6 - stub->etaRegion(), stub->phiRegion()); // 6 - etaRegion mapping stub ring with csc ring
          CSCDetId segId = cscSegments[s].cscDetId();
          std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
                    << "segId: " << segId << ", RAW: " << segId.rawId() << '\n';

          if (!(stubId == segId)) {
            continue;
          }
          std::cout << "ENDCAP found match cscdetID" << '\n';

          GlobalPoint segPos = cscGeometry_->idToDet(segId)->toGlobal(cscSegments[s].localPosition());
          // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
          double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
          if (dr2 <= bestDr2) {
            std::cout << "ENDCAP found match dR" << '\n';
            bestDr2 = dr2;
            bestSegIndex = s;
            bestInCsc = true;
          }
        }
        std::cout << "ENDCAP best segment: " << bestSegIndex << " found in csc? " << bestInCsc << '\n';
        break;

      case overlap:
        std::cout << "OVERLAP checking " << dtSegments.size() << " DT segments" << '\n';
        for (size_t s = 0; s != dtSegments.size(); ++s) {
          if (!stub->isBarrel()) {
            continue;
          }
          DTChamberId stubId = DTChamberId(stub->etaRegion(), stub->depthRegion(), stub->phiRegion() + 1);
          DTChamberId segId = dtSegments[s].chamberId();
          std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
                    << "segId: " << segId << ", RAW: " << segId.rawId() << '\n';
          if (!(stubId == segId)) {
            continue;
          }
          std::cout << "OVERLAP found match dtDetID" << '\n';

          GlobalPoint segPos = dtGeometry_->idToDet(segId)->toGlobal(dtSegments[s].localPosition());
          // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
          double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
          if (dr2 <= bestDr2) {
            std::cout << "OVERLAP found match dR in DT" << '\n';
            bestDr2 = dr2;
            bestSegIndex = s;
          }
        }

        std::cout << "OVERLAP checking also " << cscSegments.size() << " CSC segments" << '\n';
        for (size_t s = 0; s != cscSegments.size(); ++s) {
          if (!stub->isEndcap()) {
            continue;
          }
          int endcap = (l1TkMuRef->phEta() > 0) ? 1 : 2;  // endcap for CSC DetId (1 -> Forward, 2 -> Backward)
          CSCDetId stubId = CSCDetId(endcap, stub->depthRegion(), 6 - stub->etaRegion(), stub->phiRegion());
          CSCDetId segId = cscSegments[s].cscDetId();
          std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
                    << "segId: " << segId << ", RAW: " << segId.rawId() << '\n';

          if (!(stubId == segId)) {
            continue;
          }
          std::cout << "OVERLAP found match cscDetID" << '\n';

          GlobalPoint segPos = cscGeometry_->idToDet(segId)->toGlobal(cscSegments[s].localPosition());
          // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
          double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
          if (dr2 <= bestDr2) {
            std::cout << "OVERLAP found match dR in CSC" << '\n';
            bestDr2 = dr2;
            bestSegIndex = s;
            bestInCsc = true;
          }
        }
        std::cout << "OVERLAP best segment: " << bestSegIndex << " found in csc? " << bestInCsc << '\n';
        break;

      default:
        std::cout << "L1TkMu must be either barrel, endcap or overlap" << std::endl;
        break;
    }

    // auto bestSegment = (bestInCsc) ? cscSegments[bestSegIndex] : dtSegments[bestSegIndex];  // CSCSegment vs DTRecSegment4D
    // auto segId = (bestInCsc) ? bestSegment->cscDetId() : bestSegment->chamberId();
    // const GeomDet* geomDet = (bestInCsc) ? cscGeometry_->idToDet(*segId) : dtGeometry_->idToDet(*segId);
    /*
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
      std::cout << "best in CSC " << bestSegIndex << '\n';
      auto bestSegment = cscSegments[bestSegIndex];
      std::cout << "best segment initialised" << '\n';
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
    */
    theSeed = L2MuonTrajectorySeed();
  }
  return theSeed;
}

DEFINE_FWK_MODULE(Phase2L2MuonSeedCreator);
