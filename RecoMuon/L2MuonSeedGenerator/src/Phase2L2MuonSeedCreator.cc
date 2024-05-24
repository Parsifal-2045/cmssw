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
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/L1TMuonPhase2/interface/MuonStub.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Vector/ThreeVector.h"

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
      maxEtaOverlap_{pset.getParameter<double>("maximumEtaOverlap")},
      propagatorName_{pset.getParameter<string>("Propagator")} {
  // Service parameters
  edm::ParameterSet serviceParameters = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  // Services
  service_ = std::make_unique<MuonServiceProxy>(serviceParameters, consumesCollector());
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
  desc.add<string>("Propagator", "SteppingHelixPropagatorAny");

  // Service parameters
  edm::ParameterSetDescription psd0;
  psd0.addUntracked<std::vector<std::string>>("Propagators", {"SteppingHelixPropagatorAny"});
  psd0.add<bool>("RPCLayers", true);
  psd0.addUntracked<bool>("UseMuonNavigation", true);
  desc.add<edm::ParameterSetDescription>("ServiceParameters", psd0);
  descriptions.add("Phase2L2MuonSeedCreator", desc);
}

void Phase2L2MuonSeedCreator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "Muon|RecoMuon|Phase2L2MuonSeedCreator";
  MuonPatternRecoDumper debug;

  auto output = std::make_unique<L2MuonTrajectorySeedCollection>();

  auto const l1TkMuColl = iEvent.getHandle(l1TkMuCollToken_);

  auto cscHandle = iEvent.getHandle(cscSegmentCollToken_);
  auto cscSegments = *cscHandle;
  auto dtHandle = iEvent.getHandle(dtSegmentCollToken_);
  auto dtSegments = *dtHandle;

  cscGeometry_ = iSetup.getHandle(cscGeometryToken_);
  dtGeometry_ = iSetup.getHandle(dtGeometryToken_);
  magneticField_ = iSetup.getHandle(magneticFieldToken_);

  //std::cout << "Number of muons " << l1TkMuColl->size() << std::endl;

  // Loop on all L1TkMu in event
  for (size_t l1TkMuIndex = 0; l1TkMuIndex != l1TkMuColl->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuColl, l1TkMuIndex);

    // Physical info of L1TkMu
    float pt = l1TkMuRef->phPt();
    //if (pt < minMomentum_) {
    //  continue;
    //}
    float eta = l1TkMuRef->phEta();
    float theta = 2 * std::atan(std::exp(-eta));
    float phi = l1TkMuRef->phPhi();
    int charge = l1TkMuRef->phCharge();

    //std::cout << "phEta: " << eta << ", " << "phPhi: " << phi << std::endl;
    Type muonType = overlap;
    if (std::abs(eta) < maxEtaBarrel_) {
      muonType = barrel;
    } else if (std::abs(eta) > maxEtaOverlap_) {
      muonType = endcap;
    }

    // Starting seed creation
    //std::cout << "Start seed creation" << '\n';

    l1t::MuonStubRefVector stubRefs = l1TkMuRef->stubs();

    //std::cout << "Number of stubs per L1TkMu: " << stubRefs.size() << '\n';
    //std::cout << "Number of DT segments in event: " << dtSegments.size() << '\n';
    //std::cout << "Number of CSC segments in event: " << cscSegments.size() << '\n';

    double bestDr2 = dRCone_ * dRCone_;
    int bestSegIndex = -1;
    bool bestInCsc = false;

    // Loop on L1TkMu stubs to find best association to DT/CSC segments
    for (auto stub : stubRefs) {
      //stub->print();

      // Separate barrel, endcap and overlap cases
      switch (muonType) {
        case barrel:
          if (!stub->isBarrel()) {
            continue;  // skip all non-barrel stubs
          }

          //std::cout << "BARREL checking " << dtSegments.size() << " DT segments" << '\n';
          for (size_t s = 0; s != dtSegments.size(); ++s) {
            // Check for stubs and segments with same detID
            // wheel->eta, station->depth, sector->phi
            DTChamberId stubId = DTChamberId(stub->etaRegion(),
                                             stub->depthRegion(),
                                             stub->phiRegion() + 1);  // converting online sectors to offline
            DTChamberId segId = dtSegments[s].chamberId();

            //std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
            //<< "segId: " << segId << ", RAW: " << segId.rawId() << '\n';

            if (!matchingDtIds(stubId, segId)) {
              continue;  // skip segments with different detID
            }
            //std::cout << "BARREL found match dtdetID" << '\n';

            // Match stub with DT segment
            GlobalPoint segPos = dtGeometry_->idToDet(segId)->toGlobal(dtSegments[s].localPosition());
            // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
            double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
            if (dr2 <= bestDr2) {
              //std::cout << "BARREL found match dR" << '\n';
              bestDr2 = dr2;
              bestSegIndex = s;
            }
          }
          //std::cout << "BARREL best segment: " << bestSegIndex << " found in csc? " << bestInCsc << '\n';
          break;

        case endcap:
          if (!stub->isEndcap()) {
            continue;  // skip all non-endcap stubs
          }

          //std::cout << "ENDCAP checking " << cscSegments.size() << " CSC segments" << '\n';
          for (size_t s = 0; s != cscSegments.size(); ++s) {
            // Check for stubs and segments with same detID
            // Endcap for CSC DetId (1 -> Forward, 2 -> Backward)
            int endcap = (eta > 0) ? 1 : 2;
            CSCDetId stubId = CSCDetId(endcap,
                                       stub->depthRegion(),
                                       6 - stub->etaRegion(),  //  mapping stub ring with csc ring
                                       stub->phiRegion());
            CSCDetId segId = cscSegments[s].cscDetId();

            //std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
            //<< "segId: " << segId << ", RAW: " << segId.rawId() << '\n';

            if (!(stubId == segId)) {
              continue;  // skip segments with different detIDs
            }
            //std::cout << "ENDCAP found match cscdetID" << '\n';

            // Match stub with CSC segment
            GlobalPoint segPos = cscGeometry_->idToDet(segId)->toGlobal(cscSegments[s].localPosition());
            // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
            double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
            if (dr2 <= bestDr2) {
              //std::cout << "ENDCAP found match dR" << '\n';
              bestDr2 = dr2;
              bestSegIndex = s;
              bestInCsc = true;
            }
          }
          //std::cout << "ENDCAP best segment: " << bestSegIndex << " found in csc? " << bestInCsc << '\n';
          break;

        case overlap:
          // Overlap runs both DT and CSC and picks the best match

          // Check DTs
          //std::cout << "OVERLAP checking " << dtSegments.size() << " DT segments" << '\n';
          for (size_t s = 0; s != dtSegments.size(); ++s) {
            if (!stub->isBarrel()) {
              continue;
            }
            DTChamberId stubId = DTChamberId(stub->etaRegion(), stub->depthRegion(), stub->phiRegion() + 1);
            DTChamberId segId = dtSegments[s].chamberId();

            //std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
            //<< "segId: " << segId << ", RAW: " << segId.rawId() << '\n';

            if (!matchingDtIds(stubId, segId)) {
              continue;
            }
            //std::cout << "OVERLAP found match dtDetID" << '\n';

            GlobalPoint segPos = dtGeometry_->idToDet(segId)->toGlobal(dtSegments[s].localPosition());
            // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
            double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
            if (dr2 <= bestDr2) {
              //std::cout << "OVERLAP found match dR in DT" << '\n';
              bestDr2 = dr2;
              bestSegIndex = s;
            }
          }

          // Check CSCs
          //std::cout << "OVERLAP checking also " << cscSegments.size() << " CSC segments" << '\n';
          for (size_t s = 0; s != cscSegments.size(); ++s) {
            if (!stub->isEndcap()) {
              continue;
            }
            int endcap = (l1TkMuRef->phEta() > 0) ? 1 : 2;  // endcap for CSC DetId (1 -> Forward, 2 -> Backward)
            CSCDetId stubId = CSCDetId(endcap, stub->depthRegion(), 6 - stub->etaRegion(), stub->phiRegion());
            CSCDetId segId = cscSegments[s].cscDetId();
            //std::cout << "stubId: " << stubId << ", RAW: " << stubId.rawId() << "\n"
            //<< "segId: " << segId << ", RAW: " << segId.rawId() << '\n';

            if (!(stubId == segId)) {
              continue;
            }
            //std::cout << "OVERLAP found match cscDetID" << '\n';

            GlobalPoint segPos = cscGeometry_->idToDet(segId)->toGlobal(cscSegments[s].localPosition());
            // IMPLEMENT check using phi (all hits have it), refine with z (if present) keep segment with most number of hits (closer in dR if ambiguous)
            double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
            if (dr2 <= bestDr2) {
              //std::cout << "OVERLAP found match dR in CSC" << '\n';
              bestDr2 = dr2;
              bestSegIndex = s;
              bestInCsc = true;
            }
          }
          //std::cout << "OVERLAP best segment: " << bestSegIndex << " found in csc? " << bestInCsc << '\n';
          break;

        default:
          std::cout << "L1TkMu must be either barrel, endcap or overlap" << std::endl;
          break;
      }
    }  // End loop on stubs

    // Emplace seeds in output
    if (bestSegIndex == -1) {
      continue;  // skip unmatched L1TkMu
    } else {
      // Info for propagation to MB2 or ME2
      service_->update(iSetup);
      const DetLayer* detLayer = nullptr;
      float radius = 0.;

      CLHEP::Hep3Vector vec(0., 1., 0.);
      vec.setTheta(theta);
      vec.setPhi(phi);

      DetId propagateToId;

      edm::OwnVector<TrackingRecHit> container;
      if (!bestInCsc) {
        // Found a matching segment in DT -> propagate to MB2
        //std::cout << "Found a matching segment in DTs, propagating to MB2 to seed" << '\n';
        // MB2
        propagateToId = DTChamberId(0, 2, 0);
        detLayer = service_->detLayerGeometry()->idToLayer(propagateToId);
        const BoundSurface* sur = &(detLayer->surface());
        const BoundCylinder* bc = dynamic_cast<const BoundCylinder*>(sur);
        radius = std::abs(bc->radius() / sin(theta));

        // Fill seed with matched segment(s)
        auto bestSegment = dtSegments[bestSegIndex];
        container.push_back(bestSegment);
      } else if (bestInCsc) {
        // Found a matching segment in CSC -> propagate to ME2
        //std::cout << "Found a matching segment in CSCs, propagating to ME2 to seed" << '\n';
        // ME2
        propagateToId = eta > 0 ? CSCDetId(1, 2, 0, 0, 0) : CSCDetId(2, 2, 0, 0, 0);
        detLayer = service_->detLayerGeometry()->idToLayer(propagateToId);
        radius = std::abs(detLayer->position().z() / cos(theta));

        // Fill seed with matched segment(s)
        auto bestSegment = cscSegments[bestSegIndex];
        container.push_back(bestSegment);
      }
      vec.setMag(radius);
      // Get Global point and direction
      GlobalPoint pos(vec.x(), vec.y(), vec.z());
      GlobalVector mom(pt * cos(phi), pt * sin(phi), pt * cos(theta) / sin(theta));

      GlobalTrajectoryParameters param(pos, mom, charge, &*magneticField_);

      AlgebraicSymMatrix55 mat;

      mat[0][0] = bestInCsc ? (0.4 / pt) * (0.4 / pt) : (0.25 / pt) * (0.25 / pt);  // sigma^2(charge/abs_momentum)
      mat[1][1] = 0.05 * 0.05;                                                      // sigma^2(lambda)
      mat[2][2] = 0.2 * 0.2;                                                        // sigma^2(phi)
      mat[3][3] = 20. * 20.;                                                        // sigma^2(x_transverse)
      mat[4][4] = 20. * 20.;                                                        // sigma^2(y_transverse)

      CurvilinearTrajectoryError error(mat);

      const FreeTrajectoryState state(param, error);

      // Create the TrajectoryStateOnSurface
      TrajectoryStateOnSurface tsos = service_->propagator(propagatorName_)->propagate(state, detLayer->surface());

      // Transform it in a Persistent TrajectoryStateOnDet
      PTrajectoryStateOnDet const& seedTSOS = trajectoryStateTransform::persistentState(tsos, propagateToId.rawId());

      // Emplace seed in output
      //std::cout << "Emplacing a DT seed in output" << '\n';
      output->emplace_back(L2MuonTrajectorySeed(seedTSOS, container, alongMomentum, l1TkMuRef));
    }  // End seed emplacing (one seed per L1TkMu)
  }  // End loop on L1TkMu
  //std::cout << "All L1TkMu in event processed" << '\n';
  iEvent.put(std::move(output));
}

DEFINE_FWK_MODULE(Phase2L2MuonSeedCreator);
