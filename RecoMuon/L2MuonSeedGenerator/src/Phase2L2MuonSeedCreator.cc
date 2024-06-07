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
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CLHEP/Vector/ThreeVector.h"

//#define MY_LOG_DEBUG

#ifdef MY_LOG_DEBUG
std::mutex myMutex;
#define LOG(s)                                         \
  do {                                                 \
    std::lock_guard<std::mutex> lock(myMutex);         \
    std::cout << "(" << __LINE__ << ") " << s << '\n'; \
  } while (false)
#else
#define LOG(s) \
  do {         \
  } while (false)
#endif

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
      matchingPhiWindow_{pset.getParameter<double>("stubMatchDPhi")},
      matchingThetaWindow_{pset.getParameter<double>("stubMatchDTheta")},
      maxEtaBarrel_{pset.getParameter<double>("maximumEtaBarrel")},
      maxEtaOverlap_{pset.getParameter<double>("maximumEtaOverlap")},
      propagatorName_{pset.getParameter<string>("Propagator")} {
  // Service parameters
  edm::ParameterSet serviceParameters = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  // Services
  service_ = std::make_unique<MuonServiceProxy>(serviceParameters, consumesCollector());
  produces<L2MuonTrajectorySeedCollection>();
  produces<TrajectorySeedCollection>("phase2Validation");
}

void Phase2L2MuonSeedCreator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("InputObjects", edm::InputTag("l1tTkMuonsGmt"));
  desc.add<edm::InputTag>("CSCRecSegmentLabel", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("DTRecSegmentLabel", edm::InputTag("hltDt4DSegments"));
  desc.add<double>("MinPL1Tk", 3.5);
  desc.add<double>("MaxPL1Tk", 200);
  desc.add<double>("stubMatchDPhi", 0.05);
  desc.add<double>("stubMatchDTheta", 0.1);
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
  auto validationOutput = std::make_unique<TrajectorySeedCollection>();

  auto const l1TkMuColl = iEvent.getHandle(l1TkMuCollToken_);

  auto cscHandle = iEvent.getHandle(cscSegmentCollToken_);
  auto cscSegments = *cscHandle;
  auto dtHandle = iEvent.getHandle(dtSegmentCollToken_);
  auto dtSegments = *dtHandle;

  cscGeometry_ = iSetup.getHandle(cscGeometryToken_);
  dtGeometry_ = iSetup.getHandle(dtGeometryToken_);
  magneticField_ = iSetup.getHandle(magneticFieldToken_);

  LOG("Number of muons in Event: " << l1TkMuColl->size());

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

    LOG("phEta: " << eta << ", " << "phPhi: " << phi);
    Type muonType = overlap;
    if (std::abs(eta) < maxEtaBarrel_) {
      muonType = barrel;
      LOG("L1TkMu found in the barrel");
    } else if (std::abs(eta) > maxEtaOverlap_) {
      muonType = endcap;
      LOG("L1TkMu found in the endcap");
    }

    // Starting seed creation
    LOG("Start seed creation");

    l1t::MuonStubRefVector stubRefs = l1TkMuRef->stubs();

    LOG("Number of stubs per L1TkMu: " << stubRefs.size());
    LOG("Number of DT segments in event: " << dtSegments.size());
    LOG("Number of CSC segments in event: " << cscSegments.size());

    // Pairs segIndex, segQuality for matches in Barrel/Overlap/Endcap
    std::pair<int, int> matchesInBarrel[4] = {
        std::make_pair(-1, -1), std::make_pair(-1, -1), std::make_pair(-1, -1), std::make_pair(-1, -1)};
    std::pair<int, int> matchesInEndcap[4] = {
        std::make_pair(-1, -1), std::make_pair(-1, -1), std::make_pair(-1, -1), std::make_pair(-1, -1)};
    bool atLeastOneMatch = false;
    bool bestInDt = false;

    // Variables for matches in Overlap
    int totalBarrelQuality = 0;
    int totalEndcapQuality = 0;
    unsigned int nDtHits = 0;
    unsigned int nCscHits = 0;

    // Loop on L1TkMu stubs to find best association to DT/CSC segments
    for (auto stub : stubRefs) {
#ifdef MY_LOG_DEBUG
      stub->print();
#endif
      // Separate barrel, endcap and overlap cases
      switch (muonType) {
        case barrel: {
          if (!stub->isBarrel()) {
            continue;  // skip all non-barrel stubs
          }
          // Create detId for stub
          DTChamberId stubId = DTChamberId(stub->etaRegion(),       // wheel
                                           stub->depthRegion(),     // station
                                           stub->phiRegion() + 1);  // sector, online to offline

          // FIXME_ currently only considers the last stub per each station
          // could be not the best match if multiple stubs are found in the same station
          matchesInBarrel[stubId.station() - 1] = matchingStubSegment(stubId, stub, dtSegments, l1TkMuRef);
          if (matchesInBarrel[stubId.station() - 1].first != -1) {
            atLeastOneMatch = true;
            bestInDt = true;
          }
#ifdef MY_LOG_DEBUG
          LOG("BARREL best segments:");
          for (int i = 0; i != 4; ++i) {
            LOG("Station " << i + 1 << " (" << matchesInBarrel[i].first << ", " << matchesInBarrel[i].second << ")");
          }
#endif
          break;
        }  // End barrel

        case endcap: {
          if (!stub->isEndcap()) {
            continue;  // skip all non-endcap stubs
          }
          // Create detId for stub
          int endcap = (eta > 0) ? 1 : 2;  // CSC DetId endcap (1 -> Forward, 2 -> Backwards)
          CSCDetId stubId =
              CSCDetId(endcap,
                       stub->depthRegion(),              // station
                       6 - std::abs(stub->etaRegion()),  // ring, online to offline FIXME_ shouldn't need abs
                       stub->phiRegion());               // chamber

          // FIXME_ currently only considers the last stub per each station
          // could be not the best match if multiple stubs are found in the same station
          matchesInEndcap[stubId.station() - 1] = matchingStubSegment(stubId, stub, cscSegments, l1TkMuRef);
          if (matchesInEndcap[stubId.station() - 1].first != -1) {
            atLeastOneMatch = true;
          }
#ifdef MY_LOG_DEBUG
          LOG("ENDCAP best segments:");
          for (int i = 0; i != 4; ++i) {
            LOG("Station " << i + 1 << " (" << matchesInEndcap[i].first << ", " << matchesInEndcap[i].second << ")");
          }
#endif
          break;
        }  // End endcap

        case overlap: {
          // Overlap runs on both DTs and CSCs and picks the best overall match
          if (stub->isBarrel()) {
            // Check DTs
            LOG("OVERLAP stub in DTs, checking " << dtSegments.size() << " DT segments");
            // Create detId for stub
            DTChamberId stubId = DTChamberId(stub->etaRegion(),       // wheel
                                             stub->depthRegion(),     // station
                                             stub->phiRegion() + 1);  // sector, online to offline
            // FIXME_ currently only considers the last stub per each station
            // could be not the best match if multiple stubs are found in the same station
            matchesInBarrel[stubId.station() - 1] = matchingStubSegment(stubId, stub, dtSegments, l1TkMuRef);
            totalBarrelQuality += matchesInBarrel[stubId.station() - 1].second;
            auto tmpBestSegIndex = matchesInBarrel[stubId.station() - 1].first;
            if (tmpBestSegIndex != -1) {
              atLeastOneMatch = true;
              auto dtSegment = dtSegments.begin() + tmpBestSegIndex;
              nDtHits += (dtSegment->hasPhi() ? dtSegment->phiSegment()->recHits().size() : 0);
              nDtHits += (dtSegment->hasZed() ? dtSegment->zSegment()->recHits().size() : 0);
            }
#ifdef MY_LOG_DEBUG
            LOG("OVERLAP best segments in DTs:");
            for (int i = 0; i != 4; ++i) {
              LOG("Station " << i + 1 << " (" << matchesInBarrel[i].first << ", " << matchesInBarrel[i].second << ")");
            }
#endif
          } else if (stub->isEndcap()) {
            // Check CSCs
            LOG("OVERLAP stub in CSCs, checking " << cscSegments.size() << " CSC segments");
            int endcap = (eta > 0) ? 1 : 2;  // CSC DetId endcap (1 -> Forward, 2 -> Backwards)
            CSCDetId stubId =
                CSCDetId(endcap,
                         stub->depthRegion(),              // station
                         6 - std::abs(stub->etaRegion()),  // ring, online to offline FIXME_ shouldn't need abs
                         stub->phiRegion());               // chamber

            // FIXME_ currently only considers the last stub per each station
            // could be not the best match if multiple stubs are found in the same station
            matchesInEndcap[stubId.station() - 1] = matchingStubSegment(stubId, stub, cscSegments, l1TkMuRef);
            totalEndcapQuality += matchesInEndcap[stubId.station() - 1].second;
            auto tmpBestSegIndex = matchesInEndcap[stubId.station() - 1].first;
            if (tmpBestSegIndex != -1) {
              atLeastOneMatch = true;
              auto cscSegment = cscSegments.begin() + tmpBestSegIndex;
              nCscHits += cscSegment->nRecHits();
            }
#ifdef MY_LOG_DEBUG
            LOG("OVERLAP best segments in CSCs:");
            for (int i = 0; i != 4; ++i) {
              LOG("Station " << i + 1 << " (" << matchesInEndcap[i].first << ", " << matchesInEndcap[i].second << ")");
            }
#endif
          }

          LOG("OVERLAP comparing total qualities. DT: " << totalBarrelQuality << ", CSC: " << totalEndcapQuality);

          // Pick segment(s) with better quality
          bestInDt = (totalBarrelQuality > totalEndcapQuality) ? true : false;

          // Same qualities, pick higher number of hits
          if (totalBarrelQuality == totalEndcapQuality and totalBarrelQuality > -1) {
            LOG("Same quality " << totalBarrelQuality << ". Checking total number of hits");
            LOG("DT hits: " << nDtHits << ", CSC hits: " << nCscHits);
            LOG((nDtHits > nCscHits ? "More hits in DT segment" : "More hits in CSC segment"));
            bestInDt = (nDtHits >= nCscHits) ? true : false;
          }

#ifdef MY_LOG_DEBUG
          LOG("OVERLAP best segments:");
          if (bestInDt) {
            LOG("OVERLAP best match in DTs:");
            for (int i = 0; i != 4; ++i) {
              if (matchesInBarrel[i].first != -1) {
                LOG("Station " << i + 1 << " (" << matchesInBarrel[i].first << ", " << matchesInBarrel[i].second
                               << ")");
              }
            }
          } else if (!bestInDt) {
            LOG("OVERLAP best match in CSCs:");
            for (int i = 0; i != 4; ++i) {
              if (matchesInEndcap[i].first != -1) {
                LOG("Station " << i + 1 << " (" << matchesInEndcap[i].first << ", " << matchesInEndcap[i].second
                               << ")");
              }
            }
          }
#endif
          break;
        }  // End overlap

        default:
          std::cout << "L1TkMu must be either barrel, endcap or overlap" << std::endl;
          break;
      }
    }  // End loop on stubs

    // Emplace seeds in output
    if (!atLeastOneMatch) {
      LOG("No matching stub found, skipping seed");
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
      if (bestInDt) {
        // Found at least one matching segment in DT -> propagate to MB2
        LOG("Found matching segment(s) in DTs, propagating L1TkMu info to MB2 to seed");
        // MB2
        propagateToId = DTChamberId(0, 2, 0);
        detLayer = service_->detLayerGeometry()->idToLayer(propagateToId);
        const BoundSurface* sur = &(detLayer->surface());
        const BoundCylinder* bc = dynamic_cast<const BoundCylinder*>(sur);
        radius = std::abs(bc->radius() / sin(theta));

        // Fill seed with matched segment(s)
        for (int i = 0; i != 4; ++i) {
          if (matchesInBarrel[i].first != -1) {
            LOG("Adding matched DT segment in station " << i + 1 << " to the seed");
            container.push_back(dtSegments[matchesInBarrel[i].first]);
          }
        }
      } else if (!bestInDt) {
        // Found a matching segment in CSC -> propagate to ME2
        LOG("Found matching segment(s) in CSCs, propagating L1TkMu info to ME2 to seed");
        // ME2
        propagateToId = eta > 0 ? CSCDetId(1, 2, 0, 0, 0) : CSCDetId(2, 2, 0, 0, 0);
        detLayer = service_->detLayerGeometry()->idToLayer(propagateToId);
        radius = std::abs(detLayer->position().z() / cos(theta));

        // Fill seed with matched segment(s)
        for (int i = 0; i != 4; ++i) {
          if (matchesInEndcap[i].first != -1) {
            LOG("Adding matched CSC segment in station " << i + 1 << " to the seed");
            container.push_back(cscSegments[matchesInEndcap[i].first]);
          }
        }
      }
      vec.setMag(radius);
      // Get Global point and direction
      GlobalPoint pos(vec.x(), vec.y(), vec.z());
      GlobalVector mom(pt * cos(phi), pt * sin(phi), pt * cos(theta) / sin(theta));

      GlobalTrajectoryParameters param(pos, mom, charge, &*magneticField_);

      AlgebraicSymMatrix55 mat;

      mat[0][0] = bestInDt ? (0.25 / pt) * (0.25 / pt) : (0.4 / pt) * (0.4 / pt);  // sigma^2(charge/abs_momentum)
      mat[1][1] = 0.05 * 0.05;                                                     // sigma^2(lambda)
      mat[2][2] = 0.2 * 0.2;                                                       // sigma^2(phi)
      mat[3][3] = 20. * 20.;                                                       // sigma^2(x_transverse)
      mat[4][4] = 20. * 20.;                                                       // sigma^2(y_transverse)

      CurvilinearTrajectoryError error(mat);

      const FreeTrajectoryState state(param, error);

      // Create the TrajectoryStateOnSurface
      TrajectoryStateOnSurface tsos = service_->propagator(propagatorName_)->propagate(state, detLayer->surface());

      // Transform it in a Persistent TrajectoryStateOnDet
      const PTrajectoryStateOnDet& seedTSOS = trajectoryStateTransform::persistentState(tsos, propagateToId.rawId());

      // Emplace seed in output
      LOG("Emplacing seed in output");
      output->emplace_back(L2MuonTrajectorySeed(seedTSOS, container, alongMomentum, l1TkMuRef));
      validationOutput->emplace_back(TrajectorySeed(seedTSOS, container, alongMomentum));
    }  // End seed emplacing (one seed per L1TkMu)
  }  // End loop on L1TkMu
  LOG("All L1TkMu in event processed");
  iEvent.put(std::move(output));
  iEvent.put(std::move(validationOutput), "phase2Validation");
}

const bool Phase2L2MuonSeedCreator::matchingDtIds(const DTChamberId& stubId, const DTChamberId& segId) const {
  if (stubId.sector() == 4 or stubId.sector() == 10) {
    if (stubId.sector() == 4 and (segId.sector() == 4 or segId.sector() == 13)) {
      return (stubId.wheel() == segId.wheel() and stubId.station() == segId.station());
    }
    if (stubId.sector() == 10 and (segId.sector() == 10 or segId.sector() == 14)) {
      return (stubId.wheel() == segId.wheel() and stubId.station() == segId.station());
    }
  }
  return stubId == segId;
}

// Pair bestSegIndex, quality for DT segments matching
const std::pair<int, int> Phase2L2MuonSeedCreator::matchingStubSegment(const DTChamberId& stubId,
                                                                       const l1t::MuonStubRef stub,
                                                                       const DTRecSegment4DCollection& segments,
                                                                       const l1t::TrackerMuonRef l1TkMuRef) const {
  int bestSegIndex = -1;
  unsigned int nHitsPhiBest = 0;
  unsigned int nHitsThetaBest = 0;
  int quality = -1;

  LOG("Matching stub with DT segment");
  int matchingIds = 0;
  for (DTRecSegment4DCollection::const_iterator segment = segments.begin(), last = segments.end(); segment != last;
       ++segment) {
    DTChamberId segId = segment->chamberId();
    if (!matchingDtIds(stubId, segId)) {
      continue;  // skip segments with different detIds
    }
    ++matchingIds;
    // Global position of the segment
    GlobalPoint segPos = dtGeometry_->idToDet(segId)->toGlobal(segment->localPosition());

    // Check delta phi
    double deltaPhi = std::abs(segPos.phi() - stub->offline_coord1());
    LOG("deltaPhi: " << deltaPhi);
    if (deltaPhi > matchingPhiWindow_) {
      continue;
    }
    quality = 1;
    // Inside phi window -> check hit multiplicity (phi)
    unsigned int nHitsPhi = (segment->hasPhi() ? segment->phiSegment()->recHits().size() : 0);
    if (nHitsPhi == nHitsPhiBest) {
      // Same phi hit multiplicity -> check delta theta
      double deltaTheta = std::abs(segPos.theta() - 2 * std::atan(std::exp(-l1TkMuRef->phEta())));
      LOG("deltaTheta: " << deltaTheta);
      if (deltaTheta > matchingThetaWindow_) {
        continue;
      }
      quality = 2;
      // Inside theta window -> check hit multiplicity (theta)
      unsigned int nHitsTheta = (segment->hasZed() ? segment->zSegment()->recHits().size() : 0);
      if (nHitsTheta > nHitsThetaBest) {
        nHitsThetaBest = nHitsTheta;
        bestSegIndex = std::distance(segments.begin(), segment);
        quality = 4;
      }
    } else if (nHitsPhi > nHitsPhiBest) {
      nHitsPhiBest = nHitsPhi;
      bestSegIndex = std::distance(segments.begin(), segment);
      quality = 3;
    }
  }  // End loop on segments
  LOG("Looped over " << matchingIds << (matchingIds > 1 ? " segments" : " segment") << " with same DT detId as stub");
  LOG("Best DT segment: " << bestSegIndex << " with " << nHitsPhiBest + nHitsThetaBest << " hits and quality "
                          << quality);
  return std::make_pair(bestSegIndex, quality);
}

// Pair bestSegIndex, quality for CSC segments matching
const std::pair<int, int> Phase2L2MuonSeedCreator::matchingStubSegment(const CSCDetId& stubId,
                                                                       const l1t::MuonStubRef stub,
                                                                       const CSCSegmentCollection& segments,
                                                                       const l1t::TrackerMuonRef l1TkMuRef) const {
  int bestSegIndex = -1;
  unsigned int nHitsBest = 0;
  int quality = -1;
  int matchingIds = 0;

  LOG("Matching stub with CSC segment");

  for (CSCSegmentCollection::const_iterator segment = segments.begin(), last = segments.end(); segment != last;
       ++segment) {
    CSCDetId segId = segment->cscDetId();
    if (stubId != segId) {
      continue;  // skip segments with different detIds
    }
    ++matchingIds;

    // Global position of the segment
    GlobalPoint segPos = cscGeometry_->idToDet(segId)->toGlobal(segment->localPosition());

    // Check delta phi
    double deltaPhi = std::abs(segPos.phi() - stub->offline_coord1());
    LOG("deltaPhi: " << deltaPhi);
    if (deltaPhi > matchingPhiWindow_) {
      continue;
    }
    quality = 1;
    // Inside phi window -> check hit multiplicity (phi)
    unsigned int nHits = segment->nRecHits();
    if (nHits == nHitsBest) {
      // Same phi hit multiplicity -> check delta theta
      double deltaTheta = std::abs(segPos.theta() - 2 * std::atan(std::exp(-l1TkMuRef->phEta())));
      LOG("deltaTheta: " << deltaTheta);
      if (deltaTheta > matchingThetaWindow_) {
        continue;
      }
      // Inside theta window -> update bestSegment
      bestSegIndex = std::distance(segments.begin(), segment);
      quality = 3;
    } else if (nHits > nHitsBest) {
      nHitsBest = nHits;
      bestSegIndex = std::distance(segments.begin(), segment);
      quality = 4;
    }
  }  // End loop on segments
  LOG("Looped over " << matchingIds << (matchingIds > 1 ? " segments" : " segment") << " with same CSC detId as stub");
  LOG("Best CSC segment: " << bestSegIndex << " with " << nHitsBest << " hits and quality " << quality);
  return std::make_pair(bestSegIndex, quality);
}

DEFINE_FWK_MODULE(Phase2L2MuonSeedCreator);
