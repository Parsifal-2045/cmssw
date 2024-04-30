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
    : source_(pset.getParameter<edm::InputTag>("InputObjects")),
      muCollToken_(consumes(source_)),
      minMomentum_{pset.getParameter<double>("minimumSeedPt")},
      maxMomentum_{pset.getParameter<double>("maximumSeedPt")},
      defaultMomentum_{pset.getParameter<double>("defaultSeedPt")},
      debug{pset.getParameter<bool>("DebugMuonSeed")},
      sysError{pset.getParameter<double>("SeedPtSystematics")},
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

  auto const muColl = iEvent.getHandle(muCollToken_);
  LogDebug(metname) << "Number of muons " << muColl->size() << std::endl;

  iEvent.put(std::move(output));
}

/*
 * createSeed
 *
 * Note type = 1 --> CSC
 *           = 2 --> Overlap
 *           = 3 --> DT
 */
L2MuonTrajectorySeed Phase2L2MuonSeedCreator::createSeed(const int type,
                                                         const SegmentContainer& seg,
                                                         const l1t::TrackerMuon& l1TkMu,
                                                         const double& dRCone) {
  // The index of the station closest to the IP
  double l1Pt = l1TkMu.phPt();
  double sptmean = minMomentum_;
  int charge = l1TkMu.phCharge();

  if (l1Pt < minMomentum_) {
    l1Pt = minMomentum_;
    sptmean = minMomentum_;
  } else if (l1Pt > maxMomentum_) {
    l1Pt = maxMomentum_;
    sptmean = maxMomentum_ * 0.25;
  }

  l1t::MuonStubRefVector stubRefs = l1TkMu.stubs();

  AlgebraicSymMatrix mat(5, 0);
  double p_err = 0;

  int bestSeg = -1;

  edm::OwnVector<TrackingRecHit> container;

  L2MuonTrajectorySeed theSeed;

  if (type == 1) {
    // CSC
  }

  if (type == 2) {
    // Overlap
  }
  if (type == 3) {
    // DT
    // Cannot do get only available for DTRecSegmentCollection

    // for (auto stub : stubRefs) {
    // DTChamberId stubId = DTChamberId(stub->etaRegion(), stub->depthRegion(), stub->phiRegion());
    // auto [it, end] = seg->get(stubId);
    // Seg is MuonTransientTrackingRecHit
    // size_t s = 0;
    // for (; it != end; ++it, ++s) {
    //  GlobalPoint segPos = (*it)->globalPosition();
    //  // IMPLEMENT check best track matching (most number of hits)
    //  double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
    //  if (dr2 < dRCone * dRCone) {
    //    bestSeg = s;
    //    container.push_back(it->hit()->clone());
    //    break;
    //  }
    //  continue;
    //}
    for (auto stub : stubRefs) {
      for (size_t s = 0; s != seg.size(); ++s) {
        DTChamberId stubId = DTChamberId(stub->etaRegion(), stub->depthRegion(), stub->phiRegion());
        DetId segId = seg[s]->geographicalId();
        if (!(stubId.rawId() == segId.rawId())) {
          continue;
        }

        GlobalPoint segPos = seg[s]->globalPosition();
        // IMPLEMENT check best track matching (most number of hits)
        double dr2 = reco::deltaR2(segPos.eta(), segPos.phi(), stub->offline_eta1(), stub->offline_coord1());
        if (dr2 < dRCone * dRCone) {
          bestSeg = s;
          container.push_back(seg[s]->hit()->clone());
          break;
        }
      }
    }
    // Fill the local trajectory parameters
    LocalPoint segPos = seg[bestSeg]->localPosition();
    GlobalVector mom = seg[bestSeg]->globalPosition() - GlobalPoint();
    // Get the global direction
    GlobalVector polar(GlobalVector::Spherical(mom.theta(), seg[bestSeg]->globalDirection().phi(), 1.));
    // Scale magnitude of total momentum
    polar *= fabs(l1Pt) / polar.perp();
    // Trasfer into local direction
    LocalVector segDirFromPos = seg[bestSeg]->det()->toLocal(polar);
    LocalTrajectoryParameters param(segPos, segDirFromPos, charge);
    // Get error matrix
    mat = seg[bestSeg]->parametersError().similarityT(seg[bestSeg]->projectionMatrix());
    p_err = (sptmean * sptmean) / (polar.mag() * polar.mag() * l1Pt * l1Pt);
    mat[0][0] = p_err;
    LocalTrajectoryError error(asSMatrix<5>(mat));

    // Create the TrajectoryStateOnSurface
    TrajectoryStateOnSurface tsos(param, error, seg[bestSeg]->det()->surface(), &*BField);

    // Transform it in a TrajectoryStateOnSurface
    DetId segId = seg[bestSeg]->geographicalId();
    PTrajectoryStateOnDet seedTSOS = trajectoryStateTransform::persistentState(tsos, segId.rawId());
    auto dummyRef = edm::Ref<l1t::TrackerMuonCollection>();  // TrackerMuonRef is edm::Ref<TrackerMuonCollection>
    // edm::Ref<l1t::TrackerMuonCollection> mRef = l1TkMu.trkPtr(); // no conversion from edm::Ptr to edm::Ref
    theSeed = L2MuonTrajectorySeed(seedTSOS, container, alongMomentum, dummyRef);
  }
  return theSeed;
}

DEFINE_FWK_MODULE(Phase2L2MuonSeedCreator);
