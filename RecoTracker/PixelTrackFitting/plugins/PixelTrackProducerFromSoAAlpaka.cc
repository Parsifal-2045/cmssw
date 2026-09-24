#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/SiPixelClusterSoA/interface/ClusteringConstants.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackSoA/interface/TracksHost.h"
#include "DataFormats/TrackSoA/interface/alpaka/TrackUtilities.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "DataFormats/TrackingRecHitSoA/interface/OTRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/StubsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/OTRecHitsHost.h"
#include "DataFormats/TrackingRecHitSoA/interface/StubsHost.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/PixelSeeding/interface/OTHitTag.h"
#include "RecoTracker/PixelTrackFitting/interface/alpaka/FitUtils.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "TrackingTools/AnalyticalJacobians/interface/AnalyticalCurvilinearJacobian.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCurvilinear.h"
#include "TrackingTools/GeomPropagators/interface/OptimalHelixPlaneCrossing.h"
#include "TrackingTools/GeomPropagators/interface/StraightLinePlaneCrossing.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "storeTracks.h"

/**
 * This class creates "legacy" reco::Track
 * objects from the output of SoA CA.
 */

// #define GPU_DEBUG
// struct that holds two maps for detIds of the OT modules
struct DetIdMaps {
  DetIdMaps() : detIdToOTModuleId_(), detIdIsUsedOTModule_() {}

  // map from the detId of OT modules to the moduleId among the used OT modules
  // (starting from 0 for first module of first OT layer)
  std::map<uint32_t, uint32_t> detIdToOTModuleId_;
  // map from detId to bool if used as OT extension
  std::map<uint32_t, bool> detIdIsUsedOTModule_;
};

namespace {
  // State of the helix of a free state where it crosses a module plane along the momentum: what
  // AnalyticalPropagator(alongMomentum) returns for a plane, with its straight-line case and its limit on the
  // turning angle, without building a TrajectoryStateOnSurface and without its field lookup at the destination.
  struct PlaneState {
    GlobalPoint position;
    GlobalVector momentum;
    AlgebraicSymMatrix55 covariance;
    bool ok = false;
  };

  PlaneState helixStateOnPlane(FreeTrajectoryState const &fts, Plane const &plane) {
    constexpr float maxDPhi = 1.6f;  // AnalyticalPropagator's default limit on the turning angle
    PlaneState state;
    float const rho = fts.transverseCurvature();
    double s;
    if (std::abs(rho) < 1.e-10f) {
      StraightLinePlaneCrossing crossing(StraightLinePlaneCrossing::PositionType(fts.position()),
                                         StraightLinePlaneCrossing::DirectionType(fts.momentum()),
                                         alongMomentum);
      auto const [ok, path] = crossing.pathLength(plane);
      if (not ok)
        return state;
      s = path;
      state.position = GlobalPoint(crossing.position(s));
      state.momentum = fts.momentum();
    } else {
      OptimalHelixPlaneCrossing crossing(plane,
                                         HelixPlaneCrossing::PositionType(fts.position()),
                                         HelixPlaneCrossing::DirectionType(fts.momentum()),
                                         rho,
                                         alongMomentum);
      auto const [ok, path] = (*crossing).pathLength(plane);
      if (not ok)
        return state;
      s = path;
      float const dphi = float(s) * rho;
      if (dphi * dphi * fts.momentum().perp2() > maxDPhi * maxDPhi * fts.momentum().mag2())
        return state;
      state.position = GlobalPoint((*crossing).position(s));
      GlobalVector const direction((*crossing).direction(s));
      state.momentum = direction * (fts.momentum().mag() / direction.mag());
    }
    AnalyticalCurvilinearJacobian const jacobian(fts.parameters(), state.position, state.momentum, s);
    state.covariance = ROOT::Math::Similarity(jacobian.jacobian(), fts.curvilinearError().matrix());
    state.ok = true;
    return state;
  }

  // TrackExtra state on the plane of a hit, from the helix of the fit. A helix that does not reach the plane
  // (loopers) is replaced by the straight line of the momentum, with the untransported covariance and the ok
  // flag down: the state still lies on the plane of the hit, so the consumers that start from it on that surface
  // (trajectoryStateTransform::inner/outerStateOnSurface, TrackTransformer) never start off the surface.
  struct ExtraState {
    reco::TrackExtra::Point position;
    reco::TrackExtra::Vector momentum;
    reco::TrackExtra::CovarianceMatrix covariance;
    unsigned int detId;
    bool ok;
  };

  ExtraState extraStateOnHit(FreeTrajectoryState const &fts, TrackingRecHit const &hit) {
    Plane const &plane = hit.det()->surface();
    PlaneState state = helixStateOnPlane(fts, plane);
    if (not state.ok) {
      StraightLinePlaneCrossing crossing(StraightLinePlaneCrossing::PositionType(fts.position()),
                                         StraightLinePlaneCrossing::DirectionType(fts.momentum()),
                                         anyDirection);
      auto const [ok, path] = crossing.pathLength(plane);
      state.position = ok ? GlobalPoint(crossing.position(path)) : hit.globalPosition();
      state.momentum = fts.momentum();
      state.covariance = fts.curvilinearError().matrix();
    }
    return ExtraState{reco::TrackExtra::Point(state.position.x(), state.position.y(), state.position.z()),
                      reco::TrackExtra::Vector(state.momentum.x(), state.momentum.y(), state.momentum.z()),
                      state.covariance,
                      hit.geographicalId().rawId(),
                      state.ok};
  }
}  // namespace

class PixelTrackProducerFromSoAAlpaka : public edm::global::EDProducer<edm::RunCache<DetIdMaps>> {
  using TrackSoAHost = reco::TracksHost;
  using HMSstorage = std::vector<uint32_t>;
  using IndToEdm = std::vector<uint32_t>;
  using TrackHitSoA = reco::TrackHitSoA;

public:
  explicit PixelTrackProducerFromSoAAlpaka(const edm::ParameterSet &iConfig);
  ~PixelTrackProducerFromSoAAlpaka() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  std::shared_ptr<DetIdMaps> globalBeginRun(edm::Run const &, edm::EventSetup const &) const override;
  void globalEndRun(edm::Run const &, edm::EventSetup const &) const override {};

private:
  void produce(edm::StreamID streamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const override;

  // Event Data tokens
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::EDGetTokenT<TrackSoAHost> trackSoAToken_;
  const edm::EDGetTokenT<SiPixelRecHitCollectionNew> pixelRecHitsToken_;
  edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> otRecHitsToken_;
  const edm::EDGetTokenT<HMSstorage> pixelHMSToken_;
  edm::EDGetTokenT<HMSstorage> otHMSToken_;
  edm::EDGetTokenT<reco::OTRecHitsHost> otRecHitsSoAToken_;
  edm::EDGetTokenT<reco::StubsHost> stubsSoAToken_;
  // Event Setup tokens
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> idealMagneticFieldToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryTokenRun_;

  int32_t const minNumberOfHits_;
  pixelTrack::Quality const minQuality_;
  const bool useOTExtension_;
  const bool throwOnMissing_;
  const bool expandStubs_;
  const bool requireQuadsFromConsecutiveLayers_;
  const bool setAlgorithmFromIteration_;
  const bool fillTrackExtra_;
  const bool verbose_;
};

PixelTrackProducerFromSoAAlpaka::PixelTrackProducerFromSoAAlpaka(const edm::ParameterSet &iConfig)
    : beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      trackSoAToken_(consumes(iConfig.getParameter<edm::InputTag>("trackSrc"))),
      pixelRecHitsToken_(
          consumes<SiPixelRecHitCollectionNew>(iConfig.getParameter<edm::InputTag>("pixelRecHitLegacySrc"))),
      pixelHMSToken_(consumes<HMSstorage>(iConfig.getParameter<edm::InputTag>("pixelRecHitLegacySrc"))),
      idealMagneticFieldToken_(esConsumes()),
      trackerTopologyToken_(esConsumes()),
      trackerGeometryTokenRun_(esConsumes<edm::Transition::BeginRun>()),
      minNumberOfHits_(iConfig.getParameter<int>("minNumberOfHits")),
      minQuality_(pixelTrack::qualityByName(iConfig.getParameter<std::string>("minQuality"))),
      useOTExtension_(iConfig.getParameter<bool>("useOTExtension")),
      throwOnMissing_(iConfig.getParameter<bool>("throwOnMissing")),
      expandStubs_(iConfig.getParameter<bool>("expandStubs")),
      requireQuadsFromConsecutiveLayers_(iConfig.getParameter<bool>("requireQuadsFromConsecutiveLayers")),
      setAlgorithmFromIteration_(iConfig.getParameter<bool>("setAlgorithmFromIteration")),
      fillTrackExtra_(iConfig.getParameter<bool>("fillTrackExtra")),
      verbose_(iConfig.getUntrackedParameter<bool>("verbose")) {
  if (minQuality_ == pixelTrack::Quality::notQuality) {
    throw cms::Exception("PixelTrackConfiguration")
        << iConfig.getParameter<std::string>("minQuality") + " is not a pixelTrack::Quality";
  }
  if (minQuality_ < pixelTrack::Quality::dup) {
    throw cms::Exception("PixelTrackConfiguration")
        << iConfig.getParameter<std::string>("minQuality") + " not supported";
  }
  produces<TrackingRecHitCollection>();
  produces<reco::TrackExtraCollection>();
  // TrackCollection refers to TrackingRechit and TrackExtra
  // collections, need to declare its production after them to work
  // around a rare race condition in framework scheduling
  produces<reco::TrackCollection>();
  produces<IndToEdm>();

  // if useOTExtension consume the OT RecHits
  if (useOTExtension_) {
    otRecHitsToken_ =
        consumes<Phase2TrackerRecHit1DCollectionNew>(iConfig.getParameter<edm::InputTag>("outerTrackerRecHitSrc"));
    otHMSToken_ = consumes<HMSstorage>(iConfig.getParameter<edm::InputTag>("outerTrackerRecHitSoAConverterSrc"));
  }

  // if expandStubs consume the OTRecHitsSoA and StubsSoA collections
  if (expandStubs_) {
    otRecHitsSoAToken_ = consumes<reco::OTRecHitsHost>(iConfig.getParameter<edm::InputTag>("otRecHitsSoASrc"));
    stubsSoAToken_ = consumes<reco::StubsHost>(iConfig.getParameter<edm::InputTag>("stubsSoASrc"));
  }
}

std::shared_ptr<DetIdMaps> PixelTrackProducerFromSoAAlpaka::globalBeginRun(const edm::Run &iRun,
                                                                           const edm::EventSetup &iSetup) const {
  // make the maps object
  auto detIdMaps = std::make_shared<DetIdMaps>();

  // if OT RecHits are used in PixelTracks, fill the detIdToOTModuleId_ map
  if (useOTExtension_) {
    // get track geometry
    const auto &trackerGeometry = &iSetup.getData(trackerGeometryTokenRun_);

    // function to check if given module is used as OT for CA
    auto isPinPSinOTBarrel = [&](DetId detId) {
      // Select only P-hits from the OT barrel
      return (trackerGeometry->getDetectorType(detId) == TrackerGeometry::ModuleType::Ph2PSP &&
              detId.subdetId() == StripSubdetector::TOB);
    };

    // loop over all modules and fill the map detIdToOTModuleId_
    auto const &detUnits = trackerGeometry->detUnits();
    for (uint32_t otModuleId{0}; auto &detUnit : detUnits) {
      DetId detId(detUnit->geographicalId());
      // check if the module is used for OT extension
      bool isUsedOTModule = isPinPSinOTBarrel(detId);
      detIdMaps->detIdIsUsedOTModule_[detUnit->geographicalId()] = isUsedOTModule;
      if (isUsedOTModule) {
        // save the module index among the extension modules
        detIdMaps->detIdToOTModuleId_[detUnit->geographicalId()] = otModuleId;
        otModuleId++;
      }
    }
  }

  return detIdMaps;
}

void PixelTrackProducerFromSoAAlpaka::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("trackSrc", edm::InputTag("pixelTracksAlpaka"));
  desc.add<edm::InputTag>("pixelRecHitLegacySrc", edm::InputTag("siPixelRecHitsPreSplittingLegacy"));
  desc.add<edm::InputTag>("outerTrackerRecHitSrc", edm::InputTag("hltSiPhase2RecHits"));
  desc.add<edm::InputTag>("outerTrackerRecHitSoAConverterSrc", edm::InputTag("phase2OTRecHitsSoAConverter"));
  desc.add<edm::InputTag>("otRecHitsSoASrc", edm::InputTag("pixelSeedingOTRecHitsSoA"));
  desc.add<edm::InputTag>("stubsSoASrc", edm::InputTag("otStubProducer"));
  desc.add<int>("minNumberOfHits", 0);
  desc.add<std::string>("minQuality", "loose");
  desc.add<bool>("useOTExtension", false);
  desc.add<bool>("throwOnMissing", true)
      ->setComment(
          "Throw when the track SoA is absent; false makes the validation clones write empty collections for the "
          "events "
          "whose HLT paths did not run the pixel tracking");
  desc.add<bool>("expandStubs", false);
  desc.add<bool>("setAlgorithmFromIteration", false)
      ->setComment(
          "Stamp each reco::Track with the algorithm name of the CA iteration that found it, taken from the SoA "
          "iteration column: promptHighPt -> hltPixel, promptLowPt -> lowPtTripletStep, displaced -> "
          "displacedGeneralStep. None of the three names is produced anywhere else in the Phase-2 HLT menu, so the "
          "algorithm word splits a merged collection by the iteration it came from (validation labels). False (the "
          "default) leaves the algorithm at undefAlgorithm.");
  desc.add<bool>("fillTrackExtra", false)
      ->setComment(
          "Fill the inner and outer states of each TrackExtra (position, momentum, curvilinear covariance and detId "
          "at the first and last hit), from the helix of the fitted perigee state where it crosses the planes of the "
          "two hits; the outer helix bends in Bz averaged along the chord to the last hit (perigee, midpoint and last "
          "hit, weights 1:2:1). Needed when the tracks feed consumers that start from those states without a refit, "
          "such as the global muon matching and refit (L3MuonProducer), the muon identification and PF. No material "
          "is applied between the perigee and the hits: the outer momentum is the perigee momentum (no energy loss), "
          "and the covariance is the perigee covariance transported along the helix, a conservative uncertainty of "
          "the extrapolated state and not the covariance of a state smoothed with the hits (much tighter at the last "
          "hit). False (the default) leaves the TrackExtras without states.");

  // this option for removing tracks with exactly 4 hits is a temporary solution to reduce the fake rate in Phase-2
  // and is to be replaced by a smarter inclusive track selection in the CA directly
  desc.add<bool>("requireQuadsFromConsecutiveLayers", false);

  // Per-event hit-conservation diagnostic (OT-extra -> legacy rechit resolution). OFF by default.
  desc.addUntracked<bool>("verbose", false);

  descriptions.addWithDefaultLabel(desc);
}

void PixelTrackProducerFromSoAAlpaka::produce(edm::StreamID streamID,
                                              edm::Event &iEvent,
                                              const edm::EventSetup &iSetup) const {
  // enum class Quality : uint8_t { bad = 0, edup, dup, loose, strict, tight, highPurity };
  reco::TrackBase::TrackQuality recoQuality[] = {reco::TrackBase::undefQuality,
                                                 reco::TrackBase::undefQuality,
                                                 reco::TrackBase::discarded,
                                                 reco::TrackBase::loose,
                                                 reco::TrackBase::tight,
                                                 reco::TrackBase::tight,
                                                 reco::TrackBase::highPurity};
  assert(reco::TrackBase::highPurity == recoQuality[int(pixelTrack::Quality::highPurity)]);

  // CA iteration -> reco::Track algorithm, used only when setAlgorithmFromIteration_ is on. The three
  // names are existing TrackAlgorithm values that no other Phase-2 HLT module produces, so a merged
  // collection can be split by the iteration each track came from. A track united across the two arms
  // by the merger keeps the winning arm's iteration (the merger copies the column, and a twin winner
  // that absorbs its sibling keeps its own state), so it carries the winner's algorithm.
  // enum class Iteration : uint8_t { promptHighPt, promptLowPt, displaced, notIteration };
  constexpr reco::TrackBase::TrackAlgorithm recoAlgo[] = {reco::TrackBase::hltPixel,
                                                          reco::TrackBase::lowPtTripletStep,
                                                          reco::TrackBase::displacedGeneralStep,
                                                          reco::TrackBase::undefAlgorithm};
  static_assert(std::size(recoAlgo) == pixelTrack::iterationSize + 1,
                "pixelTrack::Iteration changed size: update the iteration -> algorithm mapping");

#ifdef GPU_DEBUG
  std::cout << "Converting soa helix in reco tracks" << std::endl;
#endif

  // index map: trackId(in SoA) -> trackId(in legacy edm)
  auto indToEdmP = std::make_unique<IndToEdm>();
  auto &indToEdm = *indToEdmP;

  auto const &idealField = iSetup.getData(idealMagneticFieldToken_);

  // prepare container for legacy tracks
  pixeltrackfitting::TracksWithRecHits tracks;

  // get trackerTopology
  auto const &trackerTopology = iSetup.getData(trackerTopologyToken_);

  // get the maps for the detId of the OT modules
  auto const &detIdIsUsedOTModule = runCache(iEvent.getRun().index())->detIdIsUsedOTModule_;
  auto const &detIdToOTModuleId = runCache(iEvent.getRun().index())->detIdToOTModuleId_;

  // Validation clones run for every event, including those whose HLT paths did not run the
  // pixel tracking: without the track SoA they write empty collections.
  if (not throwOnMissing_ and not iEvent.getHandle(trackSoAToken_).isValid()) {
    iEvent.put(std::make_unique<TrackingRecHitCollection>());
    iEvent.put(std::make_unique<reco::TrackExtraCollection>());
    iEvent.put(std::make_unique<reco::TrackCollection>());
    iEvent.put(std::move(indToEdmP));
    return;
  }

  // get beamspot
  const auto &bsh = iEvent.get(beamSpotToken_);
  GlobalPoint bs(bsh.x0(), bsh.y0(), bsh.z0());

  // get the module's starting indices in the hit collection
  auto const &pixelHitsModuleStart = iEvent.get(pixelHMSToken_);

  // get Pixel RecHits
  auto const &pixelRecHitsDSV = iEvent.get(pixelRecHitsToken_);
  auto const &pixelRecHits = pixelRecHitsDSV.data();
  auto const nPixelHits = pixelRecHits.size();

  // get OT RecHits if needed
  size_t nOTHits = 0;
  const Phase2TrackerRecHit1DCollectionNew *otRecHitsDSV = nullptr;
  if (useOTExtension_) {
    otRecHitsDSV = &iEvent.get(otRecHitsToken_);
    nOTHits = otRecHitsDSV->dataSize();
  }

  size_t nTotalHits = nPixelHits + nOTHits;

  // get OTRecHitsSoA and StubsSoA if stub expansion is enabled
  const reco::OTRecHitsHost *otRecHitsSoAHost = nullptr;
  const reco::StubsHost *stubsSoAHost = nullptr;
  reco::OTRecHitsConstView otRecHitsSoAView;
  reco::StubsConstView stubsSoAView;
  int32_t offsetStubs = -1;

  if (expandStubs_) {
    otRecHitsSoAHost = &iEvent.get(otRecHitsSoAToken_);
    otRecHitsSoAView = otRecHitsSoAHost->const_view().otRecHits();

    stubsSoAHost = &iEvent.get(stubsSoAToken_);
    stubsSoAView = stubsSoAHost->const_view().stubs();

    // Stubs start after pixel hits in the merged collection
    offsetStubs = static_cast<int32_t>(nPixelHits);
  }

  // hitmap to go from a unique RecHit identifier to the RecHit in the legacy collection
  // (unique hit identifier is equivalent to the position of the hit in the RecHit SoA)
  std::vector<TrackingRecHit const *> hitmap;
  hitmap.resize(nTotalHits, nullptr);

  // loop over pixel RecHits to fill the hitmap
  for (auto const &pixelHit : pixelRecHits) {
    auto const &thit = static_cast<BaseTrackerRecHit const &>(pixelHit);
    auto const detI = thit.det()->index();
    auto const &clus = thit.firstClusterRef();
    assert(clus.isPixel());

    // get hit identifier as (hit offset of the module) + (hit index in this module)
    auto const idx = pixelHitsModuleStart[detI] + clus.pixelCluster().originalId();

    assert(nullptr == hitmap[idx]);
    hitmap[idx] = &pixelHit;
  }

  // if OT RecHits are used in PixelTracks, fill the hitmap also with those
  if (useOTExtension_) {
    if (expandStubs_ && otRecHitsSoAHost != nullptr) {
      // The OT hits in the SoA are organized by StackedModuleGeometry index, not by
      // detUnit->index(). Each SoA hit stores origRecHitIdx, the flat index into the legacy
      // Phase2TrackerRecHit1DCollectionNew assigned while iterating the DetSets in legacy order,
      // so each SoA hit maps straight to its legacy RecHit by that index.
      auto const &otData = otRecHitsDSV->data();
      auto otHitsView = otRecHitsSoAHost->const_view().otRecHits();
      uint32_t nOTHitsSoA = otHitsView.metadata().size();
      for (uint32_t i = 0; i < nOTHitsSoA; ++i) {
        uint32_t flatIdx = otHitsView[i].origRecHitIdx();
        assert(flatIdx < otData.size());
        hitmap[nPixelHits + i] = &otData[flatIdx];
      }
    } else {
      // Without stub expansion: OT hits organized by detUnit->index() for Ph2PSP TOB modules.
      // The RecHits in the SoA are ordered according to the detUnit->index()
      // of the respective OT module. For this reason, we need the map from the
      // detId to the moduleId among all used OT modules. This otModuleId corresponds
      // to the module's position in the otHitsModuleStart that we get from the event.

      // get the module's starting indices in the hit collection
      auto const &otHitsModuleStart = iEvent.get(otHMSToken_);

      // perform the exact same loop of how the SoA is initially filled with OT hits
      // and get the index by counting the hits (starting from the correpondign HitStartModule)
      for (auto const &detSet : *otRecHitsDSV) {
        auto detId = detSet.detId();

        // check if module is used in extension
        if (detIdIsUsedOTModule.find(detId)->second) {
          // get the corresponding otModuleId
          auto otModuleId = detIdToOTModuleId.find(detId)->second;

          // loop over the RecHits of the module and fill the hitmap
          for (int idx = otHitsModuleStart[otModuleId]; auto const &recHit : detSet) {
            assert(nullptr == hitmap[idx]);
            hitmap[idx] = &recHit;
            idx++;
          }
        }
      }
    }
  }

  // function that returns the number of skipped layers for a given pair of RecHits
  // for the case where the inner RecHit is in the pixel barrel.
  auto getNSkippedLayersInnerInBarrel = [&](const DetId &innerDetId,
                                            const DetId &outerDetId,
                                            const TrackingRecHit *innerRecHit) {
    int nSkippedLayers = 0;
    switch (outerDetId.subdetId()) {
      case PixelSubdetector::PixelBarrel:
        nSkippedLayers = trackerTopology.pxbLayer(outerDetId) - trackerTopology.pxbLayer(innerDetId) - 1;
        break;
      case PixelSubdetector::PixelEndcap:
        nSkippedLayers = trackerTopology.pxfDisk(outerDetId) - 1;  // -1 because first disk has Id 1
        break;
      case StripSubdetector::TOB:
        // if the inner RecHit is at the edge of the barrel layer, consider the jump to the first OT layer as no skip
        if (std::abs(innerRecHit->globalPosition().z()) > 17)
          nSkippedLayers = trackerTopology.getOTLayerNumber(outerDetId) - 1;  // -1 because first barrel has Id 1
        else
          nSkippedLayers = trackerTopology.getOTLayerNumber(outerDetId) + 4 - trackerTopology.pxbLayer(innerDetId) - 1;
        break;
      case StripSubdetector::TID:
        // Pixel barrel to OT endcap disk: transition region, no skipped layers
        nSkippedLayers = 0;
        break;
    }
    return nSkippedLayers;
  };

  // function that returns the number of skipped layers for a given pair of RecHits
  // for the case where the inner RecHit is in the pixel endcap.
  auto getNSkippedLayersInnerInEndcap = [&](const DetId &innerDetId, const DetId &outerDetId) {
    int nSkippedLayers = 0;
    switch (outerDetId.subdetId()) {
      case PixelSubdetector::PixelEndcap:
        nSkippedLayers = trackerTopology.pxfDisk(outerDetId) - trackerTopology.pxfDisk(innerDetId) - 1;
        break;
      case StripSubdetector::TOB:
        nSkippedLayers = trackerTopology.getOTLayerNumber(outerDetId) - 1;  // -1 because first disk has Id 1
        break;
      case StripSubdetector::TID:
        // Pixel endcap to OT endcap disk: transition region, no skipped layers
        nSkippedLayers = 0;
        break;
    }
    return nSkippedLayers;
  };

  // function that returns the number of skipped layers for a given pair of RecHits
  // for the case where the inner RecHit is in the OT (barrel or endcap).
  auto getNSkippedLayersInnerInOT = [&](const DetId &innerDetId, const DetId &outerDetId) {
    int nSkippedLayers = 0;
    if (innerDetId.subdetId() == StripSubdetector::TOB && outerDetId.subdetId() == StripSubdetector::TOB) {
      // Both in OT barrel: compute layer difference
      nSkippedLayers = trackerTopology.getOTLayerNumber(outerDetId) - trackerTopology.getOTLayerNumber(innerDetId) - 1;
    } else if (innerDetId.subdetId() == StripSubdetector::TID && outerDetId.subdetId() == StripSubdetector::TID) {
      // Both in OT endcap: compute disk difference (same side)
      int innerDisk = trackerTopology.tidWheel(innerDetId);
      int outerDisk = trackerTopology.tidWheel(outerDetId);
      nSkippedLayers = outerDisk - innerDisk - 1;
    }
    // Barrel-to-endcap or endcap-to-barrel transitions: 0 skipped layers
    return nSkippedLayers;
  };

  // function that returns the number of skipped layers for a given pair of RecHits
  // It works only for Phase-2, as this feature does not make sense for Phase-1 due to the smaller number of layers.
  // (needed for layer-skipping quadruplet rejection)
  auto getNSkippedLayers = [&](const TrackingRecHit *innerRecHit, const TrackingRecHit *outerRecHit) {
    // get detIds and subdetectors of the hits to determine their layers
    auto innerDetId = innerRecHit->geographicalId();
    auto outerDetId = outerRecHit->geographicalId();

    int nSkippedLayers = 0;

    switch (innerDetId.subdetId()) {
      case PixelSubdetector::PixelBarrel:
        nSkippedLayers = getNSkippedLayersInnerInBarrel(innerDetId, outerDetId, innerRecHit);
        break;
      case PixelSubdetector::PixelEndcap:
        nSkippedLayers = getNSkippedLayersInnerInEndcap(innerDetId, outerDetId);
        break;
      case StripSubdetector::TOB:
      case StripSubdetector::TID:
        nSkippedLayers = getNSkippedLayersInnerInOT(innerDetId, outerDetId);
        break;
    }
    return nSkippedLayers;
  };

  // The consumers that fit the hits of a track in sequence (KF refits, the global muon refit) propagate along
  // the momentum from one hit to the next and stop at the first hit behind the state, so the hits are stored in
  // the order the track crosses them:
  // - across layers the CA order is inside-out, except for the tracks the merger unites from twins, which carry
  //   the non-shared hits of the absorbed twin after the winner's. A backward step between module centres larger
  //   than kMaxBackStep (above the spread of the module centres of one layer) flags them, and their hits are
  //   reordered by their distance from the beam spot;
  // - within a layer (the two sensors of a stub, overlapping modules) the order is not defined. Each run of
  //   consecutive hits on the same layer is sorted by the distance of the module plane from the beam spot along
  //   its normal: the planes of one layer are nearly parallel, and the distance does not depend on where the hit
  //   lies on the sensor (a strip hit sits at the strip centre).
  // A twin-merged track can also carry the same measurement twice (a raw OT hit and a stub sensor hit, or two
  // stubs sharing a sensor cluster, possibly as distinct rechit objects): a hit sharing its cluster with a hit
  // already kept in the layer run is dropped, so that no fit counts the measurement twice.
  // Both twin-merge effects are fixed on the device only in the merger's refit copy of the hit list.
  constexpr float kMaxBackStep = 5.f;  // cm
  uint32_t nTracksReordered = 0;
  uint32_t nDuplicateHits = 0;
  auto moduleDistance = [&bs](TrackingRecHit const *hit) { return (hit->det()->position() - bs).mag(); };
  auto sameLayer = [&trackerTopology](TrackingRecHit const *a, TrackingRecHit const *b) {
    DetId const ia = a->geographicalId(), ib = b->geographicalId();
    return ia.subdetId() == ib.subdetId() and trackerTopology.layer(ia) == trackerTopology.layer(ib) and
           trackerTopology.side(ia) == trackerTopology.side(ib);
  };
  auto planeDistance = [&bs](TrackingRecHit const *hit) {
    auto const &surface = hit->det()->surface();
    return std::abs(surface.normalVector().dot(surface.position() - bs));
  };
  auto sortInCrossingOrder = [&](std::vector<const TrackingRecHit *> &trackHits) {
    if (trackHits.size() < 2)
      return;
    float previous = moduleDistance(trackHits.front());
    for (auto hit = trackHits.begin() + 1; hit != trackHits.end(); ++hit) {
      float const distance = moduleDistance(*hit);
      if (distance < previous - kMaxBackStep) {
        std::stable_sort(trackHits.begin(), trackHits.end(), [&bs](TrackingRecHit const *a, TrackingRecHit const *b) {
          return (a->globalPosition() - bs).mag2() < (b->globalPosition() - bs).mag2();
        });
        ++nTracksReordered;
        break;
      }
      previous = distance;
    }
    // compacting in place: the kept hits of a run are written at `kept`, never past the hit being read
    auto kept = trackHits.begin();
    for (auto first = trackHits.begin(); first != trackHits.end();) {
      auto const last = std::find_if(
          first + 1, trackHits.end(), [&](TrackingRecHit const *hit) { return not sameLayer(*first, hit); });
      if (last - first > 1)
        std::stable_sort(first, last, [&](TrackingRecHit const *a, TrackingRecHit const *b) {
          return planeDistance(a) < planeDistance(b);
        });
      auto const runBegin = kept;
      for (auto hit = first; hit != last; ++hit) {
        if (std::any_of(runBegin, kept, [&](TrackingRecHit const *k) {
              return k == *hit or k->sharesInput(*hit, TrackingRecHit::all);
            })) {
          ++nDuplicateHits;
          continue;
        }
        *kept++ = *hit;
      }
      first = last;
    }
    trackHits.erase(kept, trackHits.end());
  };

  std::vector<const TrackingRecHit *> hits;
  hits.reserve(5);  //TODO move to a configurable parameter?

  auto const &tsoa = iEvent.get(trackSoAToken_);
  auto const quality = tsoa.view().tracks().quality();
  auto const hitOffs = tsoa.view().tracks().hitOffsets();
  // Plain column accessor for pt, used by the sort comparator below: tsoa.view().tracks()[i].pt()
  // would build a full element proxy per comparison, and that proxy's constructor builds the
  // Eigen::Map members of the layout's two Eigen columns.
  auto const trackPt = tsoa.view().tracks().pt();
  auto const hitIdxs = tsoa.view().trackHits().id();
  auto nTracks = tsoa.view().tracks().nTracks();

  tracks.reserve(nTracks);

  int32_t nt = 0;

  // sort index by pt
  std::vector<int32_t> sortIdxs(nTracks);
  std::iota(sortIdxs.begin(), sortIdxs.end(), 0);
  // sort good-quality tracks by pt, keep bad-quality tracks at the bottom
  std::sort(sortIdxs.begin(), sortIdxs.end(), [&](int32_t const i1, int32_t const i2) {
    if (quality[i1] >= minQuality_ && quality[i2] >= minQuality_)
      return trackPt[i1] > trackPt[i2];
    else
      return quality[i1] > quality[i2];
  });

  indToEdm.resize(nTracks, -1);

  // A track-hit id with caOTHitTag::kOTHitTag set is a raw OT rechit attached by the extension
  // stage; its low bits are the OT SoA row. It resolves to a legacy Phase2TrackerRecHit1D via
  // the OT portion of the hitmap (hitmap[nPixelHits + row]), populated only on the expandStubs OT
  // path. Where that map is unavailable the tagged extra is dropped, never crashing.
  const bool otTagResolvable = useOTExtension_ && expandStubs_ && otRecHitsSoAHost != nullptr;
  // Per-event diagnostic tallies of the tagged-OT-extra branch (one-shot print below).
  uint32_t nOTExtrasResolved = 0, nOTExtrasDropped = 0;

  // TrackExtras with the inner and outer states (fillTrackExtra_), one per stored track: the helix of the fitted
  // perigee state at the first and last hit, without material (the field-profile, energy-loss and scattering
  // corrections of the fit are not applied between the perigee and the hits).
  reco::TrackExtraCollection extras;
  if (fillTrackExtra_)
    extras.reserve(nTracks);
  uint32_t nExtraStatesFailed = 0;

  // loop over (sorted) tracks
  for (const auto &it : sortIdxs) {
    auto nHits = reco::nHits(tsoa.view().tracks(), it);
    assert(nHits >= 3);
    auto q = quality[it];

    // apply cuts on quality and number of hits
    if (q < minQuality_)
      // since the tracks are sorted according to quality,
      // we can break after the first track with low quality
      break;
    if (nHits < minNumberOfHits_)  //move to nLayers?
      continue;

    auto start = (it == 0) ? 0 : hitOffs[it - 1];
    auto end = hitOffs[it];
    int nRemovedHits{0};
    int nExpandedHits{0};

    // First pass: count how many hits we'll have after stub expansion
    if (expandStubs_ && stubsSoAHost != nullptr && offsetStubs >= 0) {
      for (auto iHit = start; iHit < end; ++iHit) {
        auto hitIdx = hitIdxs[iHit];
        if (caOTHitTag::isOTId(hitIdx)) {
          // Tagged raw-OT extra: one legacy rechit, no stub expansion. Neither expanded nor removed
          // when resolvable, dropped otherwise, matching the fill pass.
          const uint32_t o = caOTHitTag::otIdx(hitIdx);
          if (!(otTagResolvable && (nPixelHits + o) < nTotalHits))
            nRemovedHits++;
          continue;
        }
        if (hitIdx < nTotalHits) {
          if (hitIdx >= static_cast<uint32_t>(offsetStubs)) {
            uint32_t stubIdx = hitIdx - offsetStubs;
            if (isStub(stubsSoAView, stubIdx)) {
              nExpandedHits++;  // Regular stub expands to 2 hits, so we add 1 more
            }
            // PHitOnly stubs have only 1 hit (inner), so no expansion needed
          }
        } else {
          nRemovedHits++;
        }
      }
    } else {
      for (auto iHit = start; iHit < end; ++iHit) {
        auto hitIdx = hitIdxs[iHit];
        if (caOTHitTag::isOTId(hitIdx)) {
          // Tagged OT extra: resolvable only on the expandStubs OT path (false here) -> dropped.
          const uint32_t o = caOTHitTag::otIdx(hitIdx);
          if (!(otTagResolvable && (nPixelHits + o) < nTotalHits))
            nRemovedHits++;
          continue;
        }
        if (hitIdx >= nTotalHits) {
          nRemovedHits++;
        }
      }
    }

    hits.resize(nHits - nRemovedHits + nExpandedHits);

    int hitOutputIdx = 0;
    for (auto iHit = start; iHit < end; ++iHit) {
      auto hitIdx = hitIdxs[iHit];
      if (caOTHitTag::isOTId(hitIdx)) {
        // Tagged raw-OT extra -> its legacy Phase2TrackerRecHit1D via the OT hitmap (same lookup as
        // a stub's lower/upper sensor hit: hitmap[nPixelHits + otSoARow]). Unresolvable tagged ids
        // are dropped (counted as removed above), keeping the hits vector correctly sized.
        const uint32_t o = caOTHitTag::otIdx(hitIdx);
        if (otTagResolvable && (nPixelHits + o) < nTotalHits) {
          hits[hitOutputIdx++] = hitmap[nPixelHits + o];
          ++nOTExtrasResolved;
        } else {
          ++nOTExtrasDropped;
        }
        continue;
      }
      if (hitIdx < nTotalHits) {
        if (expandStubs_ && stubsSoAHost != nullptr && offsetStubs >= 0 &&
            hitIdx >= static_cast<uint32_t>(offsetStubs)) {
          uint32_t stubIdx = hitIdx - offsetStubs;
          uint32_t lowerHitIdx = stubsSoAView[stubIdx].lowerHitIdx();

          hits[hitOutputIdx++] = hitmap[nPixelHits + lowerHitIdx];

          // Add outer sensor hit only if not PHitOnly (PHitOnly stubs have invalid upperHitIdx)
          if (isStub(stubsSoAView, stubIdx)) {
            uint32_t upperHitIdx = stubsSoAView[stubIdx].upperHitIdx();
            hits[hitOutputIdx++] = hitmap[nPixelHits + upperHitIdx];
          }
        } else {
          hits[hitOutputIdx++] = hitmap[hitIdx];
        }
      }
      // else: removed hits are skipped
    }
    sortInCrossingOrder(hits);

    end = end - nRemovedHits;

    // implement custome requirement for quadruplets coming from consecutive layers
    if (requireQuadsFromConsecutiveLayers_ && (nHits == 4)) {
      bool skipThisTrack{false};
      // loop over layer pairs and check if they skip
      for (auto iHit = start; iHit < end - 1 and size_t(iHit - start + 1) < hits.size(); ++iHit) {
        // if the inner (iHit-start) to outer (iHit-start+1) hit layer-change skips 1 or more
        // layers skipt the track
        if (getNSkippedLayers(hits[iHit - start], hits[iHit - start + 1]) > 0) {
          skipThisTrack = true;
          break;
        }
      }
      if (skipThisTrack) {
        indToEdm[it] = pixelTrack::skippedTrack;  // mark as skipped
        continue;
      }
    }

#ifdef CA_DEBUG
    std::cout << "track soa : " << it << " with hits: ";
    for (auto iHit = start; iHit < end; ++iHit)
      std::cout << hitIdxs[iHit] << " - ";
    std::cout << std::endl;
#endif

    // store the index of the SoA:
    // indToEdm[index_SoAtrack] -> index_edmTrack (if it exists)
    indToEdm[it] = nt;
    ++nt;

    // mind: this values are respect the beamspot!
    float chi2 = tsoa.view().tracks()[it].chi2();
    float phi = reco::phi(tsoa.view().tracks(), it);

    riemannFit::Vector5d ipar, opar;
    riemannFit::Matrix5d icov, ocov;
    reco::copyToDense<riemannFit::Vector5d, riemannFit::Matrix5d>(tsoa.view().tracks(), ipar, icov, it);
    riemannFit::transformToPerigeePlane(ipar, icov, opar, ocov);

    LocalTrajectoryParameters lpar(opar(0), opar(1), opar(2), opar(3), opar(4), 1.);
    AlgebraicSymMatrix55 m;
    for (int i = 0; i < 5; ++i)
      for (int j = i; j < 5; ++j)
        m(i, j) = ocov(i, j);

    float sp = std::sin(phi);
    float cp = std::cos(phi);
    Surface::RotationType rot(sp, -cp, 0, 0, 0, -1.f, cp, sp, 0);

    Plane impPointPlane(bs, rot);
    GlobalTrajectoryParameters gp(
        impPointPlane.toGlobal(lpar.position()), impPointPlane.toGlobal(lpar.momentum()), lpar.charge(), &idealField);
    JacobianLocalToCurvilinear jl2c(impPointPlane, lpar, idealField);

    AlgebraicSymMatrix55 mo = ROOT::Math::Similarity(jl2c.jacobian(), m);

    // ndof follows the hit-count convention 2 * nhits - 5 on every path that does not expand stubs.
    // A stub-expanded track carries both outer-tracker rechits of every stub in `hits` while the fit
    // used one position per stub, so there ndof counts the positions the fit used.
    int ndof = 2 * int(hits.size()) - 5;
    if (expandStubs_) {
      constexpr int maxHitsOnTrackForFullFit = 6;  // fallback only, see below
      ndof = 2 * std::min(nHits, maxHitsOnTrackForFullFit) - 5;
      // The fit kernel stamps the degrees of freedom of the positions it fitted (2N-5 for the
      // N <= maxHitsOnTrackForFullFit positions used) into the SoA; prefer it.
      const int ndofSoA = tsoa.view().tracks()[it].ndof();
      if (ndofSoA > 0)
        ndof = ndofSoA;
    }
    chi2 = chi2 * ndof;
    GlobalPoint vv = gp.position();
    math::XYZPoint pos(vv.x(), vv.y(), vv.z());
    GlobalVector pp = gp.momentum();
    math::XYZVector mom(pp.x(), pp.y(), pp.z());

    // A device fit that failed numerically can leave finite SoA parameters that still map to a
    // non-finite reco trajectory (through transformToPerigeePlane / the local->global transform).
    // Never emit such a track: drop it and free the edm index reserved above.
    if (not(std::isfinite(chi2) and std::isfinite(mom.x()) and std::isfinite(mom.y()) and std::isfinite(mom.z()))) {
      indToEdm[it] = pixelTrack::skippedTrack;
      --nt;
      continue;
    }

    auto track = std::make_unique<reco::Track>(chi2, ndof, pos, mom, gp.charge(), CurvilinearTrajectoryError(mo));

    // bad and edup not supported as fit not present or not reliable
    auto tkq = recoQuality[int(q)];
    track->setQuality(tkq);
    // loose,tight and HP are inclusive
    if (reco::TrackBase::highPurity == tkq) {
      track->setQuality(reco::TrackBase::tight);
      track->setQuality(reco::TrackBase::loose);
    } else if (reco::TrackBase::tight == tkq) {
      track->setQuality(reco::TrackBase::loose);
    }
    track->setQuality(tkq);
    if (setAlgorithmFromIteration_) {
      const auto iter = tsoa.view().tracks()[it].iteration();
      track->setAlgorithm(recoAlgo[std::min<uint32_t>(uint32_t(iter), pixelTrack::iterationSize)]);
    }
    if (fillTrackExtra_) {
      if (hits.empty()) {
        extras.emplace_back();
      } else {
        CurvilinearTrajectoryError const perigeeError(mo);
        ExtraState const inner = extraStateOnHit(FreeTrajectoryState(gp, perigeeError), *hits.front());
        // the outer helix bends in Bz averaged along the chord from the perigee to the last hit, closer to the
        // real bending of the long path than the field at the perigee
        GlobalPoint const &x0 = gp.position();
        GlobalPoint const x1 = hits.back()->globalPosition();
        GlobalPoint const xm(0.5f * (x0.x() + x1.x()), 0.5f * (x0.y() + x1.y()), 0.5f * (x0.z() + x1.z()));
        float const bz =
            0.25f * (gp.magneticFieldInTesla().z() + 2.f * idealField.inTesla(xm).z() + idealField.inTesla(x1).z());
        GlobalTrajectoryParameters const gpAveraged(
            x0, gp.momentum(), gp.charge(), &idealField, GlobalVector(0.f, 0.f, bz));
        ExtraState const outer = extraStateOnHit(FreeTrajectoryState(gpAveraged, perigeeError), *hits.back());
        nExtraStatesFailed += (not inner.ok) + (not outer.ok);
        extras.emplace_back(outer.position,
                            outer.momentum,
                            outer.ok,
                            inner.position,
                            inner.momentum,
                            inner.ok,
                            outer.covariance,
                            outer.detId,
                            inner.covariance,
                            inner.detId,
                            alongMomentum);
      }
    }
    // filter???
    tracks.emplace_back(track.release(), hits);
  }

#ifdef GPU_DEBUG
  std::cout << "processed " << nt << " good tuples " << tracks.size() << " out of " << indToEdm.size() << std::endl;
#endif

  // Diagnostic, printed only for events with tagged OT extras: "dropped" counts unresolvable tags
  // (no OT hitmap or out-of-range row). MessageLogger rate-limits per category, and produce() is
  // const (edm::global::EDProducer), so no local counter is possible anyway.
  if (verbose_ && nOTExtrasResolved + nOTExtrasDropped > 0)
    edm::LogInfo("PixelTrackProducerFromSoAAlpaka")
        << "tagged OT extras -> legacy hits: resolved=" << nOTExtrasResolved << " dropped=" << nOTExtrasDropped;
  if (verbose_ && (nExtraStatesFailed > 0 || nTracksReordered > 0 || nDuplicateHits > 0))
    edm::LogInfo("PixelTrackProducerFromSoAAlpaka")
        << "hits reordered across layers for " << nTracksReordered << " tracks, " << nDuplicateHits
        << " repeated measurements dropped; TrackExtra states: " << nExtraStatesFailed
        << " helices not reaching the plane of the first/last hit, replaced by the straight "
        << "line (" << tracks.size() << " tracks)";

  // store tracks
  storeTracks(iEvent, tracks, trackerTopology, std::move(extras));
  iEvent.put(std::move(indToEdmP));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PixelTrackProducerFromSoAAlpaka);
