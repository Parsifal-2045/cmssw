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
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/Common/interface/RefCoreGet.h"
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

/**
 * This class creates "legacy" reco::Track
 * objects from the output of SoA CA.
 */

// #define GPU_DEBUG
// per-run data: the maps of the OT modules used by the extension, and per tracker detUnit the words the track
// loop would otherwise compute from the topology for every hit
struct DetIdMaps {
  DetIdMaps() : detIdToOTModuleId_(), detIdIsUsedOTModule_() {}

  // map from the detId of OT modules to the moduleId among the used OT modules
  // (starting from 0 for first module of first OT layer)
  std::map<uint32_t, uint32_t> detIdToOTModuleId_;
  // map from detId to bool if used as OT extension
  std::map<uint32_t, bool> detIdIsUsedOTModule_;
  // per tracker detUnit, by GeomDet::index(): the HitPattern word of a valid hit and the layer key of the
  // crossing-order sort, with the detId they were computed for
  std::vector<uint32_t> detUnitDetId_;
  std::vector<uint16_t> validHitPattern_;
  std::vector<uint32_t> layerKey_;
};

namespace {
  // subdetector, side and layer of a module in one word: two hits are on the same layer iff their words are equal
  uint32_t layerKeyOf(DetId const &id, TrackerTopology const &topology) {
    return (uint32_t(id.subdetId()) << 24) | (uint32_t(topology.side(id)) << 16) | topology.layer(id);
  }

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
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyTokenRun_;

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
      trackerTopologyTokenRun_(esConsumes<edm::Transition::BeginRun>()),
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

  // per detUnit, the HitPattern word of a valid hit (appending a hit by its word skips the topology lookups of
  // HitPattern::encode) and the layer key
  auto const &geometry = iSetup.getData(trackerGeometryTokenRun_);
  auto const &topology = iSetup.getData(trackerTopologyTokenRun_);
  // GeomDet::index() is the position in detUnits(); the range is taken from the indices themselves
  int maxIndex = -1;
  for (auto const *detUnit : geometry.detUnits())
    maxIndex = std::max(maxIndex, detUnit->index());
  detIdMaps->detUnitDetId_.assign(maxIndex + 1, 0);
  detIdMaps->validHitPattern_.assign(maxIndex + 1, 0);
  detIdMaps->layerKey_.assign(maxIndex + 1, 0);
  for (auto const *detUnit : geometry.detUnits()) {
    DetId const id = detUnit->geographicalId();
    // detUnits the topology cannot describe stay unknown to the cache (their hits take the topology path)
    if (detUnit->index() < 0 or id.det() != DetId::Tracker or id.subdetId() < 1 or id.subdetId() > 6)
      continue;
    reco::HitPattern pattern;
    pattern.appendHit(id, TrackingRecHit::valid, topology);
    detIdMaps->detUnitDetId_[detUnit->index()] = id.rawId();
    detIdMaps->validHitPattern_[detUnit->index()] = pattern.getHitPattern(reco::HitPattern::TRACK_HITS, 0);
    detIdMaps->layerKey_[detUnit->index()] = layerKeyOf(id, topology);
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

  // legacy tracks, built in place in their output collection, and their hits (all tracks in a row): the hits are
  // cloned into their collection when the products are stored
  auto outputTracks = std::make_unique<reco::TrackCollection>();
  std::vector<TrackingRecHit const *> outputHits;

  // get trackerTopology
  auto const &trackerTopology = iSetup.getData(trackerTopologyToken_);

  // get the maps for the detId of the OT modules
  auto const &runData = *runCache(iEvent.getRun().index());
  auto const &detIdIsUsedOTModule = runData.detIdIsUsedOTModule_;
  auto const &detIdToOTModuleId = runData.detIdToOTModuleId_;
  auto const &detUnitDetId = runData.detUnitDetId_;
  auto const &validHitPattern = runData.validHitPattern_;
  auto const &cachedLayerKey = runData.layerKey_;
  // index of the detUnit of a hit in the per-detUnit caches, or -1 if the caches do not know it
  auto cacheIndex = [&](TrackingRecHit const &hit) -> int {
    auto const *det = hit.det();
    if (det == nullptr)
      return -1;
    auto const index = det->index();
    return (index >= 0 and size_t(index) < detUnitDetId.size() and detUnitDetId[index] == hit.geographicalId().rawId())
               ? index
               : -1;
  };
  // append a hit (with its cache index) to the HitPattern of a track: a valid hit by the cached word of its detUnit,
  // anything else (or a detUnit the cache does not know) through the topology, with the same result
  auto appendHitPattern = [&](reco::Track &track, TrackingRecHit const &hit, int index) {
    if (hit.getType() == TrackingRecHit::valid and index >= 0)
      return track.appendHitPattern(validHitPattern[index], TrackingRecHit::valid);
    return track.appendHitPattern(hit, trackerTopology);
  };

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
  // On the stub path each OT SoA row stores origRecHitIdx, the flat index of its legacy RecHit in the
  // Phase2TrackerRecHit1DCollectionNew (assigned while iterating the DetSets in legacy order): an OT row then resolves
  // to its legacy RecHit directly when a track uses it (otHitOfRow below), and the hitmap holds the pixel hits only.
  bool const otRowsResolvedDirectly = useOTExtension_ && expandStubs_ && otRecHitsSoAHost != nullptr;
  std::vector<TrackingRecHit const *> hitmap(otRowsResolvedDirectly ? nPixelHits : nTotalHits, nullptr);

  // loop over pixel RecHits to fill the hitmap, module by module; the clusters are read by index from the data of
  // their DetSetVector (that of the first hit, fetched once), rather than through an edm::Ref per hit, which
  // resolves to the same element but takes the filling lock of the DetSetVector (an atomic compare-exchange and a
  // store) on every dereference
  SiPixelCluster const *pixelClusters = nullptr;
  edm::ProductID pixelClustersID;
  for (auto const &detSet : pixelRecHitsDSV) {
    if (detSet.empty())
      continue;
    // hit identifier = (hit offset of the module) + (hit index in this module)
    auto const moduleStart = pixelHitsModuleStart[detSet.begin()->det()->index()];
    for (auto const &pixelHit : detSet) {
      auto const &clus = pixelHit.omniClusterRef();
      assert(clus.isPixel());
      if (pixelClusters == nullptr) {
        auto const clusterRef = clus.cluster_pixel();
        assert(clusterRef.isNonnull());
        pixelClusters =
            edm::getProductWithCoreFromRef<SiPixelClusterCollectionNew>(clusterRef.refCore(), &iEvent.productGetter())
                ->data()
                .data();
        pixelClustersID = clus.id();
      }
      auto const &cluster = clus.id() == pixelClustersID ? pixelClusters[clus.index()] : clus.pixelCluster();
      auto const idx = moduleStart + cluster.originalId();
      assert(nullptr == hitmap[idx]);
      hitmap[idx] = &pixelHit;
    }
  }

  // legacy RecHit of an OT SoA row on the stub path; rows past the SoA resolve to nullptr, as unfilled map entries
  uint32_t const nOTRowsSoA = otRowsResolvedDirectly ? otRecHitsSoAView.metadata().size() : 0;
  auto otHitOfRow = [&](uint32_t row) -> TrackingRecHit const * {
    if (row >= nOTRowsSoA)
      return nullptr;
    uint32_t const flatIdx = otRecHitsSoAView[row].origRecHitIdx();
    assert(flatIdx < otRecHitsDSV->data().size());
    return &otRecHitsDSV->data()[flatIdx];
  };

  // without stub expansion, fill the hitmap also with the OT RecHits used by the extension
  if (useOTExtension_ and not otRowsResolvedDirectly) {
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
  std::vector<uint32_t> layerKeys;
  std::vector<float> planeDistances;
  auto moduleDistance = [&bs](TrackingRecHit const *hit) { return (hit->det()->position() - bs).mag(); };
  auto planeDistance = [&bs](TrackingRecHit const *hit) {
    auto const &surface = hit->det()->surface();
    return std::abs(surface.normalVector().dot(surface.position() - bs));
  };
  // the cache index of each hit is returned in cacheIndices, in the final order, for the HitPattern of the track
  auto sortInCrossingOrder = [&](std::vector<const TrackingRecHit *> &trackHits, std::vector<int> &cacheIndices) {
    if (trackHits.size() < 2) {
      cacheIndices.resize(trackHits.size());
      for (size_t i = 0; i < trackHits.size(); ++i)
        cacheIndices[i] = cacheIndex(*trackHits[i]);
      return;
    }
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
    auto const nHits = trackHits.size();
    cacheIndices.resize(nHits);
    layerKeys.resize(nHits);
    for (size_t i = 0; i < nHits; ++i) {
      cacheIndices[i] = cacheIndex(*trackHits[i]);
      layerKeys[i] = cacheIndices[i] >= 0 ? cachedLayerKey[cacheIndices[i]]
                                          : layerKeyOf(trackHits[i]->geographicalId(), trackerTopology);
    }
    // compacting in place: the kept hits of a run are written at `kept`, never past the hit being read
    size_t kept = 0;
    for (size_t first = 0; first < nHits;) {
      size_t last = first + 1;
      while (last < nHits and layerKeys[last] == layerKeys[first])
        ++last;
      if (last - first > 1) {
        // stable insertion sort by the plane distance: a run holds a few hits
        planeDistances.resize(last - first);
        for (size_t i = first; i < last; ++i)
          planeDistances[i - first] = planeDistance(trackHits[i]);
        for (size_t i = first + 1; i < last; ++i) {
          auto const *hit = trackHits[i];
          int const index = cacheIndices[i];
          float const distance = planeDistances[i - first];
          size_t j = i;
          for (; j > first and distance < planeDistances[j - 1 - first]; --j) {
            trackHits[j] = trackHits[j - 1];
            cacheIndices[j] = cacheIndices[j - 1];
            planeDistances[j - first] = planeDistances[j - 1 - first];
          }
          trackHits[j] = hit;
          cacheIndices[j] = index;
          planeDistances[j - first] = distance;
        }
      }
      size_t const runBegin = kept;
      for (size_t i = first; i < last; ++i) {
        auto const *hit = trackHits[i];
        if (std::any_of(trackHits.begin() + runBegin, trackHits.begin() + kept, [&](TrackingRecHit const *k) {
              return k == hit or
                     (k->geographicalId() == hit->geographicalId() and k->sharesInput(hit, TrackingRecHit::all));
            })) {
          ++nDuplicateHits;
          continue;
        }
        cacheIndices[kept] = cacheIndices[i];
        trackHits[kept++] = hit;
      }
      first = last;
    }
    trackHits.resize(kept);
    cacheIndices.resize(kept);
  };

  std::vector<const TrackingRecHit *> hits;
  hits.reserve(5);  //TODO move to a configurable parameter?
  // parallel to hits: the index of each hit in the per-detUnit caches, filled by sortInCrossingOrder
  std::vector<int> hitCacheIndices;

  auto const &tsoa = iEvent.get(trackSoAToken_);
  auto const quality = tsoa.view().tracks().quality();
  auto const hitOffs = tsoa.view().tracks().hitOffsets();
  // Plain column accessor for pt, used by the sort comparator below: tsoa.view().tracks()[i].pt()
  // would build a full element proxy per comparison, and that proxy's constructor builds the
  // Eigen::Map members of the layout's two Eigen columns.
  auto const trackPt = tsoa.view().tracks().pt();
  auto const hitIdxs = tsoa.view().trackHits().id();
  auto nTracks = tsoa.view().tracks().nTracks();

  outputTracks->reserve(nTracks);
  // stub expansion at most doubles the hits of a track
  outputHits.reserve(nTracks > 0 ? size_t(hitOffs[nTracks - 1]) * (expandStubs_ ? 2 : 1) : 0);

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
  // its OT SoA row (otHitOfRow), available only on the expandStubs OT
  // path. Where it is unavailable the tagged extra is dropped, never crashing.
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
          if (!(otRowsResolvedDirectly && (nPixelHits + o) < nTotalHits))
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
          if (!(otRowsResolvedDirectly && (nPixelHits + o) < nTotalHits))
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
        // a stub's lower/upper sensor hit: otHitOfRow(otSoARow)). Unresolvable tagged ids
        // are dropped (counted as removed above), keeping the hits vector correctly sized.
        const uint32_t o = caOTHitTag::otIdx(hitIdx);
        if (otRowsResolvedDirectly && (nPixelHits + o) < nTotalHits) {
          hits[hitOutputIdx++] = otHitOfRow(o);
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

          hits[hitOutputIdx++] = otHitOfRow(lowerHitIdx);

          // Add outer sensor hit only if not PHitOnly (PHitOnly stubs have invalid upperHitIdx)
          if (isStub(stubsSoAView, stubIdx)) {
            uint32_t upperHitIdx = stubsSoAView[stubIdx].upperHitIdx();
            hits[hitOutputIdx++] = otHitOfRow(upperHitIdx);
          }
        } else {
          hits[hitOutputIdx++] = hitmap[hitIdx];
        }
      }
      // else: removed hits are skipped
    }
    sortInCrossingOrder(hits, hitCacheIndices);

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
    float chi2 = tsoa.view().tracks().chi2()[it];
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
      const int ndofSoA = tsoa.view().tracks().ndof()[it];
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

    auto &track = outputTracks->emplace_back(chi2, ndof, pos, mom, gp.charge(), CurvilinearTrajectoryError(mo));

    // bad and edup not supported as fit not present or not reliable
    auto tkq = recoQuality[int(q)];
    track.setQuality(tkq);
    // loose,tight and HP are inclusive
    if (reco::TrackBase::highPurity == tkq) {
      track.setQuality(reco::TrackBase::tight);
      track.setQuality(reco::TrackBase::loose);
    } else if (reco::TrackBase::tight == tkq) {
      track.setQuality(reco::TrackBase::loose);
    }
    track.setQuality(tkq);
    if (setAlgorithmFromIteration_) {
      const auto iter = tsoa.view().tracks().iteration()[it];
      track.setAlgorithm(recoAlgo[std::min<uint32_t>(uint32_t(iter), pixelTrack::iterationSize)]);
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
    for (size_t k = 0; k < hits.size(); ++k) {
      appendHitPattern(track, *hits[k], hitCacheIndices[k]);
      outputHits.push_back(hits[k]);
    }
  }

#ifdef GPU_DEBUG
  std::cout << "processed " << nt << " good tuples " << outputTracks->size() << " out of " << indToEdm.size()
            << std::endl;
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
        << "line (" << outputTracks->size() << " tracks)";

  // store the hits, the TrackExtras and the tracks, which refer to the other two
  auto const nStored = outputTracks->size();
  auto recHits = std::make_unique<TrackingRecHitCollection>();
  recHits->reserve(outputHits.size());
  for (auto const *hit : outputHits)
    recHits->push_back(hit->clone());  // need to clone (at least if from SoA)
  edm::OrphanHandle<TrackingRecHitCollection> const recHitsHandle = iEvent.put(std::move(recHits));
  edm::RefProd<TrackingRecHitCollection> const recHitsRef(recHitsHandle);

  auto trackExtras = std::make_unique<reco::TrackExtraCollection>(std::move(extras));
  trackExtras->resize(nStored);
  LocalTrajectoryParameters const noTrajParams(AlgebraicVector5(0, 0, 0, 0, 0), 1.);
  for (unsigned int k = 0, firstHit = 0; k < nStored; ++k) {
    unsigned int const nValidHits = (*outputTracks)[k].numberOfValidHits();
    auto &extra = (*trackExtras)[k];
    extra.setHits(recHitsRef, firstHit, nValidHits);
    firstHit += nValidHits;
    extra.setTrajParams(reco::TrackExtra::TrajParams(nValidHits, noTrajParams),
                        reco::TrackExtra::Chi2sFive(nValidHits, 0));
  }
  edm::OrphanHandle<reco::TrackExtraCollection> const extrasHandle = iEvent.put(std::move(trackExtras));
  for (unsigned int k = 0; k < nStored; ++k)
    (*outputTracks)[k].setExtra(reco::TrackExtraRef(extrasHandle, k));
  iEvent.put(std::move(outputTracks));
  iEvent.put(std::move(indToEdmP));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PixelTrackProducerFromSoAAlpaka);
