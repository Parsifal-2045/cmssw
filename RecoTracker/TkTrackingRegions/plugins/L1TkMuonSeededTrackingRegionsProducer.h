#ifndef RecoTracker_TkTrackingRegions_L1TkMuonSeededTrackingRegionsProducer_h
#define RecoTracker_TkTrackingRegions_L1TkMuonSeededTrackingRegionsProducer_h

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/Record/interface/TrackerMultipleScatteringRecord.h"
#include "RecoTracker/TkMSParametrization/interface/MultipleScatteringParametrisationMaker.h"

#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"

/** class L1TkMuonSeededTrackingRegionsProducer
 *
 *   \author Luca Ferragina (INFN Bologna)
 *       based on RecoTracker/TkTrackingRegions/plugins/CandidateSeededTrackingRegionsProducer.h
 */

class L1TkMuonSeededTrackingRegionsProducer : public TrackingRegionProducer {
public:
  typedef enum { BEAM_SPOT_FIXED, BEAM_SPOT_SIGMA, VERTICES_FIXED, VERTICES_SIGMA } Mode;

  explicit L1TkMuonSeededTrackingRegionsProducer(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& iC)
      : BFieldToken_{iC.esConsumes()},
        l1TkMuToken_{iC.consumes(iConfig.getParameter<edm::InputTag>("L1TkMuonInputCollection"))},
        l1TkMuMinPt_{iConfig.getParameter<double>("L1TkMuMinPt")},
        l1TkMuMaxEta_{iConfig.getParameter<double>("L1TkMuMaxEta")},
        minPtBarrel_{iConfig.getParameter<double>("setMinPtBarrelTo")},
        minPtEndcap_{iConfig.getParameter<double>("setMinPtEndcapTo")} {
    edm::ParameterSet regPSet = iConfig.getParameter<edm::ParameterSet>("RegionPSet");

    // operation mode
    std::string modeString = regPSet.getParameter<std::string>("mode");
    if (modeString == "BeamSpotFixed")
      mode_ = BEAM_SPOT_FIXED;
    else if (modeString == "BeamSpotSigma")
      mode_ = BEAM_SPOT_SIGMA;
    else if (modeString == "VerticesFixed")
      mode_ = VERTICES_FIXED;
    else if (modeString == "VerticesSigma")
      mode_ = VERTICES_SIGMA;
    else
      edm::LogError("L1MuonSeededTrackingRegionsProducer") << "Unknown mode string: " << modeString;

    maxNVertices_ = 1;
    if (mode_ == VERTICES_FIXED || mode_ == VERTICES_SIGMA) {
      vertexToken_ = iC.consumes<reco::VertexCollection>(regPSet.getParameter<edm::InputTag>("vertexCollection"));
      maxNVertices_ = regPSet.getParameter<int>("maxNVertices");
    }

    beamSpotToken_ = iC.consumes<reco::BeamSpot>(regPSet.getParameter<edm::InputTag>("beamSpot"));
    maxNRegions_ = regPSet.getParameter<int>("maxNRegions");

    // RectangularEtaPhiTrackingRegion parameters:
    ptMin_ = regPSet.getParameter<double>("ptMin");
    originRadius_ = regPSet.getParameter<double>("originRadius");
    zErrorBeamSpot_ = regPSet.getParameter<double>("zErrorBeamSpot");
    ptRanges_ = regPSet.getParameter<std::vector<double>>("ptRanges");
    if (ptRanges_.size() < 2) {
      edm::LogError("L1TkMuonSeededTrackingRegionsProducer") << "Size of ptRanges cannot be less than 2" << std::endl;
    }
    deltaEtas_ = regPSet.getParameter<std::vector<double>>("deltaEtas");
    if (deltaEtas_.size() != ptRanges_.size() - 1) {
      edm::LogError("L1TkMuonSeededTrackingRegionsProducer")
          << "Size of deltaEtas does not match number of pt bins." << std::endl;
    }
    deltaPhis_ = regPSet.getParameter<std::vector<double>>("deltaPhis");
    if (deltaPhis_.size() != ptRanges_.size() - 1) {
      edm::LogError("L1TkMuonSeededTrackingRegionsProducer")
          << "Size of deltaPhis does not match number of pt bins." << std::endl;
    }

    precise_ = regPSet.getParameter<bool>("precise");
    if (precise_) {
      msmakerToken_ = iC.esConsumes();
    }

    whereToUseMeasurementTracker_ = RectangularEtaPhiTrackingRegion::stringToUseMeasurementTracker(
        regPSet.getParameter<std::string>("whereToUseMeasurementTracker"));
    if (whereToUseMeasurementTracker_ != RectangularEtaPhiTrackingRegion::UseMeasurementTracker::kNever) {
      measurementTrackerToken_ =
          iC.consumes<MeasurementTrackerEvent>(regPSet.getParameter<edm::InputTag>("measurementTrackerName"));
    }

    searchOpt_ = regPSet.getParameter<bool>("searchOpt");

    // mode-dependent z-halflength of tracking regions
    if (mode_ == VERTICES_SIGMA)
      nSigmaZVertex_ = regPSet.getParameter<double>("nSigmaZVertex");
    if (mode_ == VERTICES_FIXED)
      zErrorVertex_ = regPSet.getParameter<double>("zErrorVertex");
    nSigmaZBeamSpot_ = -1.;
    if (mode_ == BEAM_SPOT_SIGMA) {
      nSigmaZBeamSpot_ = regPSet.getParameter<double>("nSigmaZBeamSpot");
      if (nSigmaZBeamSpot_ < 0.)
        edm::LogError("L1TkMuonSeededTrackingRegionsProducer")
            << "nSigmaZBeamSpot must be positive for BeamSpotSigma mode!";
    }
  }

  ~L1TkMuonSeededTrackingRegionsProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("L1TkMuonInputCollection", edm::InputTag("l1tTkMuonsGmt"));
    // L1Tk muon selection parameters
    desc.add<double>("L1TkMuMinPt", 0.0);
    desc.add<double>("L1TkMuMaxEta", 2.5);
    desc.add<double>("setMinPtBarrelTo", 0.0);
    desc.add<double>("setMinPtEndcapTo", 0.0);

    // Tracking region parameters
    edm::ParameterSetDescription descRegion;
    descRegion.add<std::string>("mode", "BeamSpotSigma");
    descRegion.add<int>("maxNRegions", 10);
    descRegion.add<edm::InputTag>("beamSpot", edm::InputTag("hltOnlineBeamSpot"));
    descRegion.add<edm::InputTag>("vertexCollection", edm::InputTag("notUsed"));
    descRegion.add<int>("maxNVertices", 1);
    descRegion.add<double>("ptMin", 0.0);
    descRegion.add<double>("originRadius", 0.2);
    descRegion.add<double>("zErrorBeamSpot", 24.2);
    descRegion.add<std::vector<double>>("ptRanges", {0.0, 5.0, 10.0, 15.0, 20.0, 1e+64});
    descRegion.add<std::vector<double>>("deltaEtas", {0.015, 0.01, 0.01, 0.01, 0.01});
    descRegion.add<std::vector<double>>("deltaPhis", {0.05, 0.05, 0.02, 0.02, 0.01});
    descRegion.add<double>("nSigmaZVertex", 3.);
    descRegion.add<double>("zErrorVertex", 0.2);
    descRegion.add<double>("nSigmaZBeamSpot", 4.);
    descRegion.add<std::string>("whereToUseMeasurementTracker", "Never");
    descRegion.add<edm::InputTag>("measurementTrackerName", edm::InputTag(""));
    descRegion.add<bool>("precise", true);
    descRegion.add<bool>("searchOpt", false);
    desc.add<edm::ParameterSetDescription>("RegionPSet", descRegion);
    descriptions.addWithDefaultLabel(desc);
  }

  std::vector<std::unique_ptr<TrackingRegion>> regions(const edm::Event& iEvent,
                                                       const edm::EventSetup& iSetup) const override {
    std::vector<std::unique_ptr<TrackingRegion>> result;

    // pick up the candidate objects of interest
    edm::Handle<l1t::TrackerMuonCollection> L1TkMuCollH;
    iEvent.getByToken(l1TkMuToken_, L1TkMuCollH);
    if (L1TkMuCollH->size() == 0)
      return result;

    // always need the beam spot (as a fall back strategy for vertex modes)
    edm::Handle<reco::BeamSpot> bs;
    iEvent.getByToken(beamSpotToken_, bs);
    if (!bs.isValid())
      return result;

    // this is a default origin for all modes
    GlobalPoint default_origin(bs->x0(), bs->y0(), bs->z0());

    // vector of origin & halfLength pairs:
    std::vector<std::pair<GlobalPoint, float>> origins;

    // fill the origins and halfLengths depending on the mode
    if (mode_ == BEAM_SPOT_FIXED || mode_ == BEAM_SPOT_SIGMA) {
      origins.push_back(std::make_pair(default_origin,
                                       (mode_ == BEAM_SPOT_FIXED) ? zErrorBeamSpot_ : nSigmaZBeamSpot_ * bs->sigmaZ()));
    } else if (mode_ == VERTICES_FIXED || mode_ == VERTICES_SIGMA) {
      edm::Handle<reco::VertexCollection> vertices;
      iEvent.getByToken(vertexToken_, vertices);
      int nVertices = 0;
      for (reco::VertexCollection::const_iterator v = vertices->begin();
           v != vertices->end() && nVertices < maxNVertices_;
           ++v) {
        if (v->isFake() || !v->isValid())
          continue;

        origins.push_back(std::make_pair(GlobalPoint(v->x(), v->y(), v->z()),
                                         (mode_ == VERTICES_FIXED) ? zErrorVertex_ : nSigmaZVertex_ * v->zError()));
        ++nVertices;
      }
      // no-vertex fall-back case:
      if (origins.empty()) {
        origins.push_back(std::make_pair(default_origin,
                                         (nSigmaZBeamSpot_ > 0.) ? nSigmaZBeamSpot_ * bs->z0Error() : zErrorBeamSpot_));
      }
    }

    const MeasurementTrackerEvent* measurementTracker = nullptr;
    if (!measurementTrackerToken_.isUninitialized()) {
      measurementTracker = &iEvent.get(measurementTrackerToken_);
    }
    const auto& Bfield = iSetup.getData(BFieldToken_);
    const MultipleScatteringParametrisationMaker* msmaker = nullptr;
    if (precise_) {
      msmaker = &iSetup.getData(msmakerToken_);
    }

    // create tracking regions (maximum MaxNRegions of them) in directions of the
    // objects of interest (we expect that the collection was sorted in decreasing pt order)
    int nRegions = 0;
    for (size_t l1TkMuIndex = 0; l1TkMuIndex < L1TkMuCollH->size() && nRegions < maxNRegions_; ++l1TkMuIndex) {
      l1t::TrackerMuonRef l1TkMu(L1TkMuCollH, l1TkMuIndex);
      // Physical info of the L1Tk muon
      float pt = l1TkMu->phPt();
      const float eta = l1TkMu->phEta();
      if (pt < l1TkMuMinPt_ || std::abs(eta) > l1TkMuMaxEta_) {
        continue;
      }
      const bool barrel = eta < 1.2;
      if (barrel && pt < minPtBarrel_) {
        pt = minPtBarrel_;
      }
      if (!barrel && pt < minPtEndcap_) {
        pt = minPtEndcap_;
      }
      // set deltaEta and deltaPhi from L1Tk muon pt
      auto deltaEta = deltaEtas_.at(0);
      auto deltaPhi = deltaPhis_.at(0);
      if (pt < ptRanges_.back()) {
        auto lowEdge = std::upper_bound(ptRanges_.begin(), ptRanges_.end(), pt);
        deltaEta = deltaEtas_.at(lowEdge - ptRanges_.begin() - 1);
        deltaPhi = deltaPhis_.at(lowEdge - ptRanges_.begin() - 1);
      }
      // direction of the L1Tk muon
      const auto trk = l1TkMu->trkPtr();
      GlobalVector direction(trk->momentum().x(), trk->momentum().y(), trk->momentum().z());

      for (size_t j = 0; j < origins.size() && nRegions < maxNRegions_; ++j) {
        result.push_back(std::make_unique<RectangularEtaPhiTrackingRegion>(direction,
                                                                           origins[j].first,
                                                                           ptMin_,
                                                                           originRadius_,
                                                                           origins[j].second,
                                                                           deltaEta,
                                                                           deltaPhi,
                                                                           Bfield,
                                                                           msmaker,
                                                                           precise_,
                                                                           whereToUseMeasurementTracker_,
                                                                           measurementTracker,
                                                                           searchOpt_));
        ++nRegions;
      }
    }
    edm::LogInfo("L1TkMuonSeededTrackingRegionsProducer") << "produced " << nRegions << " regions";

    return result;
  }

private:
  Mode mode_;

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> BFieldToken_;
  edm::ESGetToken<MultipleScatteringParametrisationMaker, TrackerMultipleScatteringRecord> msmakerToken_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuToken_;
  edm::EDGetTokenT<MeasurementTrackerEvent> measurementTrackerToken_;

  std::vector<double> ptRanges_;
  std::vector<double> deltaEtas_;
  std::vector<double> deltaPhis_;

  const double l1TkMuMinPt_;
  const double l1TkMuMaxEta_;
  const double minPtBarrel_;
  const double minPtEndcap_;

  int maxNVertices_;
  int maxNRegions_;
  float ptMin_;
  float originRadius_;
  float zErrorBeamSpot_;
  float nSigmaZVertex_;
  float zErrorVertex_;
  float nSigmaZBeamSpot_;

  RectangularEtaPhiTrackingRegion::UseMeasurementTracker whereToUseMeasurementTracker_;

  bool precise_;
  bool searchOpt_;
};

#endif
