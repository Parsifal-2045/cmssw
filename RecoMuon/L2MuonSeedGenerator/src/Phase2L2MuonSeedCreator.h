#ifndef RecoMuon_L2MuonSeedCreator_Phase2L2MuonSeedCreator_H
#define RecoMuon_L2MuonSeedCreator_Phase2L2MuonSeedCreator_H

/** \class Phase2L2MuonSeedCreator
 *  
 *
 */

#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include <vector>

class RecHit;
class Plane;
class GeomDet;
class MagneticField;
class MuonTransientTrackingRecHit;

enum Type { barrel, overlap, endcap };

class Phase2L2MuonSeedCreator : public edm::stream::EDProducer<> {
public:
  typedef MuonTransientTrackingRecHit::MuonRecHitContainer SegmentContainer;

  /// Constructor
  explicit Phase2L2MuonSeedCreator(const edm::ParameterSet& pset);

  /// Destructor
  ~Phase2L2MuonSeedCreator() override = default;

  // Operations
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  // Tokens
  edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuCollToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentCollToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentCollToken_;

  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeometryToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<MuonDetLayerGeometry, MuonRecoGeometryRecord> muonLayersToken_;

  // Miminum and maximum pt momentum of a track
  double minMomentum_;
  double maxMomentum_;
  double dRCone_;
  double maxEtaBarrel_;   // barrel with |eta| < 0.7
  double maxEtaOverlap_;  // overlap with |eta| < 1.3, endcap after that

  // Handles
  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<MuonDetLayerGeometry> muonLayers_;
  edm::ESHandle<CSCGeometry> cscGeometry_;
  edm::ESHandle<DTGeometry> dtGeometry_;

  // Online sector 4 == offline sector 4 or 10, Online sector 10 == offline sector 10 or 14
  // Chambers are split due to material requirements, online doesn't have the split
  bool matchingDtIds(DTChamberId const& stubId, DTChamberId const& segId) {
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
};
#endif