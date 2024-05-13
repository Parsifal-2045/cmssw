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

  /// Create a seed from set of segments
  L2MuonTrajectorySeed createSeed(Type muonType,
                                  CSCSegmentCollection cscSegments,
                                  DTRecSegment4DCollection dtSegments,
                                  l1t::TrackerMuonRef l1TkMuRef);

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
};
#endif