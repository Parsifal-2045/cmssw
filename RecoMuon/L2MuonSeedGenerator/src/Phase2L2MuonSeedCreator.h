#ifndef RecoMuon_L2MuonSeedCreator_Phase2L2MuonSeedCreator_H
#define RecoMuon_L2MuonSeedCreator_Phase2L2MuonSeedCreator_H

/** \class Phase2L2MuonSeedCreator
 *  
 *
 */

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include <vector>

class RecHit;
class Plane;
class GeomDet;
class MagneticField;
class MuonTransientTrackingRecHit;

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

  /// Cache Magnetic Field for current event
  void setBField(const MagneticField* theField) { BField = theField; };

  /// Create a seed from set of segments
  TrajectorySeed createSeed(int type, const SegmentContainer& seg, const l1t::TrackerMuon& l1TkMu, double const& dRCone);

private:
  /// Compute pt from parameters
  std::vector<double> getPt(const std::vector<double>& vParameters, double eta, double dPhi);

  /// Scale the dPhi from segment position
  double scaledPhi(double dphi, double t1);

  // Miminum and maximum pt momentum of a track
  double minMomentum_;
  double maxMomentum_;
  double defaultMomentum_;

  // Flag for internal debugging
  bool debug;

  // Error on pt estimate which prevents weighted average from blowing up ( spt --> 0 )
  double sysError;

  // Cache Magnetic Field for current event
  const MagneticField* BField;

  // seed parameters vectors
  std::vector<double> DT12;
  std::vector<double> DT13;
  std::vector<double> DT14;
  std::vector<double> DT23;
  std::vector<double> DT24;
  std::vector<double> DT34;

  // dphi scaling factors
  std::vector<double> DT12_1;
  std::vector<double> DT12_2;
  std::vector<double> DT13_1;
  std::vector<double> DT13_2;
  std::vector<double> DT14_1;
  std::vector<double> DT14_2;
  std::vector<double> DT23_1;
  std::vector<double> DT23_2;
  std::vector<double> DT24_1;
  std::vector<double> DT24_2;
  std::vector<double> DT34_1;
  std::vector<double> DT34_2;
};
#endif