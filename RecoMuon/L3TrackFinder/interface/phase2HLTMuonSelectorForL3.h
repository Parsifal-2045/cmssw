#ifndef RecoMuon_L3TrackFinder_phase2HLTMuonSelectorForL3_H
#define RecoMuon_L3TrackFinder_phase2HLTMuonSelectorForL3_H

/**  \class phase2HLTMuonSelectorForL3
 * 
 *
 *   \author 
 */

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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"


#include <optional>

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}  // End namespace edm

class phase2HLTMuonSelectorForL3 : public edm::stream::EDProducer<> {
public:
  // Constructor
  phase2HLTMuonSelectorForL3(const edm::ParameterSet&);

  // Destructor
  ~phase2HLTMuonSelectorForL3() override = default;

  // Default values
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // Select muons
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  const edm::EDGetTokenT<std::vector<Trajectory>> l2MuTrajectoriesToken_;
  const edm::EDGetTokenT<reco::TrackCollection> l2MuCollectionToken_;
  const edm::EDGetTokenT<TrajTrackAssociationCollection> trajToL2MuAssociationMapToken_;
  const edm::EDGetTokenT<reco::TrackCollection> l3TrackCollectionToken_;
  
  //const edm::EDGetTokenT<reco::TrackCollection> l3IOTrackCollectionToken_;
  //const edm::EDGetTokenT<reco::TrackCollection> l3OITrackCollectionToken_;
  const bool IOFirst_;
  const double matchingDr_;
  const bool applyL3Filters_;
  const double maxNormalizedChi2_, maxPtDifference_;
  const int minNhits_, minNHitsMuons_;

  const std::optional<l1t::TrackerMuonRef> extractL1TkMu(
      reco::TrackRef l2MuRef,
      const std::vector<Trajectory>& trajectories,
      const TrajTrackAssociationCollection& TrajToTrackAssociationMap) const;
};

#endif
