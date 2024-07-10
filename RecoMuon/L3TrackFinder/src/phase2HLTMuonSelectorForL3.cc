/**  \class phase2HLTMuonSelectorForL3
 * 
 *
 *
 *   \author 
 */
#include "RecoMuon/L3TrackFinder/interface/phase2HLTMuonSelectorForL3.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"

#define SELECTOR_DEBUG

#ifdef SELECTOR_DEBUG
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

// Constructor
phase2HLTMuonSelectorForL3::phase2HLTMuonSelectorForL3(const edm::ParameterSet& iConfig)
    : l2MuTrajectoriesToken_(consumes<std::vector<Trajectory>>(iConfig.getParameter<edm::InputTag>("l2Muons"))),
      l2MuCollectionToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("l2MuonsUpdVtx"))),
      trajToL2MuAssociationMapToken_(
          consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("TrajToL2AssociationMap"))),
      l3TrackCollectionToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("l3Tracks"))),
      //l3IOTrackCollectionToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("l3IOTracks"))),
      //l3OITrackCollectionToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("l3OITracks"))),
      IOFirst_(iConfig.getParameter<bool>("IOFirst")),
      matchingDr_(iConfig.getParameter<double>("matchingDr")),
      applyL3Filters_(iConfig.getParameter<bool>("applyL3Filters")),
      maxNormalizedChi2_(iConfig.getParameter<double>("MaxNormalizedChi2")),
      maxPtDifference_(iConfig.getParameter<double>("MaxPtDifference")),
      minNhits_(iConfig.getParameter<int>("MinNhits")),
      minNHitsMuons_(iConfig.getParameter<int>("MinNmuonHits")) {
  produces<reco::TrackCollection>();
}

void phase2HLTMuonSelectorForL3::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l2Muons", edm::InputTag("hltL2MuonsFromL1TkMuon"));
  desc.add<edm::InputTag>("l2MuonsUpdVtx", edm::InputTag("hltL2MuonsFromL1TkMuon", "UpdatedAtVtx"));
  desc.add<edm::InputTag>("TrajToL2AssociationMap", edm::InputTag("hltL2MuonsFromL1TkMuon"));
  desc.add<edm::InputTag>("l3Tracks", edm::InputTag("hltIter2Phase2L3FromL1TkMuonMerged"));
  //desc.add<edm::InputTag>("l3IOTracks", edm::InputTag("hltIter2Phase2L3FromL1TkMuonMerged"));
  //desc.add<edm::InputTag>("l3OITracks", edm::InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"));
  desc.add<bool>("IOFirst", false);
  desc.add<double>("matchingDr", 0.02);
  desc.add<bool>("applyL3Filters", true);
  desc.add<int>("MinNhits", 1);
  desc.add<double>("MaxNormalizedChi2", 20.0);
  desc.add<int>("MinNmuonHits", 0);
  desc.add<double>("MaxPtDifference", 999.0);  //relative difference
  descriptions.add("phase2HLTMuonSelectorForL3", desc);
}

// Collection of L2 muons not already reconstructed as L3
void phase2HLTMuonSelectorForL3::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "Muon|RecoMuon|phase2HLTMuonSelectorForL3";

  // L2 Muons Trajectories
  auto const l2MuTrajectoriesH = iEvent.getHandle(l2MuTrajectoriesToken_);

  // L2 Muons collection
  auto const l2MuonsCollectionH = iEvent.getHandle(l2MuCollectionToken_);

  // Trajectory to L2 Muons association map
  // TrajTrackAssociationCollection == edm::AssociationMap<edm::OneToOne<std::vector<Trajectory>, reco::TrackCollection, unsigned short> >
  auto const trajToL2MuAssocH = iEvent.getHandle(trajToL2MuAssociationMapToken_);

  // L3 tracks (IO or OI)
  auto l3TracksCollectionH = iEvent.getHandle(l3TrackCollectionToken_);
  //    IOFirst_ ? iEvent.getHandle(l3IOTrackCollectionToken_) : iEvent.getHandle(l3OITrackCollectionToken_);

  // Output
  std::unique_ptr<reco::TrackCollection> result = std::make_unique<reco::TrackCollection>();

#ifdef SELECTOR_DEBUG
  if (IOFirst_) {
    LOG("Selector IO done first, looping on L2 muons");
  } else {
    LOG("Selector OI done first, looping on L2 muons");
  }
#endif
  // Loop on L2 Muons
  for (size_t l2MuIndex = 0; l2MuIndex != l2MuonsCollectionH->size(); ++l2MuIndex) {
    reco::TrackRef l2MuRef(l2MuonsCollectionH, l2MuIndex);
    bool reuseL2 = true;

    LOG("l2MuRef: " << l2MuRef.get());
    LOG("nHits: " << l2MuRef->recHitsSize() << " | innerOk: " << l2MuRef->innerOk()
                  << " | innerPosition: " << l2MuRef->innerPosition() << " | outerOk: " << l2MuRef->outerOk()
                  << " | outerPosition: " << l2MuRef->outerPosition());

    // Extract L1TkMu
    std::optional<l1t::TrackerMuonRef> l1TkMuRef = extractL1TkMu(l2MuRef, *l2MuTrajectoriesH, *trajToL2MuAssocH);

    if (l1TkMuRef) {
      // Loop on L3 tracks
      LOG("Looping on L3 tracks");
      for (size_t l3MuIndex = 0; l3MuIndex != l3TracksCollectionH->size(); ++l3MuIndex) {
        reco::TrackRef l3TrackRef(l3TracksCollectionH, l3MuIndex);

        float dR2 = deltaR2((*l1TkMuRef)->phEta(), (*l1TkMuRef)->phPhi(), l3TrackRef->eta(), l3TrackRef->phi());
        LOG("deltaR2: " << dR2);
        if (dR2 < matchingDr_ * matchingDr_) {
          reuseL2 = false;
          // Check L3 quality
          if (applyL3Filters_) {
            LOG("Checking L3 quality");
            reuseL2 = true;
            if (l3TrackRef->numberOfValidHits() < minNhits_)
              continue;  // cut on total number of hits
            if (l3TrackRef->normalizedChi2() > maxNormalizedChi2_)
              continue;  // cut on chi2
            if (l3TrackRef->hitPattern().numberOfValidMuonHits() < minNHitsMuons_)
              continue;  // cut on muon hits
            if (std::abs(l3TrackRef->pt() - (*l1TkMuRef)->phPt()) > maxPtDifference_ * l3TrackRef->pt())
              continue;  // cut on pt difference
            reuseL2 = false;
          }
        }
      }  // End loop on L3 Muons
    } else {
      LOG("Found L2 muon without an associated L1TkMu");
    }
    LOG("Reuse L2: " << reuseL2);
    if (reuseL2) {
      LOG("Found a L2 to be reused");
      result->push_back(*l2MuRef);
    }
  }  // End loop on L2 Muons

  LOG("Placing L2 Muons to be reused in the event");
  iEvent.put(std::move(result));
}

const std::optional<l1t::TrackerMuonRef> phase2HLTMuonSelectorForL3::extractL1TkMu(
    reco::TrackRef l2MuRef,
    const std::vector<Trajectory>& trajectories,
    const TrajTrackAssociationCollection& TrajToTrackAssociationMap) const {
  LOG("Extracting L1TkMu from L2 Trajectory");
  // Loop on association map to find the Trajectory matching with l2MuRef
  int n_associators = 1;
  for (TrajTrackAssociationCollection::const_iterator association = TrajToTrackAssociationMap.begin();
       association != TrajToTrackAssociationMap.end();
       ++association) {
    const Trajectory* traj = association->key.get();
    const reco::Track* assocL2Mu = association->val.get();

    LOG("assocL2Mu: " << assocL2Mu);
    LOG("nHits: " << assocL2Mu->recHitsSize() << " | innerOk: " << assocL2Mu->innerOk()
                  << " | innerPosition: " << assocL2Mu->innerPosition() << " | outerOk: " << assocL2Mu->outerOk()
                  << " | outerPosition: " << assocL2Mu->outerPosition());
    if (assocL2Mu->recHitsSize() == l2MuRef->recHitsSize() and
        (assocL2Mu->innerOk() and l2MuRef->innerOk() and assocL2Mu->innerPosition() == l2MuRef->innerPosition()) and
        (assocL2Mu->outerOk() and l2MuRef->outerOk() and assocL2Mu->outerPosition() == l2MuRef->outerPosition())) {
      // Extract seed from the trajectory
      edm::RefToBase<TrajectorySeed> seedRef = traj->seedRef();
      edm::Ref<L2MuonTrajectorySeedCollection> l2Seed = seedRef.castTo<edm::Ref<L2MuonTrajectorySeedCollection>>();
      // Extract L1TkMu from L2MuonTrajectorySeed
      LOG("Looped over " << n_associators << " TrajTrack associators. Returning L1TKMu");
      return std::make_optional<l1t::TrackerMuonRef>(l2Seed->l1TkMu());
    }
    ++n_associators;
  }
  LOG("Looped over " << n_associators << " TrajTrack associators. No associated L2 Muon found, returning no L1TkMu");
  return std::nullopt;
}

#undef SELECTOR_DEBUG

DEFINE_FWK_MODULE(phase2HLTMuonSelectorForL3);
