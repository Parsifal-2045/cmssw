#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <unordered_set>

/** class MuonPixelTracksSelectorFromL1TkMuon
 *
 *   \author Luca Ferragina (INFN Bologna)
 */

class MuonPixelTracksSelectorFromL1TkMuon : public edm::stream::EDProducer<> {
public:
  MuonPixelTracksSelectorFromL1TkMuon(const edm::ParameterSet& iConfig);

  ~MuonPixelTracksSelectorFromL1TkMuon() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  // Tokens
  const edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuToken_;
  const edm::EDGetTokenT<reco::TrackCollection> tkToken_;
  // Track selection parameters
  const double l1TkMuMinPt_;
  const double l1TkMuMaxEta_;
  const double tkMinPt_;
  const bool cutOnSharedPt_;
  const double maxPtDiff_;
  const double maxDr_;
};

MuonPixelTracksSelectorFromL1TkMuon::MuonPixelTracksSelectorFromL1TkMuon(const edm::ParameterSet& iConfig)
    : l1TkMuToken_{consumes<l1t::TrackerMuonCollection>(
          iConfig.getParameter<edm::InputTag>("L1TkMuonInputCollection"))},
      tkToken_{consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackInputCollection"))},
      l1TkMuMinPt_{iConfig.getParameter<double>("L1TkMuMinPt")},
      l1TkMuMaxEta_{iConfig.getParameter<double>("L1TkMuMaxEta")},
      tkMinPt_{iConfig.getParameter<double>("trackMinPt")},
      cutOnSharedPt_{iConfig.getParameter<bool>("cutOnSharedPt")},
      maxPtDiff_{iConfig.getParameter<double>("maxNormalisedPtDifference")},
      maxDr_{iConfig.getParameter<double>("maxDr")} {
  produces<reco::TrackCollection>();
}

void MuonPixelTracksSelectorFromL1TkMuon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonPixelTracksSelectorFromL1TkMuon";
  // Get input collections
  auto const l1TkMuCollectionH = iEvent.getHandle(l1TkMuToken_);
  auto const tkCollectionH = iEvent.getHandle(tkToken_);
  if (!l1TkMuCollectionH.isValid() || !tkCollectionH.isValid() || l1TkMuCollectionH->empty() ||
      tkCollectionH->empty()) {
    LogDebug(metname) << "Input collection not valid or empty! Returning empty collection." << '\n';
    return;
  }
  // Output collection
  std::unordered_set<size_t> selectedTracksIndices;
  selectedTracksIndices.reserve(l1TkMuCollectionH->size() * 4);

  std::unique_ptr<reco::TrackCollection> outputTracks = std::make_unique<reco::TrackCollection>();
  outputTracks->reserve(l1TkMuCollectionH->size() * 4);

  std::cout << "MuonPixelTracksSelectorFromL1TkMuon: Starting selection of tracks from " << tkCollectionH->size()
            << " input tracks and " << l1TkMuCollectionH->size() << " L1Tk muons" << '\n';
  // Loop over L1TkMuons
  for (size_t l1TkMuIndex = 0; l1TkMuIndex != l1TkMuCollectionH->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuCollectionH, l1TkMuIndex);
    double l1TkMuEta, l1TkMuPhi, l1TkMuPt;
    auto trkPtr = l1TkMuRef->trkPtr();

    if (!trkPtr.isNull()) {
      // Prefer tracker-based coordinates
      l1TkMuEta = trkPtr->momentum().eta();
      l1TkMuPhi = trkPtr->momentum().phi();
      l1TkMuPt = trkPtr->momentum().perp();
    } else {
      // Fallback to standalone muon coordinates
      l1TkMuEta = l1TkMuRef->phEta();
      l1TkMuPhi = l1TkMuRef->phPhi();
      l1TkMuPt = l1TkMuRef->phPt();
      std::cout << "L1TkMuon has no tracker pointer, using muon system coordinates";
    }

    // Basic kinematic selection
    if (l1TkMuPt < l1TkMuMinPt_ || std::abs(l1TkMuEta) > l1TkMuMaxEta_) {
      continue;
    }
    // Loop over tracks
    for (size_t tkIndex = 0; tkIndex != tkCollectionH->size(); ++tkIndex) {
      // Skip already selected tracks
      if (selectedTracksIndices.contains(tkIndex)) {
        continue;  // Track already selected
      }

      reco::TrackRef tkRef(tkCollectionH, tkIndex);
      // Basic kinematic selection
      if (tkRef->pt() < tkMinPt_ || std::abs(tkRef->eta()) > l1TkMuMaxEta_) {
        continue;
      }
      // Check match in dR
      float dR2 = deltaR2(l1TkMuEta, l1TkMuPhi, tkRef->eta(), tkRef->phi());
      if (dR2 < maxDr_ * maxDr_) {
        // If requested, check match in shared pT
        if (cutOnSharedPt_) {
          float normPtDiff = std::abs(l1TkMuPt - tkRef->pt()) / l1TkMuPt;
          if (normPtDiff <= maxPtDiff_) {
            selectedTracksIndices.insert(tkIndex);
            outputTracks->push_back(*tkRef);
            std::cout << "MuonPixelTracksSelectorFromL1TkMuon: Track with index " << tkIndex
                      << " matched to L1Tk muon with dR = " << sqrt(dR2)
                      << " and normalised pT difference = " << normPtDiff << '\n';
          }
        } else {
          std::cout << "MuonPixelTracksSelectorFromL1TkMuon: Track with index " << tkIndex
                    << " matched to L1Tk muon with dR = " << sqrt(dR2) << '\n'
                    << "Normalised pT difference = " << std::abs(l1TkMuPt - tkRef->pt()) / l1TkMuPt << '\n';
          selectedTracksIndices.insert(tkIndex);
          outputTracks->push_back(*tkRef);
        }
      }
    }  // End loop over tracks
  }  // End loop over L1Tk muons

  std::cout << "MuonPixelTracksSelectorFromL1TkMuon: Selected " << outputTracks->size() << " tracks out of "
            << tkCollectionH->size() << " input tracks. Using " << l1TkMuCollectionH->size() << " L1 Tracker Muons\n";

  if (outputTracks->size() < l1TkMuCollectionH->size()) {
    std::cout << "MuonPixelTracksSelectorFromL1TkMuon: Warning! Fewer tracks selected (" << outputTracks->size()
              << ") than L1TkMu in event (" << l1TkMuCollectionH->size() << ")\n";
  }
  // Put output in the event
  iEvent.put(std::move(outputTracks));
}

void MuonPixelTracksSelectorFromL1TkMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("L1TkMuonInputCollection", edm::InputTag("l1tTkMuonsGmt"));
  // L1Tk muon selection parameters
  desc.add<double>("L1TkMuMinPt", 0.9);
  desc.add<double>("L1TkMuMaxEta", 2.5);

  // Track selection parameters
  desc.add<edm::InputTag>("TrackInputCollection", edm::InputTag("generalTracks"));
  desc.add<double>("trackMinPt", 0.9);
  desc.add<bool>("cutOnSharedPt", false);
  desc.add<double>("maxNormalisedPtDifference", 0.1);
  desc.add<double>("maxDr", 0.3);
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonPixelTracksSelectorFromL1TkMuon);
