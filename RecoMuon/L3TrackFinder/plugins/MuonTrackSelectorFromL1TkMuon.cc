/**  \class MuonTracksSelectorFromL1TkMuon
 *
 *  \brief Selector of reco::Track objects based on matching to L1 Tracker Muons
 *
 *  This module implements a selector for reco::Track objects based on their matching
 *  to L1 Tracker Muons (l1t::TrackerMuon). After a loose pre-selection based on kinematic
 *  and geometric criteria, tracks are ranked according to a combined quality metric that
 *  includes deltaR, deltaZ, and curvature compatibility with the L1TkMuon.
 *  The best matches are selected up to a configurable maximum number of tracks per L1TkMuon.
 *
 *   \author Luca Ferragina (INFN BO), 2026
 */

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

class MuonTracksSelectorFromL1TkMuon : public edm::stream::EDProducer<> {
public:
  MuonTracksSelectorFromL1TkMuon(const edm::ParameterSet& iConfig);

  ~MuonTracksSelectorFromL1TkMuon() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  // Tokens
  const edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuToken_;
  const edm::EDGetTokenT<reco::TrackCollection> tkToken_;
  // Track selection parameters
  const float l1TkMuMinPt_;
  const float tkMinPt_;
  const float tkMaxEta_;
  const float tkMaxPtForMatch_;
  const float maxDz_;
  const float maxDr_;
  const float maxChi2_;
  const float maxCombinedQuality_;
  const size_t nTracksToKeep_;
};

MuonTracksSelectorFromL1TkMuon::MuonTracksSelectorFromL1TkMuon(const edm::ParameterSet& iConfig)
    : l1TkMuToken_{consumes<l1t::TrackerMuonCollection>(
          iConfig.getParameter<edm::InputTag>("L1TkMuonInputCollection"))},
      tkToken_{consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackInputCollection"))},
      l1TkMuMinPt_{static_cast<float>(iConfig.getParameter<double>("L1TkMuMinPt"))},
      tkMinPt_{static_cast<float>(iConfig.getParameter<double>("trackMinPt"))},
      tkMaxEta_{static_cast<float>(iConfig.getParameter<double>("trackMaxEta"))},
      tkMaxPtForMatch_{static_cast<float>(iConfig.getParameter<double>("trackMaxPtForMatch"))},
      maxDz_{static_cast<float>(iConfig.getParameter<double>("maxDz"))},
      maxDr_{static_cast<float>(iConfig.getParameter<double>("maxDr"))},
      maxChi2_{static_cast<float>(iConfig.getParameter<double>("maxChi2"))},
      maxCombinedQuality_{static_cast<float>(iConfig.getParameter<double>("maxCombinedQuality"))},
      nTracksToKeep_{static_cast<size_t>(iConfig.getParameter<uint>("nTracksToKeep"))} {
  produces<reco::TrackCollection>();
}

void MuonTracksSelectorFromL1TkMuon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonTracksSelectorFromL1TkMuon";

  // Output collection
  std::unique_ptr<reco::TrackCollection> outputTracks = std::make_unique<reco::TrackCollection>();

  // Get input collections
  auto const l1TkMuCollectionH = iEvent.getHandle(l1TkMuToken_);
  auto const tkCollectionH = iEvent.getHandle(tkToken_);
  if (!l1TkMuCollectionH.isValid() || !tkCollectionH.isValid() || l1TkMuCollectionH->empty() ||
      tkCollectionH->empty()) {
    LogDebug(metname) << "Input collection not valid or empty! Returning empty collection." << '\n';
    iEvent.put(std::move(outputTracks));
    return;
  }

  // Store indices of matched tracks
  std::set<size_t> tracksToKeep;

  // Loop over L1TkMuons
  for (size_t l1TkMuIndex = 0; l1TkMuIndex != l1TkMuCollectionH->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuCollectionH, l1TkMuIndex);
    float l1TkMuEta, l1TkMuPhi, l1TkMuPt, l1TkMuZ0;
    auto trkPtr = l1TkMuRef->trkPtr();

    if (!trkPtr.isNull()) {
      // Prefer tracker-based coordinates
      l1TkMuEta = trkPtr->momentum().eta();
      l1TkMuPhi = trkPtr->momentum().phi();
      l1TkMuPt = trkPtr->momentum().perp();
      l1TkMuZ0 = trkPtr->z0();
    } else {
      // Fallback to standalone muon coordinates
      edm::LogWarning(metname) << "L1TkMuon with index " << l1TkMuIndex
                               << " has no tracker pointer, using muon system coordinates";
      l1TkMuEta = l1TkMuRef->phEta();
      l1TkMuPhi = l1TkMuRef->phPhi();
      l1TkMuPt = l1TkMuRef->phPt();
      l1TkMuZ0 = l1TkMuRef->phZ0();
    }

    // Basic kinematic selection
    if (l1TkMuPt < l1TkMuMinPt_) {
      continue;
    }

    // Pair of quality, track index for all candidate matches
    std::vector<std::pair<float, size_t>> potentialMatches;

    for (size_t tkIndex = 0; tkIndex != tkCollectionH->size(); ++tkIndex) {
      reco::TrackRef tkRef(tkCollectionH, tkIndex);

      // Loose pre-selection
      float dZ = std::abs(tkRef->dz() - l1TkMuZ0);
      if (dZ > 5.0)
        continue;

      float dR2 = deltaR2(l1TkMuEta, l1TkMuPhi, tkRef->eta(), tkRef->phi());
      if (dR2 > 0.4 * 0.4)
        continue;

      // Curvature compatibility
      float l1TkMuQoverPt = float(l1TkMuRef->phCharge()) / l1TkMuPt;
      float tkQoverPt = float(tkRef->charge()) / tkRef->pt();
      float tkQoverPtError = tkRef->ptError() / (tkRef->pt() * tkRef->pt());

      float curvDiff = l1TkMuQoverPt - tkQoverPt;
      float chi2Curv = (curvDiff * curvDiff) / (tkQoverPtError * tkQoverPtError);

      // Combine quality metrics
      float termDr = dR2 / (maxDr_ * maxDr_);        // (dR / maxDr)^2
      float termDz = (dZ * dZ) / (maxDz_ * maxDz_);  // (dZ / maxDz)^2
      float termChi2 = chi2Curv / maxChi2_;          // chi2 / maxChi2

      float combinedQuality = tkRef->pt() < tkMaxPtForMatch_ ? termDr + termDz + termChi2 : termDr + termDz;

      if (combinedQuality < maxCombinedQuality_) {
        potentialMatches.emplace_back(combinedQuality, tkIndex);
      }
    }  // End Track Loop

    // Sort by quality (lowest is best)
    std::partial_sort(potentialMatches.begin(),
                      potentialMatches.begin() + std::min(nTracksToKeep_, potentialMatches.size()),
                      potentialMatches.end());

    // Keep at most nTracksToKeep_ best matches per L1TkMuon
    for (size_t i = 0; i < std::min(potentialMatches.size(), nTracksToKeep_); ++i) {
      reco::TrackRef tkRef(tkCollectionH, potentialMatches[i].second);
      tracksToKeep.insert(potentialMatches[i].second);

      LogDebug(metname) << "L1Tk muon with index " << l1TkMuIndex << " matched to track with index "
                        << potentialMatches[i].second << " with combinedQuality = " << potentialMatches[i].first << '\n'
                        << "    -> dR = " << std::sqrt(deltaR2(l1TkMuEta, l1TkMuPhi, tkRef->eta(), tkRef->phi()))
                        << ", dZ =  " << std::abs(tkRef->dz() - l1TkMuZ0) << '\n'
                        << "   L1TkMu: pT = " << l1TkMuPt << ", eta = " << l1TkMuEta << ", phi = " << l1TkMuPhi
                        << ", z0 = " << l1TkMuZ0 << '\n'
                        << "Track:  pT = " << tkRef->pt() << ", eta = " << tkRef->eta() << ", phi = " << tkRef->phi()
                        << ", z0 = " << tkRef->dz();
    }
  }  // End L1TkMuon Loop

  outputTracks->reserve(tracksToKeep.size());
  for (size_t idx : tracksToKeep) {
    reco::TrackRef tkRef(tkCollectionH, idx);
    outputTracks->push_back(*tkRef);
  }

  iEvent.put(std::move(outputTracks));
}

void MuonTracksSelectorFromL1TkMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("L1TkMuonInputCollection", edm::InputTag("l1tTkMuonsGmt"))
      ->setComment("Input collection of L1TkMuons to match to");
  // L1Tk muon selection parameters
  desc.add<double>("L1TkMuMinPt", 0)->setComment("Minimum pT for L1TkMuons to be considered for matching");

  // Track selection parameters
  desc.add<edm::InputTag>("TrackInputCollection", edm::InputTag("generalTracks"))
      ->setComment("Input track collection to be filtered based on L1TkMuon matching");
  desc.add<double>("trackMinPt", 0.9)->setComment("Minimum pT for tracks to be considered for matching");
  desc.add<double>("trackMaxEta", 3.0)->setComment("Maximum |eta| for tracks to be considered for matching");
  desc.add<double>("trackMaxPtForMatch", 50.0)
      ->setComment("Maximum pT to apply curvature compatibility in matching (above this pT, only dR and dZ are used)");
  desc.add<double>("maxChi2", 9)
      ->setComment("Expected maximum chi2 for curvature compatibility between track and L1TkMuon");
  desc.add<double>("maxDr", 0.1)->setComment("Expected maximum deltaR for matching between track and L1TkMuon");
  desc.add<double>("maxDz", 1.0)->setComment("Expected maximum deltaZ for matching between track and L1TkMuon");
  desc.add<double>("maxCombinedQuality", 25.0)
      ->setComment("Maximum combined quality metric for matching (lower is better)");
  desc.add<uint32_t>("nTracksToKeep", 1)
      ->setComment("Maximum number of tracks to keep per L1TkMuon based on quality ranking");
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonTracksSelectorFromL1TkMuon);
