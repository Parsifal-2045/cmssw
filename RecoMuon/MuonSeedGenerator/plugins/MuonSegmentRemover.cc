/* \class MuonSegmentRemover
 *
 * \brief EDProducer that removes segments used by prompt muon reconstruction, producing
 * filtered segment collections for an iterative approach to displaced/non-prompt muon reconstruction
 *
 * \author Luca Ferragina (INFN BO), 2026
 */

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <set>
#include <memory>

class MuonSegmentRemover final : public edm::global::EDProducer<> {
public:
  MuonSegmentRemover(const edm::ParameterSet& iConfig);
  ~MuonSegmentRemover() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const override;

  // Helper function to check if a segment matches a hit
  bool segmentMatchesHit(const LocalPoint& segPos,
                         const LocalVector& segDir,
                         const LocalPoint& hitPos,
                         const float maxDistance) const {
    double deltaX = hitPos.x() - segPos.x();
    double deltaY = hitPos.y() - segPos.y();
    double distance2 = deltaX * deltaX + deltaY * deltaY;
    LogDebug("RecoMuon|MuonSeedGenerator|MuonSegmentRemover")
        << "Segment at (" << segPos.x() << ", " << segPos.y() << ") with direction (" << segDir.x() << ", "
        << segDir.y() << ") compared to hit at (" << hitPos.x() << ", " << hitPos.y()
        << ") has distance2: " << distance2;
    return distance2 < maxDistance * maxDistance;
  };

  // Input tokens
  const edm::EDGetTokenT<reco::TrackCollection> standaloneMuonsToken_;
  const edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsToken_;
  const edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentsToken_;

  // Quality filtering parameters
  const float maxMuonChi2_;
  const float maxSegmentChi2_;
  const int minNumberOfSegments_;
  const float segmentMatchDistance_;
};

MuonSegmentRemover::MuonSegmentRemover(const edm::ParameterSet& iConfig)
    : standaloneMuonsToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("standaloneMuons"))),
      cscSegmentsToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"))),
      dtSegmentsToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dtSegments"))),
      maxMuonChi2_(static_cast<float>(iConfig.getParameter<double>("maxMuonChi2"))),
      maxSegmentChi2_(static_cast<float>(iConfig.getParameter<double>("maxSegmentChi2"))),
      minNumberOfSegments_(iConfig.getParameter<int>("minNumberOfSegments")),
      segmentMatchDistance_(static_cast<float>(iConfig.getParameter<double>("segmentMatchDistance"))) {
  // Produce filtered segment collections
  produces<CSCSegmentCollection>();
  produces<DTRecSegment4DCollection>();
}

void MuonSegmentRemover::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {
  const std::string metname = "RecoMuon|MuonSeedGenerator|MuonSegmentRemover";

  // Get input collections
  edm::Handle<reco::TrackCollection> standaloneMuons;
  iEvent.getByToken(standaloneMuonsToken_, standaloneMuons);

  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(cscSegmentsToken_, cscSegments);

  edm::Handle<DTRecSegment4DCollection> dtSegments;
  iEvent.getByToken(dtSegmentsToken_, dtSegments);

  if (!standaloneMuons.isValid()) {
    edm::LogWarning("MuonSegmentRemover") << "Standalone muon collection not found, skipping segment removal";
    // If no standalone muons, just copy input segments to output
    auto filteredCSCSegments = std::make_unique<CSCSegmentCollection>(*cscSegments);
    auto filteredDTSegments = std::make_unique<DTRecSegment4DCollection>(*dtSegments);
    iEvent.put(std::move(filteredCSCSegments));
    iEvent.put(std::move(filteredDTSegments));
    return;
  }

  // Sets to track which segments to remove (using chamber ID and index within chamber)
  std::set<std::pair<uint32_t, int>> usedCSCSegments;
  std::set<std::pair<uint32_t, int>> usedDTSegments;

  // Loop over standalone muon tracks to identify used segments
  for (const auto& muonTrack : *standaloneMuons) {
    // Apply quality cuts
    LogDebug(metname) << "Processing muon with pT: " << muonTrack.pt() << " GeV, eta: " << muonTrack.eta()
                      << ", phi: " << muonTrack.phi() << " reduced chi2: " << muonTrack.normalizedChi2();
    if (muonTrack.normalizedChi2() > maxMuonChi2_) {
      continue;
    }

    // Count number of muon segments (DT + CSC)
    int nMuonSegments = 0;
    for (auto const& hit : muonTrack.recHits()) {
      if (!hit->isValid())
        continue;
      DetId detId = hit->geographicalId();
      LogDebug(metname) << "Found hit/segment with DetID: " << detId.rawId() << " in (det/subdet): " << detId.det()
                        << "/" << detId.subdetId() << " with " << hit->recHits().size() << " hits";
      if (detId.det() != DetId::Muon) {
        edm::LogWarning(metname) << "Found non-muon hit: " << detId.rawId() << " in det " << detId.det() << " skipping";
        continue;
      }
      if (detId.subdetId() == MuonSubdetId::DT || detId.subdetId() == MuonSubdetId::CSC) {
        if (!hit->recHits().empty()) {  //&& hit->recHits().size() > 1) {
          LogDebug(metname) << "Found a DT/CSC muon segment with " << hit->recHits().size() << " hits";
          nMuonSegments++;
        }
      }
    }
    if (nMuonSegments == 0) {
      std::cout << "-------------WARNING: Muon has no valid DT/CSC segments-------------\n";
    }
    LogDebug(metname) << "Standalone Muon has " << nMuonSegments << " DT+CSC segments";

    // Apply minimum segment requirement
    if (nMuonSegments < minNumberOfSegments_)
      continue;

    // Loop over RecHits and identify matching segments
    for (auto const& hit : muonTrack.recHits()) {
      if (!hit->isValid())
        continue;

      DetId detId = hit->geographicalId();
      if (detId.det() != DetId::Muon)
        continue;

      int subdet = detId.subdetId();
      LocalPoint hitPos = hit->localPosition();

      // Process DT segments
      if (subdet == MuonSubdetId::DT) {
        DTChamberId chamberId(detId.rawId());
        DTRecSegment4DCollection::range range = dtSegments->get(chamberId);
        int segmentIndex = 0;
        LogDebug(metname) << "Processing DT segments in chamber " << chamberId << " for hit at position (" << hitPos.x()
                          << ", " << hitPos.y() << ")";

        for (auto segment = range.first; segment != range.second; ++segment, ++segmentIndex) {
          // Check chi2 cut
          LogDebug(metname) << "Checking DT segment " << segmentIndex << " with chi2: " << segment->chi2();
          if (segment->chi2() > maxSegmentChi2_)
            continue;

          // Match segment to hit
          if (segmentMatchesHit(segment->localPosition(), segment->localDirection(), hitPos, segmentMatchDistance_)) {
            LogDebug(metname) << "Matched DT segment reduced chi2: " << segment->chi2() / segment->degreesOfFreedom();
            usedDTSegments.insert(std::make_pair(chamberId.rawId(), segmentIndex));
          }
        }
      }

      // Process CSC segments
      else if (subdet == MuonSubdetId::CSC) {
        CSCDetId cscId(detId.rawId());
        CSCSegmentCollection::range range = cscSegments->get(cscId);
        int segmentIndex = 0;
        LogDebug(metname) << "Processing CSC segments in chamber " << cscId << " for hit at position (" << hitPos.x()
                          << ", " << hitPos.y() << ")";

        for (auto segment = range.first; segment != range.second; ++segment, ++segmentIndex) {
          // Check chi2 cut
          LogDebug(metname) << "Checking CSC segment " << segmentIndex << " with chi2: " << segment->chi2();
          if (segment->chi2() > maxSegmentChi2_)
            continue;

          // Match segment to hit
          if (segmentMatchesHit(segment->localPosition(), segment->localDirection(), hitPos, segmentMatchDistance_)) {
            LogDebug(metname) << "Matched CSC segment reduced chi2: " << segment->chi2() / segment->degreesOfFreedom();
            usedCSCSegments.insert(std::make_pair(cscId.rawId(), segmentIndex));
          }
        }
      }
    }  // end loop over RecHits
  }  // end loop over standalone muons

  // Create output collections with unused segments
  auto filteredCSCSegments = std::make_unique<CSCSegmentCollection>();
  auto filteredDTSegments = std::make_unique<DTRecSegment4DCollection>();

  // Filter CSC segments - collect by chamber
  std::map<CSCDetId, std::vector<CSCSegment>> cscSegmentsByDetId;
  for (auto chamberIt = cscSegments->id_begin(); chamberIt != cscSegments->id_end(); ++chamberIt) {
    CSCDetId cscId(*chamberIt);
    CSCSegmentCollection::range range = cscSegments->get(cscId);
    int segmentIndex = 0;

    for (auto segment = range.first; segment != range.second; ++segment, ++segmentIndex) {
      // Keep segment if it's not in the used set
      if (usedCSCSegments.find(std::make_pair(cscId.rawId(), segmentIndex)) == usedCSCSegments.end()) {
        cscSegmentsByDetId[cscId].push_back(*segment);
      }
    }
  }

  // Insert collected segments per chamber
  for (const auto& [detId, segments] : cscSegmentsByDetId) {
    if (!segments.empty()) {
      filteredCSCSegments->put(detId, segments.begin(), segments.end());
    }
  }

  // Filter DT segments - collect by chamber
  std::map<DTChamberId, std::vector<DTRecSegment4D>> dtSegmentsByDetId;
  for (auto chamberIt = dtSegments->id_begin(); chamberIt != dtSegments->id_end(); ++chamberIt) {
    DTChamberId chamberId(*chamberIt);
    DTRecSegment4DCollection::range range = dtSegments->get(chamberId);
    int segmentIndex = 0;

    for (auto segment = range.first; segment != range.second; ++segment, ++segmentIndex) {
      // Keep segment if it's not in the used set
      if (usedDTSegments.find(std::make_pair(chamberId.rawId(), segmentIndex)) == usedDTSegments.end()) {
        dtSegmentsByDetId[chamberId].push_back(*segment);
      }
    }
  }

  // Insert collected segments per chamber
  for (const auto& [detId, segments] : dtSegmentsByDetId) {
    if (!segments.empty()) {
      filteredDTSegments->put(detId, segments.begin(), segments.end());
    }
  }

  assert(cscSegments->size() == usedCSCSegments.size() + filteredCSCSegments->size());
  assert(dtSegments->size() == usedDTSegments.size() + filteredDTSegments->size());

  LogDebug(metname) << "CSC segments: " << cscSegments->size() << " total, " << usedCSCSegments.size() << " removed, "
                    << filteredCSCSegments->size() << " remaining";
  LogDebug(metname) << "DT segments: " << dtSegments->size() << " total, " << usedDTSegments.size() << " removed, "
                    << filteredDTSegments->size() << " remaining";

  iEvent.put(std::move(filteredCSCSegments));
  iEvent.put(std::move(filteredDTSegments));
}

void MuonSegmentRemover::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("standaloneMuons", edm::InputTag("standAloneMuons"));
  desc.add<edm::InputTag>("cscSegments", edm::InputTag("cscSegments"));
  desc.add<edm::InputTag>("dtSegments", edm::InputTag("dt4DSegments"));
  desc.add<double>("maxMuonChi2", 10.0);
  desc.add<double>("maxSegmentChi2", 25.0);
  desc.add<int>("minNumberOfSegments", 1);
  desc.add<double>("segmentMatchDistance", 0.1);  // 1mm

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonSegmentRemover);
