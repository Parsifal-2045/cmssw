/**  \class MuonSeedsSelectorFromL1TkMuon
 *
 *  \brief Selector of TrajectorySeed objects based on matching to L1 Tracker Muons
 *
 *  This module implements a selector for TrajectorySeed objects based on their matching
 *  to L1 Tracker Muons (l1t::TrackerMuon). The selection strategy is inspired by the
 *  MuonIdProducer arbitration logic:
 *
 *  1. Loose geometric pre-selection (dR cone)
 *  2. Quality evaluation based on seed intrinsic properties (hit counting, layer quality)
 *  3. Hit-by-hit consistency check with L1TkMuon tracker stubs (OT hits only)
 *  4. Arbitration: keep only the N best seeds per L1TkMuon based on a combined quality metric
 *
 *  Adapted from MuonTracksSelectorFromL1TkMuon
 *  
 *  \author Luca Ferragina (INFN BO), 2026
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
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include <set>
#include <vector>
#include <algorithm>

class MuonSeedsSelectorFromL1TkMuon : public edm::stream::EDProducer<> {
public:
  MuonSeedsSelectorFromL1TkMuon(const edm::ParameterSet& iConfig);
  ~MuonSeedsSelectorFromL1TkMuon() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  struct SeedQuality {
    size_t seedIndex;
    float quality;
    int nHits;
    int nPixelHits;
    int nOTHits;
    int nConsistentHits;
    float dR2;
    float normalizedPtDiff;

    // Lower quality is better
    bool operator<(const SeedQuality& other) const { return quality < other.quality; }
  };

  // Evaluate seed quality based on hits
  SeedQuality evaluateSeedQuality(const TrajectorySeed& seed,
                                  const size_t seedIndex,
                                  const float deltaR2,
                                  const float l1TkMuPt,
                                  const float normPtDiff,
                                  const std::vector<DetId>& stubDetIds,
                                  const GlobalTrackingGeometry* geometry,
                                  const MagneticField* magField,
                                  const TrackerTopology* tkTopo) const;

  bool isHitConsistent(const TrackingRecHit& hit,
                       const std::vector<DetId>& stubDetIds,
                       const TrackerTopology* tkTopo) const;

  // Tokens
  const edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuToken_;
  const edm::EDGetTokenT<TrajectorySeedCollection> seedToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geometryToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tkTopoToken_;

  // Selection parameters
  const float l1TkMuMinPt_;
  // Kinematics matching parameters
  const float maxDrForPreselection_;
  const float pTCompatibility_;
  const float maxSeedPtForCompatibilityCheck_;
  // Hits parameters
  const int minNHits_;
  const int minNPixelHits_;
  const int minNConsistentHits_;

  const float maxCombinedQuality_;
  const size_t nSeedsToKeep_;
};

MuonSeedsSelectorFromL1TkMuon::MuonSeedsSelectorFromL1TkMuon(const edm::ParameterSet& iConfig)
    : l1TkMuToken_{consumes<l1t::TrackerMuonCollection>(
          iConfig.getParameter<edm::InputTag>("L1TkMuonInputCollection"))},
      seedToken_{consumes<TrajectorySeedCollection>(iConfig.getParameter<edm::InputTag>("SeedInputCollection"))},
      magFieldToken_{esConsumes<MagneticField, IdealMagneticFieldRecord>()},
      geometryToken_{esConsumes<GlobalTrackingGeometry, GlobalTrackingGeometryRecord>()},
      tkTopoToken_{esConsumes<TrackerTopology, TrackerTopologyRcd>()},
      l1TkMuMinPt_{static_cast<float>(iConfig.getParameter<double>("L1TkMuMinPt"))},
      maxDrForPreselection_{static_cast<float>(iConfig.getParameter<double>("maxDrForPreselection"))},
      pTCompatibility_{static_cast<float>(iConfig.getParameter<double>("pTCompatibility"))},
      maxSeedPtForCompatibilityCheck_{
          static_cast<float>(iConfig.getParameter<double>("maxSeedPtForCompatibilityCheck"))},
      minNHits_{iConfig.getParameter<int>("minNHits")},
      minNPixelHits_{iConfig.getParameter<int>("minNPixelHits")},
      minNConsistentHits_{iConfig.getParameter<int>("minNConsistentHits")},
      maxCombinedQuality_{static_cast<float>(iConfig.getParameter<double>("maxCombinedQuality"))},
      nSeedsToKeep_{static_cast<size_t>(iConfig.getParameter<uint>("nSeedsToKeep"))} {
  produces<TrajectorySeedCollection>();
}

bool MuonSeedsSelectorFromL1TkMuon::isHitConsistent(const TrackingRecHit& hit,
                                                    const std::vector<DetId>& stubDetIds,
                                                    const TrackerTopology* tkTopo) const {
  if (!hit.isValid())
    return false;
  DetId hitId = hit.geographicalId();
  if (hitId.subdetId() != StripSubdetector::TOB && hitId.subdetId() != StripSubdetector::TID) {
    // Only consider strip hits for matching to L1TkMu stubs
    return false;
  }
  auto stackId = tkTopo->stack(hitId);

  bool consistent = std::find(stubDetIds.begin(), stubDetIds.end(), stackId) != stubDetIds.end();
  if (consistent) {
    LogDebug(metname) << "    Hit with detId " << hitId << " is consistent with L1TkMu stub in stack " << stackId;
  }
  return consistent;
}

MuonSeedsSelectorFromL1TkMuon::SeedQuality MuonSeedsSelectorFromL1TkMuon::evaluateSeedQuality(
    const TrajectorySeed& seed,
    const size_t seedIndex,
    const float deltaR2,
    const float l1TkMuPt,
    const float normPtDiff,
    const std::vector<DetId>& stubDetIds,
    const GlobalTrackingGeometry* geometry,
    const MagneticField* magField,
    const TrackerTopology* tkTopo) const {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonSeedsSelectorFromL1TkMuon";

  SeedQuality squal;
  squal.seedIndex = seedIndex;
  squal.nHits = 0;
  squal.nPixelHits = 0;
  squal.nOTHits = 0;
  squal.nConsistentHits = 0;
  squal.quality = 0.0;
  squal.dR2 = deltaR2;
  squal.normalizedPtDiff = normPtDiff;

  // Get seed state
  const GeomDet* seedDet = geometry->idToDet(seed.startingState().detId());
  if (!seedDet) {
    // Invalid seed
    squal.quality = 1e6;
    return squal;
  }

  TrajectoryStateOnSurface seedTSOS =
      trajectoryStateTransform::transientState(seed.startingState(), &(seedDet->surface()), magField);

  if (!seedTSOS.isValid()) {
    // Invalid seed
    squal.quality = 1e6;
    return squal;
  }

  // Hit counting and consistency evaluation
  int nPixel = 0, nOT = 0;
  for (auto const& recHit : seed.recHits()) {
    DetId hitId = recHit.geographicalId();
    if (hitId.subdetId() == PixelSubdetector::PixelBarrel || hitId.subdetId() == PixelSubdetector::PixelEndcap) {
      LogDebug(metname) << " Found pixel hit " << hitId.rawId() << " in seed " << seedIndex;
      nPixel++;
    } else if (hitId.subdetId() == StripSubdetector::TOB || hitId.subdetId() == StripSubdetector::TID) {
      LogDebug(metname) << " Found OT hit " << hitId.rawId() << " in seed " << seedIndex;
      nOT++;
      if (isHitConsistent(recHit, stubDetIds, tkTopo)) {
        squal.nConsistentHits++;
      }
    } else {
      edm::LogWarning(metname) << "Seed has hit with unexpected subdetector ID " << hitId.subdetId() << " (det "
                               << hitId.det() << ")";
    }
  }

  squal.nHits = nPixel + nOT;
  squal.nPixelHits = nPixel;
  squal.nOTHits = nOT;

  // Chi-square-like quality metric (lower is better, 0 is perfect)

  // Term 1: deltaR normalized to maximum allowed
  float termDR = squal.dR2 / (maxDrForPreselection_ * maxDrForPreselection_);

  // Term 2: Hit consistency (only for OT hits)
  float termConsistency = squal.nOTHits > 0 ? static_cast<float>(squal.nConsistentHits) / squal.nOTHits : 0.0f;

  // Term 3: Total hit quality
  float termNHits = squal.nHits > 0 ? 1.0f / squal.nHits : 1e6;

  // Term 4: normalized pT difference (only for low pT L1Tk Muons)
  float termPtDiff = seedTSOS.globalMomentum().perp() < maxSeedPtForCompatibilityCheck_
                         ? squal.normalizedPtDiff / pTCompatibility_
                         : 0.0f;

  // Combined quality: sum of terms (lower is better)
  squal.quality = termDR - termConsistency + termNHits + termPtDiff;

  LogDebug(metname) << "Seed " << seedIndex << " quality evaluation:"
                    << "\n  nHits=" << squal.nHits << " nPixelHits=" << squal.nPixelHits
                    << " nConsistentHits=" << squal.nConsistentHits << " dR2=" << squal.dR2
                    << "\n  combined quality=" << squal.quality;

  return squal;
}

void MuonSeedsSelectorFromL1TkMuon::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonSeedsSelectorFromL1TkMuon";

  // Output collection
  std::unique_ptr<TrajectorySeedCollection> outputSeeds = std::make_unique<TrajectorySeedCollection>();

  // Get input collections
  auto const l1TkMuCollectionH = iEvent.getHandle(l1TkMuToken_);
  auto const seedCollectionH = iEvent.getHandle(seedToken_);

  if (!l1TkMuCollectionH.isValid() || !seedCollectionH.isValid()) {
    edm::LogWarning(metname) << "Input collection not valid! Returning empty collection.";
    iEvent.put(std::move(outputSeeds));
    return;
  }

  if (l1TkMuCollectionH->empty() || seedCollectionH->empty()) {
    LogDebug(metname) << "Input collection empty! Returning empty collection.";
    iEvent.put(std::move(outputSeeds));
    return;
  }

  // Get event setup
  const MagneticField* magField = &iSetup.getData(magFieldToken_);
  const GlobalTrackingGeometry* geometry = &iSetup.getData(geometryToken_);
  const TrackerTopology* tkTopo = &iSetup.getData(tkTopoToken_);

  // Store indices of matched seeds
  std::set<size_t> seedsToKeep;

  LogDebug(metname) << "Processing " << l1TkMuCollectionH->size() << " L1TkMuons and " << seedCollectionH->size()
                    << " seeds";

  // Loop over L1TkMuons
  for (size_t l1TkMuIndex = 0; l1TkMuIndex != l1TkMuCollectionH->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuCollectionH, l1TkMuIndex);

    // Get L1TkMu parameters
    float l1TkMuEta, l1TkMuPhi, l1TkMuPt;
    auto trkPtr = l1TkMuRef->trkPtr();
    std::vector<DetId> stubDetIds;
    stubDetIds.reserve(trkPtr->getStubRefs().size());

    if (!trkPtr.isNull()) {
      // Prefer tracker-based coordinates
      l1TkMuEta = trkPtr->momentum().eta();
      l1TkMuPhi = trkPtr->momentum().phi();
      l1TkMuPt = trkPtr->momentum().perp();
      // Skip L1TkMu below pT threshold
      if (l1TkMuPt < l1TkMuMinPt_) {
        LogDebug(metname) << "Skipping L1TkMuon " << l1TkMuIndex << " with pT=" << l1TkMuPt << " below threshold of "
                          << l1TkMuMinPt_;
        continue;
      }

      LogDebug(metname) << "Processing L1TkMuon " << l1TkMuIndex << " with eta=" << l1TkMuEta << ", phi=" << l1TkMuPhi
                        << ", pT=" << l1TkMuPt << ". Found TTTrack with " << trkPtr->getStubRefs().size() << " stubs:";

      for (const auto& stubRef : trkPtr->getStubRefs()) {
        auto stubDetId = stubRef->getDetId();
        stubDetIds.push_back(stubDetId);
        LogDebug(metname) << "  Stub detId: " << stubDetId.rawId() << " (subdet " << stubDetId.subdetId() << ")\n";
      }

    } else {
      // Fallback to standalone muon coordinates
      // TODO figure out what to do here (no L1Tk match)
      edm::LogWarning(metname) << "L1TkMuon with index " << l1TkMuIndex
                               << " has no tracker pointer, using muon system coordinates";
      l1TkMuEta = l1TkMuRef->phEta();
      l1TkMuPhi = l1TkMuRef->phPhi();
      l1TkMuPt = l1TkMuRef->phPt();
    }

    // Vector to hold seed quality information for this L1TkMu
    std::vector<SeedQuality> candidateSeedsWithPriority;
    std::vector<SeedQuality> candidateSeeds;

    for (size_t seedIndex = 0; seedIndex != seedCollectionH->size(); ++seedIndex) {
      const TrajectorySeed& seed = (*seedCollectionH)[seedIndex];

      // Skip seeds with too few hits
      if (static_cast<int>(seed.nHits()) < minNHits_) {
        continue;
      }

      // Convert seed state to TSOS
      const GeomDet* seedDet = geometry->idToDet(seed.startingState().detId());
      if (!seedDet) {
        continue;
      }

      TrajectoryStateOnSurface seedTSOS =
          trajectoryStateTransform::transientState(seed.startingState(), &(seedDet->surface()), magField);

      if (!seedTSOS.isValid()) {
        continue;
      }

      // Get seed kinematics
      float seedEta = seedTSOS.globalMomentum().eta();
      float seedPhi = seedTSOS.globalMomentum().phi();
      float seedPt = seedTSOS.globalMomentum().perp();

      // Loose dR cone for pre-selection
      float dR2 = deltaR2(l1TkMuEta, l1TkMuPhi, seedEta, seedPhi);

      if (dR2 > maxDrForPreselection_ * maxDrForPreselection_) {
        continue;
      }

      float normPtDiff = std::abs(seedPt - l1TkMuPt) / l1TkMuPt;
      LogDebug(metname) << "Processing seed " << seedIndex << ": eta=" << seedEta << " phi=" << seedPhi
                        << " pT=" << seedPt << " dR2 to L1TkMu " << l1TkMuIndex << " is " << dR2
                        << ", normalised pT difference=" << normPtDiff;

      if (seedPt < maxSeedPtForCompatibilityCheck_ && normPtDiff > pTCompatibility_) {
        LogDebug(metname) << "Seed " << seedIndex
                          << " rejected: pT difference with L1TkMu is too large (|seedPt - l1TkMuPt|/l1TkMuPt = "
                          << normPtDiff << " > " << pTCompatibility_;
        continue;
      }

      // Evaluate seed quality
      SeedQuality squal =
          evaluateSeedQuality(seed, seedIndex, dR2, l1TkMuPt, normPtDiff, stubDetIds, geometry, magField, tkTopo);

      if (squal.quality > maxCombinedQuality_) {
        LogDebug(metname) << "Seed " << seedIndex << " rejected: poor combined quality (" << squal.quality << " > "
                          << maxCombinedQuality_ << ")\n";
        continue;
      }

      // Selection passed, add to candidate list for arbitration
      // Prefer seeds with pixel hits and compatibility with L1TkMu stubs
      if (squal.nPixelHits > minNPixelHits_ && squal.nConsistentHits >= minNConsistentHits_) {
        candidateSeedsWithPriority.push_back(squal);
      } else {
        candidateSeeds.push_back(squal);
      }
      LogDebug(metname) << "Seed " << seedIndex << " passed pre-selection for L1TkMu " << l1TkMuIndex
                        << " with quality " << squal.quality;
    }  // End loop on seeds

    // Arbitration: first take best seeds with pixels and OT hits consistent with L1 stubs
    std::vector<SeedQuality> finalSeeds;
    finalSeeds.reserve(nSeedsToKeep_);

    if (!candidateSeedsWithPriority.empty()) {
      // Partial sort to get only the best nSeedsToKeep_
      size_t nToSort = std::min(nSeedsToKeep_, candidateSeedsWithPriority.size());
      std::partial_sort(candidateSeedsWithPriority.begin(),
                        candidateSeedsWithPriority.begin() + nToSort,
                        candidateSeedsWithPriority.end());

      LogDebug(metname) << "L1TkMu " << l1TkMuIndex << ": found " << candidateSeedsWithPriority.size()
                        << " candidate seeds with priority, keeping best " << nToSort;

      for (size_t i = 0; i < nToSort; ++i) {
        finalSeeds.push_back(candidateSeedsWithPriority[i]);
        auto candidate = candidateSeedsWithPriority[i];

        LogDebug(metname) << "  Keeping priority seed " << candidate.seedIndex << " with quality=" << candidate.quality
                          << " (nHits=" << candidate.nHits << ", nPixel=" << candidate.nPixelHits
                          << ", nConsistent=" << candidate.nConsistentHits << ", dR2=" << candidate.dR2
                          << ", normPtDiff=" << candidate.normalizedPtDiff << ")";
      }
    }
    if (finalSeeds.size() < nSeedsToKeep_ and not candidateSeeds.empty()) {
      size_t nRemaining = nSeedsToKeep_ - finalSeeds.size();
      size_t nToSort = std::min(nRemaining, candidateSeeds.size());
      std::partial_sort(candidateSeeds.begin(), candidateSeeds.begin() + nToSort, candidateSeeds.end());

      LogDebug(metname) << "L1TkMu " << l1TkMuIndex << ": found " << candidateSeeds.size()
                        << " additional candidate seeds, keeping best " << nToSort;

      for (size_t i = 0; i < nToSort; ++i) {
        finalSeeds.push_back(candidateSeeds[i]);
        auto candidate = candidateSeeds[i];

        LogDebug(metname) << "  Keeping additional seed " << candidate.seedIndex
                          << " with quality=" << candidate.quality << " (nHits=" << candidate.nHits
                          << ", nPixel=" << candidate.nPixelHits << ", nConsistent=" << candidate.nConsistentHits
                          << ", dR2=" << candidate.dR2 << ", normPtDiff=" << candidate.normalizedPtDiff << ")";
      }
    }
    if (finalSeeds.empty()) {
      std::cout << "WARNING - L1TkMu " << l1TkMuIndex << ": no candidate seeds passed selection\n";
    }
    // Add selected seeds to output
    for (const auto& seed : finalSeeds) {
      seedsToKeep.insert(seed.seedIndex);

      LogDebug(metname) << "L1TkMu " << l1TkMuIndex << " keeping seed " << seed.seedIndex
                        << " with quality=" << seed.quality << " (nHits=" << seed.nHits
                        << ", nPixel=" << seed.nPixelHits << ", nConsistent=" << seed.nConsistentHits
                        << ", dR2=" << seed.dR2 << ", normPtDiff=" << seed.normalizedPtDiff << ")";
    }

  }  // End L1TkMuon loop

  LogDebug(metname) << "Seed selection from " << l1TkMuCollectionH->size() << " L1TkMuons and "
                    << seedCollectionH->size() << " input seeds: selected " << seedsToKeep.size() << " seeds";

  // Fill output collection with selected seeds
  outputSeeds->reserve(seedsToKeep.size());
  for (size_t idx : seedsToKeep) {
    outputSeeds->push_back((*seedCollectionH)[idx]);
  }
  iEvent.put(std::move(outputSeeds));
}

void MuonSeedsSelectorFromL1TkMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  // Input collections
  desc.add<edm::InputTag>("L1TkMuonInputCollection", edm::InputTag("l1tTkMuonsGmt"))
      ->setComment("Input collection of L1TkMuons to match to");
  desc.add<edm::InputTag>("SeedInputCollection", edm::InputTag("hltInitialStepTrajectorySeedsLST"))
      ->setComment("Input seed collection to be filtered based on L1TkMuon matching");

  // L1TkMuon selection
  desc.add<double>("L1TkMuMinPt", 0.0)->setComment("Minimum pT for L1TkMuons to be considered for matching");

  // Geometric matching parameters
  desc.add<double>("maxDrForPreselection", 0.4)
      ->setComment("Maximum deltaR for loose pre-selection of seeds (cone around L1TkMuon)");
  desc.add<double>("pTCompatibility", 0.5)
      ->setComment("Maximum allowed relative pT difference between seed and L1TkMuon (|seedPt - l1TkMuPt|/l1TkMuPt)");
  desc.add<double>("maxSeedPtForCompatibilityCheck", 25.0)
      ->setComment("Maximum L1Tk Muon pT for which to apply the pT compatibility check with seeds");
  // Seed quality requirements
  desc.add<int>("minNHits", 4)->setComment("Minimum number of hits in seed");
  desc.add<int>("minNPixelHits", 2)->setComment("Minimum number of pixel hits in seed");
  desc.add<int>("minNConsistentHits", 4)
      ->setComment("Minimum number of OT hits consistent with L1TkMuon TTTrack stubs");

  // Quality metric parameters
  desc.add<double>("maxCombinedQuality", 2.5)
      ->setComment("Maximum combined quality metric for seed to be considered (lower is better)");

  // Arbitration
  desc.add<uint32_t>("nSeedsToKeep", 2)
      ->setComment("Maximum number of seeds to keep per L1TkMuon based on quality ranking");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonSeedsSelectorFromL1TkMuon);
