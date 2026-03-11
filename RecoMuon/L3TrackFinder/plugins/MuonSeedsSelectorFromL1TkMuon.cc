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
 *  3. Hit-by-hit consistency check with L1TkMuon direction
 *  4. Arbitration: keep only the N best seeds per L1TkMuon
 *
 *  Adapted from MuonTracksSelectorFromL1TkMuon
 *  Inspired by MuonIdProducer::fillArbitrationInfo
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
    int nConsistentHits;
    float dR2;

    // Lower quality is better
    bool operator<(const SeedQuality& other) const { return quality < other.quality; }
  };

  // Evaluate seed quality based on hits
  SeedQuality evaluateSeedQuality(const TrajectorySeed& seed,
                                  size_t seedIndex,
                                  float l1TkMuEta,
                                  float l1TkMuPhi,
                                  const GlobalTrackingGeometry* geometry,
                                  const MagneticField* magField) const;

  // Check if a hit is consistent with L1TkMu direction
  bool isHitConsistent(const TrackingRecHit& hit,
                       float l1TkMuEta,
                       float l1TkMuPhi,
                       const GlobalTrackingGeometry* geometry) const;

  // Count hits in different detector categories
  void categorizeHits(const TrajectorySeed& seed, int& nPixel, int& nStrip) const;

  // Tokens
  const edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuToken_;
  const edm::EDGetTokenT<TrajectorySeedCollection> seedToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geometryToken_;

  // Selection parameters
  const float l1TkMuMinPt_;
  // Wide cone for initial pre-selection
  const float maxDrForPreselection_;
  // Tigheter requirement for hit consistency with L1TkMu direction
  const float maxDrForHitConsistency_;
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
      l1TkMuMinPt_{static_cast<float>(iConfig.getParameter<double>("L1TkMuMinPt"))},
      maxDrForPreselection_{static_cast<float>(iConfig.getParameter<double>("maxDrForPreselection"))},
      maxDrForHitConsistency_{static_cast<float>(iConfig.getParameter<double>("maxDrForHitConsistency"))},
      minNHits_{iConfig.getParameter<int>("minNHits")},
      minNPixelHits_{iConfig.getParameter<int>("minNPixelHits")},
      minNConsistentHits_{iConfig.getParameter<int>("minNConsistentHits")},
      maxCombinedQuality_{static_cast<float>(iConfig.getParameter<double>("maxCombinedQuality"))},
      nSeedsToKeep_{static_cast<size_t>(iConfig.getParameter<uint>("nSeedsToKeep"))} {
  produces<TrajectorySeedCollection>();
}

void MuonSeedsSelectorFromL1TkMuon::categorizeHits(const TrajectorySeed& seed, int& nPixel, int& nStrip) const {
  nPixel = 0;
  nStrip = 0;

  for (auto const& recHit : seed.recHits()) {
    if (!recHit.isValid())
      continue;

    DetId hitId = recHit.geographicalId();
    if (hitId.subdetId() == PixelSubdetector::PixelBarrel || hitId.subdetId() == PixelSubdetector::PixelEndcap) {
      nPixel++;
    } else if (hitId.subdetId() == StripSubdetector::TOB ||
               hitId.subdetId() == StripSubdetector::TID) {  //(hitId.det() == DetId::Tracker) {
      nStrip++;
    } else {
      std::cout << "WARNING: Seed has hit with unexpected subdetector ID " << hitId.subdetId() << " (det "
                << hitId.det() << ")\n";
    }
  }
}

bool MuonSeedsSelectorFromL1TkMuon::isHitConsistent(const TrackingRecHit& hit,
                                                    float l1TkMuEta,
                                                    float l1TkMuPhi,
                                                    const GlobalTrackingGeometry* geometry) const {
  if (!hit.isValid())
    return false;

  const GeomDet* hitDet = geometry->idToDet(hit.geographicalId());
  if (!hitDet)
    return false;

  GlobalPoint hitPos = hitDet->toGlobal(hit.localPosition());
  float hitDr2 = reco::deltaR2(l1TkMuEta, l1TkMuPhi, hitPos.eta(), hitPos.phi());

  return (hitDr2 < maxDrForHitConsistency_ * maxDrForHitConsistency_);
}

MuonSeedsSelectorFromL1TkMuon::SeedQuality MuonSeedsSelectorFromL1TkMuon::evaluateSeedQuality(
    const TrajectorySeed& seed,
    size_t seedIndex,
    float l1TkMuEta,
    float l1TkMuPhi,
    const GlobalTrackingGeometry* geometry,
    const MagneticField* magField) const {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonSeedsSelectorFromL1TkMuon";

  SeedQuality squal;
  squal.seedIndex = seedIndex;
  squal.nHits = 0;
  squal.nPixelHits = 0;
  squal.nConsistentHits = 0;
  squal.quality = 0.0;

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

  float seedEta = seedTSOS.globalMomentum().eta();
  float seedPhi = seedTSOS.globalMomentum().phi();
  squal.dR2 = reco::deltaR2(l1TkMuEta, l1TkMuPhi, seedEta, seedPhi);

  // Categorize hits
  int nPixel = 0, nStrip = 0;
  categorizeHits(seed, nPixel, nStrip);
  squal.nHits = nPixel + nStrip;
  squal.nPixelHits = nPixel;

  // Count hits consistent with L1TkMu direction
  for (auto const& recHit : seed.recHits()) {
    if (isHitConsistent(recHit, l1TkMuEta, l1TkMuPhi, geometry)) {
      squal.nConsistentHits++;
    }
  }
  // Chi-square-like quality metric (lower is better, 0 is perfect)
  // Note: pixel hits are NOT included here - they're used for categorization instead

  // Term 1: deltaR normalized to maximum allowed
  float termDR = squal.dR2 / (maxDrForPreselection_ * maxDrForPreselection_);

  // Term 2: Hit consistency (expect 100% consistency)
  float fractionConsistent = squal.nHits > 0 ? static_cast<float>(squal.nConsistentHits) / squal.nHits : 0.0f;
  float termConsistency = (1.0f - fractionConsistent) * (1.0f - fractionConsistent);

  // Term 3: Total hit quality (more hits -> better quality)
  float termNHits = squal.nHits > 0 ? 1.0f / squal.nHits : 1e6;

  // Combined quality: sum of terms (lower is better)
  squal.quality = termDR + termConsistency + termNHits;

  /*
  // Attempt at chi-2 like quality metric
  // Term 1: deltaR normalized to maximum allowed
  float termDR = squal.dR2 / (maxDrForPreselection_ * maxDrForPreselection_);

  // Term 2: Hit consistency as a "chi-square"
  float fractionConsistent = squal.nHits > 0 ? static_cast<float>(squal.nConsistentHits) / squal.nHits : 0.0f;
  float termConsistency = (1.0f - fractionConsistent) * (1.0f - fractionConsistent);

  // Term 3: Number of hits (more hits -> better quality)
  float termNHits = 1.0f / squal.nHits;

  // Term 4: Number of pixel hits (more pixel hits -> better quality)
  float termNPixelHits = 1.0f / (1.0f + squal.nPixelHits);

  // Combined quality: sum of chi-square-like terms
  // Good seeds have quality close to 0, bad seeds have large quality
  squal.quality = termDR + termConsistency + termNHits + termNPixelHits;
*/
  // Original combined quality metric (lower is better)
  // More hits -> better
  // More pixel hits -> better
  // More consistent hits -> better
  // Smaller dR2 -> better
  //squal.quality = -weightNHits_ * squal.nHits - weightNPixelHits_ * squal.nPixelHits -
  //                weightNConsistentHits_ * squal.nConsistentHits + weightDR_ * squal.dR2;

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

  // Store indices of matched seeds
  std::set<size_t> seedsToKeep;

  std::cout << "Processing " << l1TkMuCollectionH->size() << " L1TkMuons and " << seedCollectionH->size() << " seeds\n";

  // Loop over L1TkMuons
  for (size_t l1TkMuIndex = 0; l1TkMuIndex != l1TkMuCollectionH->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuCollectionH, l1TkMuIndex);

    // Get L1TkMu parameters
    float l1TkMuEta, l1TkMuPhi, l1TkMuPt;
    auto trkPtr = l1TkMuRef->trkPtr();

    if (!trkPtr.isNull()) {
      // Prefer tracker-based coordinates
      l1TkMuEta = trkPtr->momentum().eta();
      l1TkMuPhi = trkPtr->momentum().phi();
      l1TkMuPt = trkPtr->momentum().perp();
    } else {
      // Fallback to standalone muon coordinates
      edm::LogWarning(metname) << "L1TkMuon with index " << l1TkMuIndex
                               << " has no tracker pointer, using muon system coordinates";
      l1TkMuEta = l1TkMuRef->phEta();
      l1TkMuPhi = l1TkMuRef->phPhi();
      l1TkMuPt = l1TkMuRef->phPt();
    }

    // Basic kinematic selection on L1TkMu
    if (l1TkMuPt < l1TkMuMinPt_) {
      continue;
    }

    // Vector to hold seed quality information for this L1TkMu
    std::vector<SeedQuality> candidateSeedsWithPixels;
    std::vector<SeedQuality> candidateSeedsNoPixels;

    // Loop over seeds
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

      // Loose dR cone for pre-selection
      float dR2 = deltaR2(l1TkMuEta, l1TkMuPhi, seedEta, seedPhi);
      if (dR2 > maxDrForPreselection_ * maxDrForPreselection_) {
        continue;
      }

      // Evaluate seed quality
      SeedQuality squal = evaluateSeedQuality(seed, seedIndex, l1TkMuEta, l1TkMuPhi, geometry, magField);

      if (squal.nConsistentHits < minNConsistentHits_) {
        LogDebug(metname) << "Seed " << seedIndex << " rejected: insufficient consistent hits ("
                          << squal.nConsistentHits << " < " << minNConsistentHits_ << ")";
        continue;
      }

      if (squal.quality > maxCombinedQuality_) {
        LogDebug(metname) << "Seed " << seedIndex << " rejected: poor combined quality (" << squal.quality << " > "
                          << maxCombinedQuality_ << ")";
        continue;
      }

      // Selection passed, add to candidate list for arbitration
      if (squal.nPixelHits > minNPixelHits_) {
        candidateSeedsWithPixels.push_back(squal);
      } else {
        candidateSeedsNoPixels.push_back(squal);
      }
      LogDebug(metname) << "Seed " << seedIndex << " passed pre-selection for L1TkMu " << l1TkMuIndex
                        << " with quality " << squal.quality;
    }  // End seed loop

    // Arbitration: first take best seeds with pixels, then fill remaining slots with no-pixel seeds
    std::vector<SeedQuality> finalSeeds;
    finalSeeds.reserve(nSeedsToKeep_);

    if (!candidateSeedsWithPixels.empty()) {
      // Partial sort to get only the best nSeedsToKeep_
      size_t nToSort = std::min(nSeedsToKeep_, candidateSeedsWithPixels.size());
      std::partial_sort(
          candidateSeedsWithPixels.begin(), candidateSeedsWithPixels.begin() + nToSort, candidateSeedsWithPixels.end());

      LogDebug(metname) << "L1TkMu " << l1TkMuIndex << ": found " << candidateSeedsWithPixels.size()
                        << " candidate seeds with pixel hits, keeping best " << nToSort;

      for (size_t i = 0; i < nToSort; ++i) {
        finalSeeds.push_back(candidateSeedsWithPixels[i]);
        auto candidate = candidateSeedsWithPixels[i];

        LogDebug(metname) << "  Keeping pixel seed " << candidate.seedIndex << " with quality=" << candidate.quality
                          << " (nHits=" << candidate.nHits << ", nPixel=" << candidate.nPixelHits
                          << ", nConsistent=" << candidate.nConsistentHits << ", dR2=" << candidate.dR2 << ")";
      }
    }
    if (finalSeeds.size() < nSeedsToKeep_ and not candidateSeedsNoPixels.empty()) {
      size_t nRemaining = nSeedsToKeep_ - finalSeeds.size();
      size_t nToSort = std::min(nRemaining, candidateSeedsNoPixels.size());
      std::partial_sort(
          candidateSeedsNoPixels.begin(), candidateSeedsNoPixels.begin() + nToSort, candidateSeedsNoPixels.end());

      LogDebug(metname) << "L1TkMu " << l1TkMuIndex << ": found " << candidateSeedsNoPixels.size()
                        << " candidate seeds without pixel hits, keeping best " << nToSort;

      for (size_t i = 0; i < nToSort; ++i) {
        finalSeeds.push_back(candidateSeedsNoPixels[i]);
        auto candidate = candidateSeedsNoPixels[i];

        LogDebug(metname) << "  Keeping no-pixel seed " << candidate.seedIndex << " with quality=" << candidate.quality
                          << " (nHits=" << candidate.nHits << ", nPixel=" << candidate.nPixelHits
                          << ", nConsistent=" << candidate.nConsistentHits << ", dR2=" << candidate.dR2 << ")";
      }
    }
    // Add selected seeds to output
    for (const auto& squal : finalSeeds) {
      seedsToKeep.insert(squal.seedIndex);

      std::cout << "L1TkMu " << l1TkMuIndex << " keeping seed " << squal.seedIndex << " with quality=" << squal.quality
                << " (nHits=" << squal.nHits << ", nPixel=" << squal.nPixelHits
                << ", nConsistent=" << squal.nConsistentHits << ", dR2=" << squal.dR2 << ")\n";
    }
  }  // End L1TkMuon loop

  std::cout << "Seed selection from " << l1TkMuCollectionH->size() << " L1TkMuons and " << seedCollectionH->size()
            << " input seeds: selected " << seedsToKeep.size() << " seeds\n";

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
  desc.add<double>("maxDrForHitConsistency", 0.1)
      ->setComment("Maximum deltaR for individual hit consistency with L1TkMuon direction");

  // Seed quality requirements
  desc.add<int>("minNHits", 4)->setComment("Minimum number of hits in seed");
  desc.add<int>("minNPixelHits", 2)->setComment("Minimum number of pixel hits in seed");
  desc.add<int>("minNConsistentHits", 1)
      ->setComment("Minimum number of hits consistent with L1TkMuon direction (within maxDrForHitConsistency)");

  // Quality metric parameters
  desc.add<double>("maxCombinedQuality", 5.0)
      ->setComment("Maximum combined quality metric for seed to be considered (lower is better)");

  // Arbitration
  desc.add<uint32_t>("nSeedsToKeep", 5)
      ->setComment("Maximum number of seeds to keep per L1TkMuon based on quality ranking");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonSeedsSelectorFromL1TkMuon);
