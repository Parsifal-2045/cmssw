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
 *  For high-pT L1Tk Muons and seeds only dR check is performed and all seeds with pT above
 *  threshold found in the cone are kept (no arbitration or other limits)
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
#include <unordered_set>
#include <vector>
#include <algorithm>

class MuonSeedsSelectorFromL1TkMuon : public edm::stream::EDProducer<> {
public:
  MuonSeedsSelectorFromL1TkMuon(const edm::ParameterSet& iConfig);
  ~MuonSeedsSelectorFromL1TkMuon() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  struct SeedInfo {
    // DetIds of the OT stacks crossed by the seed hits, used for consistency check with L1TkMu stubs
    std::vector<DetId> otStackIds;

    // Index of the seed in the input collection
    size_t seedIndex;

    // Kinematics
    float eta;
    float phi;
    float pt;

    // Hits
    int nHits;
    int nPixelHits;
    int nOTHits;
    int nConsistentHits;

    // Flag seeds to be skipped based on basic criteria (e.g. invalid TSOS)
    bool skip = false;
  };

  struct SeedQuality {
    // Index of the seed in the input collection
    size_t seedIndex;

    // Kinematics matching
    float dR2;
    float normalizedPtDiff;

    // Hits
    int nHits;
    int nPixelHits;
    int nOTHits;
    int nConsistentHits;

    // Quality terms
    float termDr;
    float termConsistency;
    float termNHits;
    float termPtDiff;
    float quality;

    bool operator<(const SeedQuality& other) const { return quality < other.quality; }
  };

  // Cache seed information to avoid redundant calculations during quality evaluation
  SeedInfo fillSeedInfo(const TrajectorySeed& seed,
                        const size_t seedIndex,
                        const GlobalTrackingGeometry* geometry,
                        const MagneticField* magField,
                        const TrackerTopology* tkTopo) const;

  // Check how many OT hits in a seed are consistent with L1TkMu stubs
  int countConsistentHits(const std::vector<DetId>& hitIds, const std::unordered_set<DetId>& stubDetIds) const;

  // Compute seed quality based on kinematic matching and hit consistency with L1TkMu stubs
  SeedQuality evaluateSeedQuality(const SeedInfo& sinfo,
                                  const float deltaR2,
                                  const float normPtDiff,
                                  const std::unordered_set<DetId>& stubDetIds) const;

  // Tokens
  const edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuToken_;
  const edm::EDGetTokenT<TrajectorySeedCollection> seedToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geometryToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tkTopoToken_;

  // Selection parameters
  const float l1TkMuMinPt_;
  const float maxDrForPreselection_;
  const float pTCompatibility_;
  const float maxPtForCompatibilityCheck_;

  // Arbitration parameters
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
      maxPtForCompatibilityCheck_{static_cast<float>(iConfig.getParameter<double>("maxPtForCompatibilityCheck"))},
      minNPixelHits_{iConfig.getParameter<int>("minNPixelHits")},
      minNConsistentHits_{iConfig.getParameter<int>("minNConsistentHits")},
      maxCombinedQuality_{static_cast<float>(iConfig.getParameter<double>("maxCombinedQuality"))},
      nSeedsToKeep_{static_cast<size_t>(iConfig.getParameter<uint>("nSeedsToKeep"))} {
  produces<TrajectorySeedCollection>();
}

MuonSeedsSelectorFromL1TkMuon::SeedInfo MuonSeedsSelectorFromL1TkMuon::fillSeedInfo(
    const TrajectorySeed& seed,
    const size_t seedIndex,
    const GlobalTrackingGeometry* geometry,
    const MagneticField* magField,
    const TrackerTopology* tkTopo) const {
  SeedInfo info;
  info.seedIndex = seedIndex;
  info.nHits = seed.nHits();
  info.otStackIds.reserve(seed.nHits());

  // Check if seed has valid starting state and can be converted to global kinematics
  const GeomDet* seedDet = geometry->idToDet(seed.startingState().detId());
  if (!seedDet) {
    info.skip = true;
    return info;
  }

  TrajectoryStateOnSurface seedTSOS =
      trajectoryStateTransform::transientState(seed.startingState(), &(seedDet->surface()), magField);
  if (!seedTSOS.isValid()) {
    info.skip = true;
    return info;
  }

  info.eta = seedTSOS.globalMomentum().eta();
  info.phi = seedTSOS.globalMomentum().phi();
  info.pt = seedTSOS.globalMomentum().perp();

  info.nPixelHits = 0;
  info.nOTHits = 0;
  info.nConsistentHits = -1;

  for (auto const& recHit : seed.recHits()) {
    if (!recHit.isValid())
      continue;

    DetId hitId = recHit.geographicalId();

    if (hitId.subdetId() == PixelSubdetector::PixelBarrel || hitId.subdetId() == PixelSubdetector::PixelEndcap) {
      info.nPixelHits++;
    } else if (hitId.subdetId() == StripSubdetector::TOB || hitId.subdetId() == StripSubdetector::TID) {
      info.nOTHits++;
      // Compute and insert stack ID for consistency check with L1TkMu stubs
      // TTStub detId corresponds to the stack (PS or SS) where the stub was created
      // The recHit detId corresponds to the specific sensor where the hit was recorded (P or S)
      // We need to convert the recHit detId to the corresponding stack detId to compare with L1TkMu stubs
      info.otStackIds.push_back(tkTopo->stack(hitId));
    } else {
      edm::LogWarning("MuonSeedsSelectorFromL1TkMuon")
          << "Seed has hit with unexpected subdetector ID " << hitId.subdetId() << " (det " << hitId.det() << ")";
    }
  }
  if (info.nHits != info.nPixelHits + info.nOTHits) {
    edm::LogWarning("MuonSeedsSelectorFromL1TkMuon")
        << "Seed hit counting inconsistency: total nHits=" << info.nHits << " but pixelHits=" << info.nPixelHits
        << " and otHits=" << info.nOTHits;
  }

  return info;
}

int MuonSeedsSelectorFromL1TkMuon::countConsistentHits(const std::vector<DetId>& hitIds,
                                                       const std::unordered_set<DetId>& stubDetIds) const {
  int nConsistent = 0;
  for (const auto& hitId : hitIds) {
    if (stubDetIds.count(hitId) > 0) {
      nConsistent++;
    }
  }

  return nConsistent;
}

MuonSeedsSelectorFromL1TkMuon::SeedQuality MuonSeedsSelectorFromL1TkMuon::evaluateSeedQuality(
    const SeedInfo& sinfo,
    const float deltaR2,
    const float normPtDiff,
    const std::unordered_set<DetId>& stubDetIds) const {
  const std::string metname = "RecoMuon|L3TrackFinder|MuonSeedsSelectorFromL1TkMuon";

  SeedQuality squal;
  squal.seedIndex = sinfo.seedIndex;

  squal.nHits = sinfo.nHits;
  squal.nPixelHits = sinfo.nPixelHits;
  squal.nOTHits = sinfo.nOTHits;
  squal.nConsistentHits = countConsistentHits(sinfo.otStackIds, stubDetIds);

  squal.dR2 = deltaR2;
  squal.normalizedPtDiff = normPtDiff;

  // Quality metric (lower is better, 0 is perfect)
  // Term 1: deltaR normalized to maximum allowed
  float termDr = squal.dR2 / (maxDrForPreselection_ * maxDrForPreselection_);
  squal.termDr = termDr;

  // Term 2: Hit consistency (only for OT hits)
  float termConsistency = squal.nOTHits > 0 ? static_cast<float>(squal.nConsistentHits) / squal.nOTHits : 0.0f;
  squal.termConsistency = termConsistency;

  // Term 3: Total hit quality
  float termNHits = squal.nHits > 0 ? 1.0f / squal.nHits : 1e6;
  squal.termNHits = termNHits;

  // Term 4: normalized pT difference (only for low pT L1Tk Muons)
  float termPtDiff = squal.normalizedPtDiff / pTCompatibility_;
  squal.termPtDiff = termPtDiff;

  // Combined quality: sum of terms (lower is better)
  squal.quality = sinfo.pt <= maxPtForCompatibilityCheck_ ? termDr - termConsistency + termNHits + termPtDiff : termDr;

  LogDebug(metname) << "Seed " << squal.seedIndex << " quality evaluation:"
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

  // Cache seeds info
  std::vector<SeedInfo> seedsInfo;
  seedsInfo.reserve(seedCollectionH->size());
  for (size_t seedIndex = 0; seedIndex != seedCollectionH->size(); ++seedIndex) {
    SeedInfo info = fillSeedInfo((*seedCollectionH)[seedIndex], seedIndex, geometry, magField, tkTopo);
    if (!info.skip) {
      seedsInfo.push_back(info);
    } else {
      LogDebug(metname) << "Skipping seed " << seedIndex << " due to invalid TSOS or geometry";
    }
  }

  // Loop over L1TkMuons
  for (size_t l1TkMuIndex = 0; l1TkMuIndex != l1TkMuCollectionH->size(); ++l1TkMuIndex) {
    l1t::TrackerMuonRef l1TkMuRef(l1TkMuCollectionH, l1TkMuIndex);

    // Get L1TkMu parameters
    float l1TkMuEta, l1TkMuPhi, l1TkMuPt;
    auto trkPtr = l1TkMuRef->trkPtr();
    std::unordered_set<DetId> stubDetIds;

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
        stubDetIds.insert(stubDetId);
        LogDebug(metname) << "  Stub detId: " << stubDetId.rawId() << " (subdet " << stubDetId.subdetId() << ")";
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

    //Loop over cached seed information
    for (const auto& sinfo : seedsInfo) {
      // Loose dR cone for pre-selection
      float dR2 = deltaR2(l1TkMuEta, l1TkMuPhi, sinfo.eta, sinfo.phi);

      if (dR2 > maxDrForPreselection_ * maxDrForPreselection_) {
        continue;
      }
      if (sinfo.pt > maxPtForCompatibilityCheck_) {
        LogDebug(metname)
            << "Seed " << sinfo.seedIndex << " geometrically compatible with L1TkMu " << l1TkMuIndex
            << " with pT above threshold for compatibility check, skipping all other compatibility requirements";
        seedsToKeep.insert(sinfo.seedIndex);
        continue;
      }

      float normPtDiff = std::abs(sinfo.pt - l1TkMuPt) / l1TkMuPt;

      LogDebug(metname) << "Processing seed " << sinfo.seedIndex << ": eta=" << sinfo.eta << " phi=" << sinfo.phi
                        << " pT=" << sinfo.pt << " dR2 to L1TkMu " << l1TkMuIndex << " is " << dR2
                        << ", normalised pT difference=" << normPtDiff;

      if (normPtDiff > pTCompatibility_) {
        LogDebug(metname) << "Seed " << sinfo.seedIndex
                          << " rejected: pT difference with L1TkMu is too large (|seedPt - l1TkMuPt|/l1TkMuPt = "
                          << normPtDiff << " > " << pTCompatibility_;
        continue;
      }

      LogDebug(metname) << "Seed " << sinfo.seedIndex << " passed pre-selection for L1TkMu " << l1TkMuIndex
                        << " with dR2=" << dR2 << " and normalized pT difference=" << normPtDiff;

      // Evaluate seed quality
      SeedQuality squal = evaluateSeedQuality(sinfo, dR2, normPtDiff, stubDetIds);

      if (squal.quality > maxCombinedQuality_) {
        LogDebug(metname) << "Seed " << sinfo.seedIndex << " rejected: poor combined quality (" << squal.quality
                          << " > " << maxCombinedQuality_ << ")";
        continue;
      }

      // Selection passed, add to candidate list for arbitration
      // Prefer seeds with pixel hits and compatibility with L1TkMu stubs
      if (squal.nPixelHits >= minNPixelHits_ && squal.nConsistentHits >= minNConsistentHits_) {
        candidateSeedsWithPriority.push_back(squal);
      } else {
        candidateSeeds.push_back(squal);
      }
      LogDebug(metname) << "Seed " << sinfo.seedIndex << " passed selection for L1TkMu " << l1TkMuIndex
                        << " with quality " << squal.quality << " (nHits=" << squal.nHits
                        << ", nPixelHits=" << squal.nPixelHits << ", nConsistentHits=" << squal.nConsistentHits
                        << ", dR2=" << squal.dR2 << ", normPtDiff=" << squal.normalizedPtDiff << ")";
    }  // End loop on seeds

    // Only arbitrate low pT seeds, high pT seeds already added to output collection
    if (l1TkMuPt <= maxPtForCompatibilityCheck_) {
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

          LogDebug(metname) << "  Keeping priority seed " << candidate.seedIndex
                            << " with quality=" << candidate.quality << " (nHits=" << candidate.nHits
                            << ", nPixel=" << candidate.nPixelHits << ", nConsistent=" << candidate.nConsistentHits
                            << ", dR2=" << candidate.dR2 << ", normPtDiff=" << candidate.normalizedPtDiff << ")";
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
        LogDebug(metname) << "L1TkMu " << l1TkMuIndex << ": no candidate seed passed selection";
      }
      // Add selected seeds to output
      for (const auto& seed : finalSeeds) {
        seedsToKeep.insert(seed.seedIndex);

        LogDebug(metname) << "L1TkMu " << l1TkMuIndex << " keeping seed " << seed.seedIndex
                          << " with quality=" << seed.quality << " (nHits=" << seed.nHits
                          << ", nPixel=" << seed.nPixelHits << ", nConsistent=" << seed.nConsistentHits
                          << ", dR2=" << seed.dR2 << ", normPtDiff=" << seed.normalizedPtDiff << ")";
      }
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

  // Matching parameters
  desc.add<double>("maxDrForPreselection", 0.4)
      ->setComment("Maximum deltaR for loose pre-selection of seeds (cone around L1TkMuon)");
  desc.add<double>("pTCompatibility", 0.5)
      ->setComment("Maximum allowed relative pT difference between seed and L1TkMuon (|seedPt - l1TkMuPt|/l1TkMuPt)");
  desc.add<double>("maxPtForCompatibilityCheck", 25.0)
      ->setComment(
          "Maximum pT for which to apply the full compatibility check with L1TkMu. Above threshold keep all seeds in "
          "dR window.");

  // Quality
  desc.add<double>("maxCombinedQuality", 2.0)
      ->setComment("Maximum combined quality metric for seed to be considered (lower is better)");

  // Arbitration
  desc.add<int>("minNPixelHits", 2)
      ->setComment("Minimum number of pixel hits for seed to be prioritised in arbitration");
  desc.add<int>("minNConsistentHits", 4)
      ->setComment(
          "Minimum number of OT hits consistent with L1TkMuon TTTrack stubs for seed to be prioritised in arbitration");
  desc.add<uint32_t>("nSeedsToKeep", 2)
      ->setComment(
          "Maximum number of seeds to keep per L1TkMuon based on quality ranking (does not apply to L1TkMu above "
          "maxPtForCompatibilityCheck, for which all seeds passing pre-selection are kept)");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonSeedsSelectorFromL1TkMuon);
