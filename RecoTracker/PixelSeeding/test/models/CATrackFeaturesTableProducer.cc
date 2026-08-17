// CATrackFeaturesTableProducer - nano flat table with the in-kernel track DNN's 12-feature
// vector (CATrackDNN ABI) plus extra hit/stub columns, computed HOST-SIDE from the final track
// SoA hit lists via the shared CATrackFeatures.h -- the same single-source-of-truth function
// Kernel_classifyTracks evaluates on device. Pair with Trk<X>Truth (join offline by phi).
// Rows: ALL tracks in SoA order (same row space as Trk<X>Full); <3 hits / corrupt hit lists
// get NaN features (drop offline).
//
// The four DEPLOYED extras (rzChi2 / meanStubKappa / leverArm / rMax = PixelTrackFeaturesSoA
// columns 28-31) come from fill()'s rzKappaOut out-parameter, so they are bit-identical to the HP
// selector's feature kernel. Only meanClusterY / nTilted / tiltedFrac still come from a separate
// host study walk (skips tagged OT extras, which index the OT SoA and are never stubs).
//
// OPTIONAL MERGED-COLLECTION PROVENANCE (emitMergedProvenance=True, OFF by default): four extra
// track-level columns iteration / ndof / nOTExtra / nAttached. The merged collection is the only
// place iteration is meaningful (neither this table nor PixelTrackSoATab carries it otherwise).
//
// OPTIONAL PIXEL-CLUSTER AGGREGATES (emitClusterFeatures=True, OFF by default): eight extra
// track-level columns from the SiPixelCluster charge/shape on the track's PIXEL hits in the
// merged rechit SoA. Every input is a DataFormats/TrackingRecHitSoA column the device kernel has
// in hand, so the block is device-portable as written.
//
// OPTIONAL ATTACH-PURITY TABLE (emitHitTruth=True, OFF by default). A SECOND FlatTable, one row
// PER (track, hit) in SoA hit-list order, carrying per-hit truth so the purity of the in-CA fit
// extension's attached OT rechits can be measured against the track's own TrackingParticle. The
// track's own TP is the plurality TP over the track's resolvable (stub + OT-extra) hits. Emit
// ownTpNShared so a shared-hit threshold can be applied offline.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <vector>

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackSoA/interface/TracksHost.h"
#include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
#include "DataFormats/TrackingRecHitSoA/interface/OTRecHitsHost.h"
#include "DataFormats/TrackingRecHitSoA/interface/OTRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/StubsHost.h"
#include "DataFormats/TrackingRecHitSoA/interface/StubsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsHost.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "RecoTracker/PixelSeeding/interface/CATrackFeatures.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "RecoTracker/PixelSeeding/test/models/TrackerLayerId.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

class CATrackFeaturesTableProducer : public edm::global::EDProducer<> {
public:
  explicit CATrackFeaturesTableProducer(const edm::ParameterSet& params)
      : tableName_(params.getParameter<std::string>("tableName")),
        tracksToken_(consumes<reco::TracksHost>(params.getParameter<edm::InputTag>("trackSrc"))),
        hitsToken_(consumes<reco::TrackingRecHitHost>(params.getParameter<edm::InputTag>("mergedHitsSrc"))),
        otHitsToken_(consumes<reco::OTRecHitsHost>(params.getParameter<edm::InputTag>("otRecHitsSoASrc"))),
        emitMergedProvenance_(params.getParameter<bool>("emitMergedProvenance")),
        emitClusterFeatures_(params.getParameter<bool>("emitClusterFeatures")),
        pixelBarrelModuleEnd_(params.getParameter<unsigned int>("pixelBarrelModuleEnd")),
        lowChargeThreshold_(float(params.getParameter<double>("lowChargeThreshold"))),
        emitHitTruth_(params.getParameter<bool>("emitHitTruth")),
        hitTableName_(params.getParameter<std::string>("hitTableName")),
        minSharedForOwnTP_(params.getParameter<int>("minSharedForOwnTP")) {
    produces<nanoaod::FlatTable>(tableName_);
    if (emitHitTruth_) {
      stubsToken_ = consumes<reco::StubsHost>(params.getParameter<edm::InputTag>("stubsSrc"));
      otRecHitCollToken_ =
          consumes<Phase2TrackerRecHit1DCollectionNew>(params.getParameter<edm::InputTag>("otRecHitSrc"));
      tpToken_ = consumes<std::vector<TrackingParticle>>(params.getParameter<edm::InputTag>("trackingParticleSrc"));
      topoToken_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();
      geomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();
      hitAssocConfig_.emplace(params.getParameter<edm::ParameterSet>("hitAssociatorConfig"), consumesCollector());
      tpChargedOnly_ = params.getParameter<bool>("TP_chargedOnly");
      tpSignalOnly_ = params.getParameter<bool>("TP_signalOnly");
      tpMinPt_ = params.getParameter<double>("TP_minPt");
      tpMaxEta_ = params.getParameter<double>("TP_maxEta");
      tpMaxTip_ = params.getParameter<double>("TP_maxTip");
      tpMaxVtxZ_ = params.getParameter<double>("TP_maxVtxZ");
      produces<nanoaod::FlatTable>(hitTableName_);
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("tableName", "TrkDispCA");
    desc.add<edm::InputTag>("trackSrc", edm::InputTag("hltPhase2PixelTracksSoALowPt"));
    desc.add<edm::InputTag>("mergedHitsSrc", edm::InputTag("hltPhase2PixelRecHitsStubsMerger"));
    // Raw OT-rechit SoA: resolves tagged OT extras attached by the in-CA fit extension so the
    // deployed 12-feature vector matches the device Kernel_classifyTracks on OT-extended tracks.
    desc.add<edm::InputTag>("otRecHitsSoASrc", edm::InputTag("hltPixelSeedingOTRecHitsSoA"));

    // Merged-collection provenance columns (iteration/ndof/nOTExtra/nAttached). OFF by default so
    // the per-arm datasets already produced stay byte-identical; the merged dump turns it on, where
    // iteration is the only place the two arms are still distinguishable (TracksSoA.h).
    desc.add<bool>("emitMergedProvenance", false)
        ->setComment("Also emit per-track provenance columns (iteration/ndof/nOTExtra/nAttached)");

    // Pixel-CLUSTER shape/charge aggregates (OFF by default). Everything they use lives in
    // DataFormats/TrackingRecHitSoA (TrackingRecHitsSoA.h: chargeAndStatus() -- a
    // SiPixelHitStatusAndCharge with a 24-bit `charge` field in ELECTRONS plus the CPE status
    // bitfield -- clusterSizeX(), clusterSizeY(), detectorIndex()), so a device kernel can rebuild
    // them from the same merged-hit SoA the feature kernel already reads.
    desc.add<bool>("emitClusterFeatures", false)
        ->setComment("Also emit per-track pixel-cluster charge/shape aggregates (no vertex information)");
    // Pixel BARREL/ENDCAP split, needed only for the path-length normalisation of the charge.
    // Modules indexed in DetId order: pixel barrel occupies [0, 864) for the Phase-2 pixel
    // topology (4 TBPX layers). Exposed as a parameter so a different geometry needs no rebuild.
    desc.add<unsigned int>("pixelBarrelModuleEnd", 864)
        ->setComment("First pixel ENDCAP module index (detectorIndex < this == barrel)");
    desc.add<double>("lowChargeThreshold", 7000.0)
        ->setComment("Path-length-normalised cluster charge (electrons) below which a pixel hit counts as low-charge");

    // --- Optional attach-purity per-hit truth table (OFF by default; test-config only). ---
    desc.add<bool>("emitHitTruth", false)
        ->setComment("Also emit a per-(track,hit) truth table (isOTExtra/isStub/isTrueForOwnTP) for the purity study");
    desc.add<std::string>("hitTableName", "TrkDispCAHit");
    desc.add<int>("minSharedForOwnTP", 3)
        ->setComment("(diagnostic only) # shared resolvable hits for the own-TP plurality to count as matched");
    desc.add<edm::InputTag>("stubsSrc", edm::InputTag("hltOTStubProducer"))
        ->setComment("StubsHost (lowerHitIdx/upperHitIdx -> OTRecHitsSoA indices)");
    desc.add<edm::InputTag>("otRecHitSrc", edm::InputTag("hltSiPhase2RecHits"))
        ->setComment("Phase2TrackerRecHit1DCollectionNew for OT truth matching");
    desc.add<edm::InputTag>("trackingParticleSrc", edm::InputTag("mix", "MergedTrackTruth"));
    desc.add<bool>("TP_signalOnly", false);
    desc.add<bool>("TP_chargedOnly", true);
    desc.add<double>("TP_minPt", 0.0);
    desc.add<double>("TP_maxEta", 4.5);
    desc.add<double>("TP_maxTip", 60.0);
    desc.add<double>("TP_maxVtxZ", 60.0);
    edm::ParameterSetDescription hitAssocDesc;
    hitAssocDesc.add<bool>("associatePixel", true);
    hitAssocDesc.add<bool>("associateStrip", true);
    hitAssocDesc.add<bool>("usePhase2Tracker", true);
    hitAssocDesc.add<bool>("associateRecoTracks", false);
    hitAssocDesc.add<bool>("associateHitbySimTrack", false);
    hitAssocDesc.add<edm::InputTag>("phase2TrackerSimLinkSrc", edm::InputTag("simSiPixelDigis", "Tracker"));
    hitAssocDesc.add<edm::InputTag>("pixelSimLinkSrc", edm::InputTag("simSiPixelDigis", "Pixel"));
    hitAssocDesc.add<edm::InputTag>("stripSimLinkSrc", edm::InputTag("simSiStripDigis"));
    hitAssocDesc.add<std::vector<std::string>>("ROUList",
                                               {"TrackerHitsPixelBarrelLowTof",
                                                "TrackerHitsPixelBarrelHighTof",
                                                "TrackerHitsPixelEndcapLowTof",
                                                "TrackerHitsPixelEndcapHighTof",
                                                "TrackerHitsTIBLowTof",
                                                "TrackerHitsTIBHighTof",
                                                "TrackerHitsTIDLowTof",
                                                "TrackerHitsTIDHighTof",
                                                "TrackerHitsTOBLowTof",
                                                "TrackerHitsTOBHighTof",
                                                "TrackerHitsTECLowTof",
                                                "TrackerHitsTECHighTof"});
    desc.add<edm::ParameterSetDescription>("hitAssociatorConfig", hitAssocDesc);

    descriptions.addWithDefaultLabel(desc);
  }

private:
  void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const override {
    const auto& tracksHost = iEvent.get(tracksToken_);
    const auto& mergedHits = iEvent.get(hitsToken_);
    const auto& otHits = iEvent.get(otHitsToken_);
    const auto tracks = tracksHost.const_view().tracks();
    const auto trackHits = tracksHost.const_view().trackHits();
    const auto hh = mergedHits.const_view().trackingHits();
    const auto otView = otHits.const_view().otRecHits();
    const int nHitsTot = hh.metadata().size();
    const uint32_t nOTHits = otView.metadata().size();
    const ::reco::OTRecHitsConstView* otViewPtr = (nOTHits > 0u) ? &otView : nullptr;
    const int nTracks = tracks.nTracks();

    // Column names == caTrackFeatures::fill order == train_disp_nano.py FEATS (12 feats).
    static const char* kNames[caTrackFeatures::kNFeat] = {
        "fitChi2", "psFrac", "r0", "nPS", "nh", "spanZ", "nStubs", "nl", "logChi2Stub", "kErr", "dcaEst", "nBarrel"};

    constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
    std::vector<std::vector<float>> cols(caTrackFeatures::kNFeat, std::vector<float>(nTracks, kNaN));
    std::vector<float> phiCol(nTracks, kNaN);  // join key vs the truth table
    std::vector<int> qualityCol(nTracks, -1);  // SoA quality enum value

    // --- Columns 0-3 are DEPLOYED (PixelTrackFeaturesSoA columns 28-31, read by the HP selector);
    // columns 4-6 are host-only candidates: the CA features capture stub-INTERNAL consistency
    // (logChi2Stub, kErr), these add extrapolation/shape. Order == kExtraNames. nTilted/tiltedFrac
    // count barrel stubs on TILTED modules (StubFlags !isFlat), whose biased bend produces
    // near-beamline fakes at |eta|~1-2.
    static const char* kExtraNames[7] = {
        "meanStubKappa", "leverArm", "rMax", "rzChi2", "meanClusterY", "nTilted", "tiltedFrac"};
    std::vector<std::vector<float>> ex(7, std::vector<float>(nTracks, kNaN));

    // Merged-collection provenance (emitMergedProvenance_ only): the SoA columns the two table
    // producers otherwise drop. -1 = the track's hit span is empty/corrupt (no row content at all).
    std::vector<int> iterCol, ndofCol, nOTExtraCol, nAttachedCol;
    if (emitMergedProvenance_) {
      iterCol.assign(nTracks, -1);
      ndofCol.assign(nTracks, -1);
      nOTExtraCol.assign(nTracks, -1);
      nAttachedCol.assign(nTracks, -1);
    }

    // Pixel-cluster charge/shape aggregates (emitClusterFeatures_ only). Computed over the track's
    // PIXEL hits only: stub rows and OT extras carry no cluster information (the stubs merger and
    // the OT converter both zero charge/clusterSizeX/clusterSizeY explicitly). -1 on every column
    // when the track has no pixel hit with cluster information at all.
    // Order == kClusterNames.
    static const char* kClusterNames[8] = {
        "nPixHits", "minCharge", "meanCharge", "minChargeNorm", "maxSizeY", "meanSizeY", "maxSizeX", "nLowCharge"};
    std::vector<std::vector<float>> cl;
    if (emitClusterFeatures_)
      cl.assign(8, std::vector<float>(nTracks, -1.f));
    const uint32_t offsetStubsMain = hh.offsetStubs();

    float feat[caTrackFeatures::kNFeat];
    for (int it = 0; it < nTracks; ++it) {
      const uint32_t start = (it == 0) ? 0 : tracks[it - 1].hitOffsets();
      const uint32_t end = tracks[it].hitOffsets();
      if (end <= start || end > uint32_t(trackHits.metadata().size()))
        continue;
      const auto* hitsBegin = trackHits.id().data() + start;
      const auto* hitsEnd = trackHits.id().data() + end;

      // Provenance is independent of whether the feature walk succeeds, so it is filled on the
      // whole valid-span row space (an all-NaN feature row still says which iteration produced it).
      if (emitMergedProvenance_) {
        iterCol[it] = int(tracks[it].iteration());
        ndofCol[it] = int(tracks[it].ndof());
        int nOT = 0, nAtt = 0;
        for (uint32_t k = start; k < end; ++k) {
          if (caOTHitTag::isOTId(trackHits[k].id()))
            ++nOT;  // raw-OT extra attached by the extension walk
          if (trackHits[k].attached() == 1)
            ++nAtt;  // any hit added by the extension stage (superset of the OT extras)
        }
        nOTExtraCol[it] = nOT;
        nAttachedCol[it] = nAtt;
      }

      // --- Pixel-cluster aggregates. Filled on the whole valid-span row space (like provenance),
      // so a row whose feature walk fails still carries them. Path-length normalisation: state()(3)
      // = cotan(theta). A barrel sensor's normal is radial -> path ~ t/|sin(theta)|, so the normalised
      // charge is Q*|sin(theta)|; an endcap sensor's normal is along z -> Q*|cos(theta)|. Same as the
      // offline pixel cluster-charge cuts.
      if (emitClusterFeatures_) {
        const float cot = tracks[it].state()(3);
        const float invHyp = 1.f / std::sqrt(1.f + cot * cot);
        const float absSinTheta = invHyp;                  // barrel path factor
        const float absCosTheta = std::abs(cot) * invHyp;  // endcap path factor
        int nPix = 0, nLow = 0;
        float qMin = 0.f, qSum = 0.f, qnMin = 0.f;
        float syMax = 0.f, sySum = 0.f, sxMax = 0.f;
        for (uint32_t k = start; k < end; ++k) {
          const uint32_t h = trackHits[k].id();
          if (caOTHitTag::isOTId(h))
            continue;  // raw OT extra: indexes the OT SoA, no cluster information there
          if (h >= uint32_t(nHitsTot) || h >= offsetStubsMain)
            continue;  // stub row (charge/sizes zeroed by the merger) or corrupt index
          const float q = float(hh[h].chargeAndStatus().charge);
          if (!(q > 0.f))
            continue;  // defensive: a pixel row with no charge carries no usable cluster
          const bool barrel = uint32_t(hh[h].detectorIndex()) < pixelBarrelModuleEnd_;
          const float qn = q * (barrel ? absSinTheta : absCosTheta);
          const float sx = float(hh[h].clusterSizeX());
          const float sy = float(hh[h].clusterSizeY());
          if (nPix == 0) {
            qMin = q;
            qnMin = qn;
          } else {
            qMin = std::min(qMin, q);
            qnMin = std::min(qnMin, qn);
          }
          qSum += q;
          sySum += sy;
          syMax = std::max(syMax, sy);
          sxMax = std::max(sxMax, sx);
          if (qn < lowChargeThreshold_)
            ++nLow;
          ++nPix;
        }
        if (nPix > 0) {
          cl[0][it] = float(nPix);
          cl[1][it] = qMin;
          cl[2][it] = qSum / float(nPix);
          cl[3][it] = qnMin;
          cl[4][it] = syMax;
          cl[5][it] = sySum / float(nPix);
          cl[6][it] = sxMax;
          cl[7][it] = float(nLow);
        }
      }

      // The four DEPLOYED extras (rzChi2/meanStubKappa/leverArm/rMax) come from fill()'s rzKappaOut,
      // NOT from the host study walk below: the study walk skips bit30-tagged OT extras, while the
      // deployed kernel includes them, so a host recomputation would be a train-vs-serve skew on
      // every OT-extended track (exactly the tail a merged-collection model is meant to judge).
      // Sentinels match the kernel (rzChi2 = -1 when undefined, not NaN).
      float rzk[4] = {-1.f, 0.f, 0.f, 0.f};
      const bool ok = caTrackFeatures::fill(
          hitsBegin, hitsEnd, hh, nHitsTot, float(tracks[it].nLayers()), tracks[it].chi2(), feat, rzk, otViewPtr);
      if (!ok)
        continue;
      for (int f = 0; f < caTrackFeatures::kNFeat; ++f)
        cols[f][it] = feat[f];
      phiCol[it] = tracks[it].state()(0);
      qualityCol[it] = int(tracks[it].quality());
      ex[0][it] = rzk[1];  // meanStubKappa
      ex[1][it] = rzk[2];  // leverArm
      ex[2][it] = rzk[3];  // rMax
      ex[3][it] = rzk[0];  // rzChi2 (-1 = undefined, as the kernel writes it)

      // Study-only host walk over the same hit list. Tagged OT extras index the OT SoA, not hh, and
      // are never stubs -> skip them here (do NOT break, which would truncate at the first OT hit).
      int nClY = 0, nTilted = 0;
      double sumClY = 0.0;
      for (const uint32_t* ph = hitsBegin; ph != hitsEnd; ++ph) {
        const uint32_t h = *ph;
        if (caOTHitTag::isOTId(h))
          continue;
        if (h >= uint32_t(nHitsTot))
          break;
        const short cly = hh[h].clusterSizeY();
        if (cly > 0) {
          sumClY += cly;
          ++nClY;
        }
        if (::reco::isStub(hh, h)) {
          const auto flags = hh[h].stubFlags();
          if (::reco::StubFlags::isBarrel(flags) && !::reco::StubFlags::isFlat(flags))
            ++nTilted;  // tilted barrel module (|eta|~1-2 transition)
        }
      }
      ex[4][it] = (nClY > 0) ? float(sumClY / nClY) : 0.f;  // meanClusterY
      ex[5][it] = float(nTilted);                           // nTilted
      ex[6][it] = nTilted / std::max(feat[6], 1.f);         // tiltedFrac (feat[6]=nStubs)
    }

    auto table = std::make_unique<nanoaod::FlatTable>(nTracks, tableName_, /*singleton*/ false, /*extension*/ false);
    for (int f = 0; f < caTrackFeatures::kNFeat; ++f)
      table->addColumn<float>(kNames[f], cols[f], "CA track classifier feature", -1);
    for (int f = 0; f < 7; ++f)
      table->addColumn<float>(kExtraNames[f],
                              ex[f],
                              f < 4 ? "deployed CA track feature (PixelTrackFeaturesSoA columns 28-31, from "
                                      "caTrackFeatures::fill's rzKappaOut -- OT-aware, bit-identical to the kernel)"
                                    : "candidate track feature under study (extrapolation/shape, host-only)",
                              -1);
    table->addColumn<float>("phi", phiCol, "track phi (join key vs the truth table)", -1);
    table->addColumn<int>("quality", qualityCol, "SoA quality enum value", -1);
    if (emitMergedProvenance_) {
      table->addColumn<int>(
          "iteration", iterCol, "pixelTrack::Iteration that produced the track (-1 = empty/corrupt hit span)", -1);
      table->addColumn<int>("ndof", ndofCol, "fitted degrees of freedom (0 = never fitted, -1 = no hit span)", -1);
      table->addColumn<int>("nOTExtra", nOTExtraCol, "# raw-OT extras on the track (caOTHitTag-tagged hit ids)", -1);
      table->addColumn<int>(
          "nAttached", nAttachedCol, "# hits added by the extension stage (trackHits.attached()==1)", -1);
    }
    if (emitClusterFeatures_) {
      static const char* kClusterDocs[8] = {
          "# pixel hits on the track carrying cluster information (charge > 0)",
          "min SiPixelCluster charge over the track's pixel hits (electrons)",
          "mean SiPixelCluster charge over the track's pixel hits (electrons)",
          "min path-length-normalised cluster charge (Q*|sin(theta)| barrel, Q*|cos(theta)| endcap, electrons)",
          "max cluster sizeY over the track's pixel hits",
          "mean cluster sizeY over the track's pixel hits",
          "max cluster sizeX over the track's pixel hits",
          "# pixel hits whose normalised charge is below lowChargeThreshold"};
      for (int f = 0; f < 8; ++f)
        table->addColumn<float>(kClusterNames[f], cl[f], kClusterDocs[f], -1);
    }
    iEvent.put(std::move(table), tableName_);

    if (emitHitTruth_)
      produceHitTruthTable(iEvent, iSetup, tracks, trackHits, hh, otView, nOTHits, nTracks);
  }

  // ---- Attach-purity per-(track,hit) truth table -----------------------------------------------
  void produceHitTruthTable(edm::Event& iEvent,
                            const edm::EventSetup& iSetup,
                            const ::reco::TrackSoAConstView& tracks,
                            const ::reco::TrackHitSoAConstView& trackHits,
                            const ::reco::TrackingRecHitConstView& hh,
                            const ::reco::OTRecHitsConstView& otView,
                            uint32_t nOTHits,
                            int nTracks) const {
    const auto& tTopo = iSetup.getData(topoToken_);
    const auto& tGeom = iSetup.getData(geomToken_);
    const auto& stubsHost = iEvent.get(stubsToken_);
    const auto stubsView = stubsHost.const_view().stubs();
    const uint32_t nStubs = stubsView.metadata().size();
    const uint32_t offsetStubs = hh.offsetStubs();
    const auto& otRecHitColl = iEvent.get(otRecHitCollToken_);
    edm::Handle<std::vector<TrackingParticle>> tpH;
    iEvent.getByToken(tpToken_, tpH);

    TrackerHitAssociator hitAssociator(iEvent, *hitAssocConfig_);

    // (simTrackId, eventId) -> TP index
    std::map<SimHitIdpr, uint32_t> simTrackToTP;
    for (uint32_t iTP = 0; iTP < tpH->size(); ++iTP)
      for (const auto& g4 : (*tpH)[iTP].g4Tracks())
        simTrackToTP[{g4.trackId(), g4.eventId()}] = iTP;

    // TP selection (broad; the purity join wants any real association not to be dropped).
    std::vector<char> tpSel(tpH->size(), 0);
    for (uint32_t iTP = 0; iTP < tpH->size(); ++iTP) {
      const auto& tp = (*tpH)[iTP];
      if (tpChargedOnly_ && tp.charge() == 0)
        continue;
      if (tpSignalOnly_ && (tp.eventId().bunchCrossing() != 0 || tp.eventId().event() != 0))
        continue;
      if (tp.pt() < tpMinPt_ || std::abs(tp.eta()) > tpMaxEta_ || std::abs(tp.d0()) > tpMaxTip_ ||
          std::abs(tp.vz()) > tpMaxVtxZ_)
        continue;
      tpSel[iTP] = 1;
    }

    // flat Phase2 rechit index -> {set<selected TP index>, layerId} (same flat iteration the
    // OTRecHitsSoA converter used to assign origRecHitIdx).
    std::unordered_map<uint32_t, std::set<uint32_t>> origIdxToTPs;
    std::unordered_map<uint32_t, int> origIdxToLayer;
    {
      uint32_t flatIdx = 0;
      for (const auto& detSet : otRecHitColl) {
        const int layer =
            int(tracknano::getLayerId<pixelTopology::Phase2, uint16_t>(DetId(detSet.detId()), &tTopo, &tGeom));
        for (const auto& recHit : detSet) {
          origIdxToLayer[flatIdx] = layer;
          std::vector<SimHitIdpr> ids;
          hitAssociator.associatePhase2TrackerRecHit(&recHit, ids);
          for (const auto& id : ids) {
            auto it = simTrackToTP.find(id);
            if (it != simTrackToTP.end() && tpSel[it->second])
              origIdxToTPs[flatIdx].insert(it->second);
          }
          ++flatIdx;
        }
      }
    }
    auto soaTPs = [&](uint32_t soaIdx) -> const std::set<uint32_t>* {
      if (soaIdx >= nOTHits)
        return nullptr;
      auto it = origIdxToTPs.find(otView[soaIdx].origRecHitIdx());
      return (it != origIdxToTPs.end()) ? &it->second : nullptr;
    };
    auto soaLayer = [&](uint32_t soaIdx) -> int {
      if (soaIdx >= nOTHits)
        return -1;
      auto it = origIdxToLayer.find(otView[soaIdx].origRecHitIdx());
      return (it != origIdxToLayer.end()) ? it->second : -1;
    };

    // Per-(track,hit) rows.
    const uint32_t* hitIds = trackHits.id().data();
    std::vector<int> cTrackIdx, cLayer, cIsStub, cIsOTExtra, cIsTrue, cOwnTp, cOwnN, cHasTP, cHitId, cHitNTP, cHitTpKey;
    for (int it = 0; it < nTracks; ++it) {
      const uint32_t start = (it == 0) ? 0 : tracks[it - 1].hitOffsets();
      const uint32_t end = tracks[it].hitOffsets();
      if (end <= start || end > uint32_t(trackHits.metadata().size()))
        continue;

      // Resolve every hit's TP set + layer + kind first (one short list per track).
      const uint32_t nh = end - start;
      std::vector<std::set<uint32_t>> hitTPs(nh);
      std::vector<int> hitLayer(nh, -1), hitIsStub(nh, 0), hitIsOT(nh, 0), hitHasTP(nh, 0), hitNTP(nh, 0),
          hitTpKey(nh, -1);
      std::map<uint32_t, int> votes;  // TP index -> # of the track's resolvable hits sharing it
      for (uint32_t k = 0; k < nh; ++k) {
        const uint32_t h = hitIds[start + k];
        std::set<uint32_t> tps;
        if (caOTHitTag::isOTId(h)) {
          hitIsOT[k] = 1;
          const uint32_t o = caOTHitTag::otIdx(h);
          hitLayer[k] = soaLayer(o);
          if (const auto* p = soaTPs(o))
            tps = *p;
        } else if (h >= offsetStubs) {
          hitIsStub[k] = 1;
          const uint32_t stubIdx = h - offsetStubs;
          if (stubIdx < nStubs) {
            const uint32_t lowerIdx = stubsView[stubIdx].lowerHitIdx();
            const uint32_t upperIdx = stubsView[stubIdx].upperHitIdx();
            hitLayer[k] = soaLayer(lowerIdx);
            const auto* lo = soaTPs(lowerIdx);
            const auto* up = soaTPs(upperIdx);
            const bool paired = ::reco::isStub(stubsView, stubIdx) && upperIdx < nOTHits;
            if (lo) {
              if (paired) {
                if (up)
                  std::set_intersection(
                      lo->begin(), lo->end(), up->begin(), up->end(), std::inserter(tps, tps.begin()));
              } else {
                tps = *lo;  // PHitOnly: lower sensor only
              }
            }
          }
        }
        // pixel hits (h < offsetStubs, not tagged): left unresolved (tps empty, layer -1).
        if (!tps.empty()) {
          hitHasTP[k] = 1;
          for (uint32_t tp : tps)
            ++votes[tp];
        }
        hitNTP[k] = int(tps.size());  // raw associator TP-set size, BEFORE the plurality vote (shared-hit multiplicity)
        // Per-HIT TP key: the primary (smallest-index) selected TP this individual hit resolves to,
        // independent of the track's own-TP plurality. Lets an offline tool join hit -> TP -> layerId
        // across ALL tracks, e.g. to reconstruct each TP's true OT layer set. hitNTP>1 hits are ambiguous.
        hitTpKey[k] = tps.empty() ? -1 : int(*tps.begin());
        hitTPs[k] = std::move(tps);
      }

      // own TP = plurality over resolvable hits (ties -> smallest index, since std::map is ordered).
      int ownTp = -1, ownN = 0;
      for (const auto& kv : votes)
        if (kv.second > ownN) {
          ownN = kv.second;
          ownTp = int(kv.first);
        }

      for (uint32_t k = 0; k < nh; ++k) {
        cTrackIdx.push_back(it);
        cHitId.push_back(int(hitIds[start + k]));
        cLayer.push_back(hitLayer[k]);
        cIsStub.push_back(hitIsStub[k]);
        cIsOTExtra.push_back(hitIsOT[k]);
        cOwnTp.push_back(ownTp);
        cOwnN.push_back(ownN);
        cHasTP.push_back(hitHasTP[k]);
        cHitNTP.push_back(hitNTP[k]);
        cHitTpKey.push_back(hitTpKey[k]);
        // isTrueForOwnTP: only for resolvable (OT-extra / stub) hits on a RELIABLY matched track
        // (own-TP plurality shared by >= minSharedForOwnTP_ resolvable hits). ownTpNShared is also
        // emitted so a different threshold can be applied offline.
        int isTrue = -1;
        if ((hitIsOT[k] || hitIsStub[k]) && ownTp >= 0 && ownN >= minSharedForOwnTP_)
          isTrue = hitTPs[k].count(uint32_t(ownTp)) ? 1 : 0;
        cIsTrue.push_back(isTrue);
      }
    }

    const unsigned nRows = cTrackIdx.size();
    auto t = std::make_unique<nanoaod::FlatTable>(nRows, hitTableName_, /*singleton*/ false, /*extension*/ false);
    t->addColumn<int>("trackIdx", cTrackIdx, "SoA track index this hit belongs to", -1);
    t->addColumn<int>(
        "hitId", cHitId, "raw merged track-hit id (pixel/stub index into mergedHits; OT extras tagged bit30)", -1);
    t->addColumn<int>("layerId", cLayer, "getLayerId (Phase2 V1) OT layer 28-53, -1 if pixel/unresolved", -1);
    t->addColumn<int>("isStub", cIsStub, "1 if merged stub-region hit", -1);
    t->addColumn<int>("isOTExtra", cIsOTExtra, "1 if attached raw-OT extra (tag bit set)", -1);
    t->addColumn<int>("isTrueForOwnTP", cIsTrue, "1 hit shares track's own TP, 0 not, -1 unresolved/no ownTP", -1);
    t->addColumn<int>("ownTpKey", cOwnTp, "plurality TP index over the track's resolvable hits (-1 none)", -1);
    t->addColumn<int>("ownTpNShared", cOwnN, "# resolvable hits sharing the plurality TP", -1);
    t->addColumn<int>("hitHasTP", cHasTP, "1 if this hit resolved to >=1 selected TP", -1);
    t->addColumn<int>(
        "hitNTP",
        cHitNTP,
        "# distinct selected TrackingParticles the associator matched to this hit, BEFORE the "
        "track's plurality vote (0 = unresolved/pixel or no TP; used for the shared-hit multiplicity study)",
        -1);
    t->addColumn<int>("hitTpKey",
                      cHitTpKey,
                      "primary (smallest-index) selected TP this hit resolves to, independent of the track's ownTpKey "
                      "(-1 = unresolved/pixel or no TP). Join hit->TP->layerId across tracks for the Phase-0 per-layer "
                      "budget/captured/lost decomposition; hitNTP>1 flags an ambiguous (multi-TP) hit",
                      -1);
    iEvent.put(std::move(t), hitTableName_);
  }

  const std::string tableName_;
  const edm::EDGetTokenT<reco::TracksHost> tracksToken_;
  const edm::EDGetTokenT<reco::TrackingRecHitHost> hitsToken_;
  const edm::EDGetTokenT<reco::OTRecHitsHost> otHitsToken_;

  // Optional merged-collection provenance columns on the main table.
  const bool emitMergedProvenance_;

  // Optional pixel-cluster charge/shape aggregates on the main table.
  const bool emitClusterFeatures_;
  const unsigned int pixelBarrelModuleEnd_;
  const float lowChargeThreshold_;

  // Optional attach-purity table members (used only when emitHitTruth_).
  const bool emitHitTruth_;
  const std::string hitTableName_;
  const int minSharedForOwnTP_;
  edm::EDGetTokenT<reco::StubsHost> stubsToken_;
  edm::EDGetTokenT<Phase2TrackerRecHit1DCollectionNew> otRecHitCollToken_;
  edm::EDGetTokenT<std::vector<TrackingParticle>> tpToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  std::optional<TrackerHitAssociator::Config> hitAssocConfig_;
  bool tpChargedOnly_ = true;
  bool tpSignalOnly_ = false;
  double tpMinPt_ = 0.0;
  double tpMaxEta_ = 4.5;
  double tpMaxTip_ = 60.0;
  double tpMaxVtxZ_ = 60.0;
};

DEFINE_FWK_MODULE(CATrackFeaturesTableProducer);
