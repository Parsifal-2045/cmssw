#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

class MuonPixelTracksDNNSelector : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>> {
public:
  explicit MuonPixelTracksDNNSelector(const edm::ParameterSet&, const cms::Ort::ONNXRuntime*);
  ~MuonPixelTracksDNNSelector() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);
  static std::unique_ptr<cms::Ort::ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const cms::Ort::ONNXRuntime*);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // Feature extraction
  std::vector<float> extractFeatures(const reco::Track& track, const l1t::TrackerMuonCollection& l1TkMuons) const;

  // Input tokens
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<l1t::TrackerMuonCollection> l1TkMuonsToken_;

  // Scaler parameters (loaded from config)
  std::vector<int> prunedFeatureIndices_;

  // Model parameters
  const float decisionThreshold_;
  const bool usePrunedFeatures_;
  const int nFeatures_;
};

MuonPixelTracksDNNSelector::MuonPixelTracksDNNSelector(const edm::ParameterSet& iConfig,
                                                       const cms::Ort::ONNXRuntime* cache)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      l1TkMuonsToken_(consumes<l1t::TrackerMuonCollection>(iConfig.getParameter<edm::InputTag>("l1TkMuons"))),
      decisionThreshold_(iConfig.getParameter<double>("decisionThreshold")),
      usePrunedFeatures_(iConfig.getParameter<bool>("usePrunedFeatures")),
      nFeatures_(iConfig.getParameter<int>("nFeatures")) {
  if (usePrunedFeatures_) {
    auto indices = iConfig.getParameter<std::vector<int>>("prunedFeatureIndices");
    prunedFeatureIndices_ = indices;
  }

  produces<reco::TrackCollection>();
  produces<std::vector<float>>("scores");
}

std::unique_ptr<cms::Ort::ONNXRuntime> MuonPixelTracksDNNSelector::initializeGlobalCache(
    const edm::ParameterSet& iConfig) {
  edm::FileInPath modelPath(iConfig.getParameter<std::string>("modelPath"));
  return std::make_unique<cms::Ort::ONNXRuntime>(modelPath.fullPath());
}

void MuonPixelTracksDNNSelector::globalEndJob(const cms::Ort::ONNXRuntime* cache) {}

std::vector<float> MuonPixelTracksDNNSelector::extractFeatures(const reco::Track& track,
                                                               const l1t::TrackerMuonCollection& l1TkMuons) const {
  constexpr float kEpsilon = 1e-6f;

  // Extract basic track features (before log transform)
  float p = track.p();
  float pt = track.pt();
  float ptErr = track.ptError();
  float eta = track.eta();
  float etaErr = track.etaError();
  float phi = track.phi();
  float phiErr = track.phiError();
  float chi2 = track.chi2();
  float normalizedChi2 = track.normalizedChi2();

  // Hit information
  int nPixelHits = track.hitPattern().numberOfValidPixelHits();
  int nTrkLays = track.hitPattern().trackerLayersWithMeasurement();
  int nFoundHits = track.numberOfValidHits();
  int nLostHits = track.numberOfLostHits();

  // Track parameters and errors
  float dsz = track.dsz();
  float dszErr = track.dszError();
  float dxy = track.dxy();
  float dxyErr = track.dxyError();
  float dz = track.dz();
  float dzErr = track.dzError();
  float qoverp = track.qoverp();
  float qoverpErr = track.qoverpError();
  float lambdaErr = track.lambdaError();

  // Build full feature vector
  std::vector<float> features;
  features.reserve(32);  // Total number of features before pruning

  // Original features with log transforms
  features.push_back(std::log10(p + kEpsilon));               // 0: p
  features.push_back(std::log10(pt + kEpsilon));              // 1: pt
  features.push_back(std::log10(ptErr + kEpsilon));           // 2: ptErr
  features.push_back(eta);                                    // 3: eta
  features.push_back(std::log10(etaErr + kEpsilon));          // 4: etaErr
  features.push_back(phi);                                    // 5: phi
  features.push_back(std::log10(phiErr + kEpsilon));          // 6: phiErr
  features.push_back(std::log10(chi2 + kEpsilon));            // 7: chi2
  features.push_back(std::log10(normalizedChi2 + kEpsilon));  // 8: normalizedChi2
  features.push_back(nPixelHits);                             // 9: nPixelHits
  features.push_back(nTrkLays);                               // 10: nTrkLays
  features.push_back(nFoundHits);                             // 11: nFoundHits
  features.push_back(nLostHits);                              // 12: nLostHits
  features.push_back(dsz);                                    // 13: dsz
  features.push_back(std::log10(dszErr + kEpsilon));          // 14: dszErr
  features.push_back(dxy);                                    // 15: dxy
  features.push_back(std::log10(dxyErr + kEpsilon));          // 16: dxyErr
  features.push_back(dz);                                     // 17: dz
  features.push_back(std::log10(dzErr + kEpsilon));           // 18: dzErr
  features.push_back(qoverp);                                 // 19: qoverp
  features.push_back(std::log10(qoverpErr + kEpsilon));       // 20: qoverpErr
  features.push_back(std::log10(lambdaErr + kEpsilon));       // 21: lambdaErr

  // Derived features
  float sigmaPtOverPt = ptErr / pt;
  features.push_back(std::log10(sigmaPtOverPt + kEpsilon));  // 22: sigmaPtOverPt

  float hitEfficiency = static_cast<float>(nFoundHits) / std::max(nFoundHits + nLostHits, 1);
  features.push_back(hitEfficiency);  // 23: hitEfficiency

  float chi2PerHit = chi2 / std::max(nFoundHits, 1);
  features.push_back(std::log10(chi2PerHit + kEpsilon));  // 24: chi2PerHit

  float impact3D = dxy * dxy + dz * dz;
  features.push_back(std::log10(impact3D + kEpsilon));  // 25: impact3D

  float dxySignificance = dxy / dxyErr;
  float dzSignificance = dz / dzErr;
  float impactSignificance = std::sqrt(dxySignificance * dxySignificance + dzSignificance * dzSignificance);
  features.push_back(std::log10(impactSignificance + kEpsilon));  // 26: impactSignificance

  float relUncertaintyProduct = sigmaPtOverPt * (qoverpErr / std::abs(qoverp));
  features.push_back(std::log10(relUncertaintyProduct + kEpsilon));  // 27: relUncertaintyProduct

  // L1TkMuon matching features
  float minDR2 = 999.0f;
  float l1PtMatched = -1.0f;

  for (auto it = l1TkMuons.begin(); it != l1TkMuons.end(); ++it) {
    auto l1TTrack = it->trkPtr();
    float dr2 = reco::deltaR2(eta, phi, l1TTrack->momentum().eta(), l1TTrack->momentum().phi());
    if (dr2 < minDR2) {
      minDR2 = dr2;
      l1PtMatched = l1TTrack->momentum().perp();
    }
  }

  float dPtNorm = (l1PtMatched > 0) ? std::abs(pt - l1PtMatched) / l1PtMatched : 999.0f;
  float ptRatio = (l1PtMatched > 0) ? pt / l1PtMatched : 999.0f;
  float matchingScore = minDR2 * (1.0f + dPtNorm);

  features.push_back(std::log10(minDR2 + kEpsilon));         // 28: L1TkMu_dR2min
  features.push_back(std::log10(dPtNorm + kEpsilon));        // 29: L1TkMu_dPtNorm
  features.push_back(std::log10(ptRatio + kEpsilon));        // 30: L1TkMu_ptRatio
  features.push_back(std::log10(matchingScore + kEpsilon));  // 31: L1TkMu_matchingScore

  if (usePrunedFeatures_) {
    std::vector<float> prunedFeatures;
    prunedFeatures.reserve(prunedFeatureIndices_.size());
    for (int idx : prunedFeatureIndices_) {
      prunedFeatures.push_back(features[idx]);
    }
    return prunedFeatures;
  }

  return features;
}

void MuonPixelTracksDNNSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto selectedTracks = std::make_unique<reco::TrackCollection>();
  auto scores = std::make_unique<std::vector<float>>();

  // Get input collections
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<l1t::TrackerMuonCollection> l1TkMuons;
  iEvent.getByToken(l1TkMuonsToken_, l1TkMuons);

  if (!tracks.isValid() || tracks->empty()) {
    iEvent.put(std::move(selectedTracks));
    iEvent.put(std::move(scores), "scores");
    return;
  }

  // Prepare batch input
  std::vector<float> inputData;
  inputData.reserve(tracks->size() * nFeatures_);

  for (const auto& track : *tracks) {
    auto features = extractFeatures(track, *l1TkMuons);
    inputData.insert(inputData.end(), features.begin(), features.end());
  }

  // Run inference

  std::vector<std::vector<int64_t>> inputShapes = {
      {static_cast<int64_t>(tracks->size()), static_cast<int64_t>(nFeatures_)}};
  cms::Ort::FloatArrays inputTensor({inputData});

  auto outputs = globalCache()->run({"input"}, inputTensor, inputShapes, {"output"}, 1);

  // Apply sigmoid and threshold
  const auto& logits = outputs[0];
  for (size_t i = 0; i < tracks->size(); ++i) {
    float logit = logits[i];
    float prob = 1.0f / (1.0f + std::exp(-logit));
    scores->push_back(prob);

    if (prob >= decisionThreshold_) {
      selectedTracks->push_back((*tracks)[i]);
    }
  }

  std::cout << "MuonPixelTracksDNNSelector - Selected " << selectedTracks->size() << " out of " << tracks->size()
            << " tracks." << std::endl;

  iEvent.put(std::move(selectedTracks));
  iEvent.put(std::move(scores), "scores");
}

void MuonPixelTracksDNNSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("hltPhase2L3FromL1TkMuonPixelTracks"));
  desc.add<edm::InputTag>("l1TkMuons", edm::InputTag("l1tTkMuonsGmt"));
  desc.add<std::string>("modelPath", "RecoMuon/L3TrackFinder/data/muon_pixeltrack_selector.onnx");
  desc.add<std::vector<int>>("prunedFeatureIndices", {});
  desc.add<double>("decisionThreshold", 0.5);
  desc.add<bool>("usePrunedFeatures", true);
  desc.add<int>("nFeatures", 0);
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonPixelTracksDNNSelector);
