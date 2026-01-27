#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "DataFormats/L1Trigger/interface/EGamma.h"
typedef BXVectorSimpleFlatTableProducer<l1t::EGamma> SimpleTriggerL1EGFlatTableProducer;

#include "DataFormats/L1Trigger/interface/Jet.h"
typedef BXVectorSimpleFlatTableProducer<l1t::Jet> SimpleTriggerL1JetFlatTableProducer;

#include "DataFormats/L1Trigger/interface/Tau.h"
typedef BXVectorSimpleFlatTableProducer<l1t::Tau> SimpleTriggerL1TauFlatTableProducer;

#include "DataFormats/L1Trigger/interface/Muon.h"
typedef BXVectorSimpleFlatTableProducer<l1t::Muon> SimpleTriggerL1MuonFlatTableProducer;

#include "DataFormats/L1Trigger/interface/EtSum.h"
typedef BXVectorSimpleFlatTableProducer<l1t::EtSum> SimpleTriggerL1EtSumFlatTableProducer;

#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
typedef SimpleFlatTableProducer<l1t::TrackerMuon> SimpleTriggerL1TkMuonFlatTableProducer;

class L1TkMuStubFlatTableProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TkMuStubFlatTableProducer(const edm::ParameterSet &iConfig)
      : src_(consumes<l1t::TrackerMuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
        name_(iConfig.getParameter<std::string>("name")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) override {
    // Input collection
    edm::Handle<l1t::TrackerMuonCollection> muons;
    iEvent.getByToken(src_, muons);

    // Output features
    std::vector<int> muIdxs;
    std::vector<int> etaRegions, phiRegions, depthRegions, types, qualities;

    if (muons.isValid() && !muons->empty()) {
      muIdxs.reserve(muons->size() * 4);
      etaRegions.reserve(muons->size() * 4);
      phiRegions.reserve(muons->size() * 4);
      depthRegions.reserve(muons->size() * 4);
      types.reserve(muons->size() * 4);
      qualities.reserve(muons->size() * 4);

      for (unsigned int iMu = 0; iMu < muons->size(); ++iMu) {
        const auto &mu = muons->at(iMu);
        //std::cout << "\nProcessing muon " << iMu << " with " << mu.stubs().size() << " stubs." << std::endl;
        //std::cout << "L1TkMuon Pt: " << mu.phPt() << ", Eta: " << mu.phEta() << ", Phi: " << mu.phPhi() << std::endl;
        for (const auto &stubRef : mu.stubs()) {
          if (stubRef.isNull())
            continue;

          const auto &stub = *stubRef;

          //std::cout << "Stub details - Parent L1TkMuon: " << iMu << " \nEtaRegion: " << stub.etaRegion()
          //          << ", PhiRegion: " << stub.phiRegion() << ", DepthRegion: " << stub.depthRegion()
          //          << ", Type: " << stub.type() << ", Quality: " << stub.quality() << std::endl;
          //std::cout << "Stub offline coordinates - offlineCoord1: " << stub.offline_coord1()
          //          << ", offlineCoord2: " << stub.offline_coord2() << ", offline eta1: " << stub.offline_eta1()
          //          << ", offline eta2: " << stub.offline_eta2() << std::endl;
          //std::cout << "---------------------------------\n";

          muIdxs.push_back(iMu);
          etaRegions.push_back(stub.etaRegion());
          phiRegions.push_back(stub.phiRegion());
          depthRegions.push_back(stub.depthRegion());
          types.push_back(stub.type());
          qualities.push_back(stub.quality());
        }
        //std::cout << "=============================================\n";
      }
    }

    auto out = std::make_unique<nanoaod::FlatTable>(etaRegions.size(), name_, false);
    out->setDoc("Stubs belonging to L1TkMuons");
    out->addColumn<int>("parentL1TkMu", muIdxs, "Index of the parent L1TkMuon");
    out->addColumn<int>("etaRegion", etaRegions, "Eta region (Wheel/Ring)");
    out->addColumn<int>("phiRegion", phiRegions, "Phi region (Sector/Chamber)");
    out->addColumn<int>("depthRegion", depthRegions, "Station/Depth");
    out->addColumn<int>("type", types, "Stub type: 0=DT, 1=RPC Barrel, 2=CSC, 3=RPC Endcap");
    out->addColumn<int>("quality", qualities, "Stub quality");

    iEvent.put(std::move(out));
  }

private:
  edm::EDGetTokenT<l1t::TrackerMuonCollection> src_;
  std::string name_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SimpleTriggerL1EGFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTriggerL1JetFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTriggerL1MuonFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTriggerL1TauFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTriggerL1EtSumFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTriggerL1TkMuonFlatTableProducer);
DEFINE_FWK_MODULE(L1TkMuStubFlatTableProducer);
