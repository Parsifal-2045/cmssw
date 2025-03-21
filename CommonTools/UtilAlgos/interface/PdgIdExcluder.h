#ifndef CommonTools_UtilAlgos_PdgIdSelector_h
#define CommonTools_UtilAlgos_PdgIdSelector_h
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/Utils/interface/PdgIdExcluder.h"

namespace reco {
  namespace modules {

    template <>
    struct ParameterAdapter<PdgIdExcluder> {
      static PdgIdExcluder make(const edm::ParameterSet& cfg, edm::ConsumesCollector& iC) {
        return PdgIdExcluder(cfg.getParameter<std::vector<int> >("pdgId"));
      }

      static void fillPSetDescription(edm::ParameterSetDescription& desc) { desc.add<std::vector<int> >("pdgId", {}); }
    };

  }  // namespace modules
}  // namespace reco

#endif
