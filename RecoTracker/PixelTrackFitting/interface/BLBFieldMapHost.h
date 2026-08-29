#ifndef RecoTracker_PixelTrackFitting_interface_BLBFieldMapHost_h
#define RecoTracker_PixelTrackFitting_interface_BLBFieldMapHost_h

#include <algorithm>
#include <array>

#include <alpaka/alpaka.hpp>

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "RecoTracker/PixelTrackFitting/interface/BLBFieldMap.h"

// Host-resident normalized (Bz,Br) r-z lattice (layout in BLBFieldMap.h). Sampled once per IOV from the
// MagneticField EventSetup product by BLBFieldMapESProducerAlpaka; copied to the device once per IOV
// (CopyToDevice<BLBFieldMapHost> in interface/alpaka/BLBFieldMapCollection.h).
class BLBFieldMapHost {
public:
  using Buffer = cms::alpakatools::host_buffer<float[]>;
  using ConstBuffer = cms::alpakatools::const_host_buffer<float[]>;

  explicit BLBFieldMapHost(std::array<float, blBFieldMap::kNValues> const& values)
      : buffer_(cms::alpakatools::make_host_buffer<float[]>(blBFieldMap::kNValues)) {
    std::copy_n(values.data(), blBFieldMap::kNValues, buffer_.data());
  }

  // non-copyable
  BLBFieldMapHost(BLBFieldMapHost const&) = delete;
  BLBFieldMapHost& operator=(BLBFieldMapHost const&) = delete;

  // movable (required by the ESProducer ownership / CopyToDevice machinery)
  BLBFieldMapHost(BLBFieldMapHost&&) = default;
  BLBFieldMapHost& operator=(BLBFieldMapHost&&) = default;

  ~BLBFieldMapHost() = default;

  // access the buffer (the CopyToDevice source)
  ConstBuffer buffer() const { return buffer_; }

  float const* data() const { return buffer_.data(); }

private:
  Buffer buffer_;
};

#endif  // RecoTracker_PixelTrackFitting_interface_BLBFieldMapHost_h
