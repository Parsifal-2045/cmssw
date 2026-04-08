#include <alpaka/alpaka.hpp>
#include <cstdint>

#include "FWCore/Utilities/interface/EDMException.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"

namespace cms::alpakatools {
  void throwSharedMemoryLimitExceeded(const size_t nElements,
                                      const uint32_t nBlocks,
                                      const size_t requiredSharedMem,
                                      const size_t sharedMemLimit) {
    throw cms::Exception("SharedMemoryLimitExceeded")
        << "OneToManyAssoc: Shared memory limit exceeded for prefix scan of " << nElements << " elements in " << nBlocks
        << " blocks. Required shared memory: " << requiredSharedMem << " bytes. Shared memory limit: " << sharedMemLimit
        << " bytes.";
  }

  void throwIterativePrefixScanMaxLevelsExceeded(const uint32_t nLevels) {
    throw cms::Exception("IterativePrefixScanMaxLevelsExceeded")
        << "OneToManyAssoc: Requested an iterative prefix scan with " << nLevels
        << " levels, which is above the compile time constant " << cms::alpakatools::iterativePrefixScanMaxLevels
        << ". Consider increasing the value of iterativePrefixScanMaxLevels or reducing the problem's size.";
  }

}  // namespace cms::alpakatools
