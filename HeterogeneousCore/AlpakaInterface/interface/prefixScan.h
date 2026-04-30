#ifndef HeterogeneousCore_AlpakaInterface_interface_prefixScan_h
#define HeterogeneousCore_AlpakaInterface_interface_prefixScan_h

#include <alpaka/alpaka.hpp>

#include "FWCore/Utilities/interface/CMSUnrollLoop.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
namespace cms::alpakatools {
  template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
  constexpr bool isPowerOf2(T v) {
    // returns true iif v has only one bit set.
    while (v) {
      if (v & 1)
        return !(v >> 1);
      else
        v >>= 1;
    }
    return false;
  }

  // Decoupled lookback tile descriptor helpers
  // A tile descriptor is a 64-bit word packing a 32-bit status in the upper half
  // and a 32-bit value (aggregate or inclusive prefix, cast from T) in the lower half.
  // Atomic 64-bit CAS/Exch on the packed word guarantees the reader always sees a
  // consistent (status, value) pair, avoiding the race where a value could be
  // overwritten between a status read and a value read.
  // Restricted to integral T with sizeof(T) <= 4 (uint16_t, int32_t, uint32_t, etc.).

  constexpr uint32_t kTileStatusInvalid = 0u;
  constexpr uint32_t kTileStatusAggregate = 1u;
  constexpr uint32_t kTileStatusPrefix = 2u;

  template <typename T>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE uint64_t packTileDescriptor(uint32_t status, T value) {
    static_assert(std::is_integral_v<T> && sizeof(T) <= sizeof(uint32_t),
                  "Packed tile descriptors require integral types of at most 32 bits");
    return (static_cast<uint64_t>(status) << 32) | static_cast<uint64_t>(static_cast<uint32_t>(value));
  }

  template <typename T>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE T unpackTileValue(uint64_t packed) {
    return static_cast<T>(static_cast<uint32_t>(packed));
  }

  ALPAKA_FN_ACC ALPAKA_FN_INLINE uint32_t unpackTileStatus(uint64_t packed) {
    return static_cast<uint32_t>(packed >> 32);
  }

  // Portable count-trailing-zeros for ballot masks in the warp-cooperative lookback.
  template <typename TWord>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE int ctzBallot(TWord mask) {
    int n = 0;
    while (mask && !(mask & TWord(1))) {
      mask >>= 1;
      ++n;
    }
    return n;
  }

  template <alpaka::concepts::Acc TAcc, typename T>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void warpPrefixScan(
      const TAcc& acc, int32_t laneId, T const* ci, T* co, uint32_t i, bool active = true) {
    // ci and co may be the same
    T x = active ? ci[i] : 0;
    CMS_UNROLL_LOOP
    for (int32_t offset = 1; offset < alpaka::warp::getSize(acc); offset <<= 1) {
      // Force the exact type for integer types otherwise the compiler will find the template resolution ambiguous.
      using dataType = std::conditional_t<std::is_floating_point_v<T>, T, std::int32_t>;
      T y = alpaka::warp::shfl(acc, static_cast<dataType>(x), laneId - offset);
      if (laneId >= offset)
        x += y;
    }
    if (active)
      co[i] = x;
  }

  template <alpaka::concepts::Acc TAcc, typename T>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void warpPrefixScan(
      const TAcc& acc, int32_t laneId, T* c, uint32_t i, bool active = true) {
    warpPrefixScan(acc, laneId, c, c, i, active);
  }

  // limited to warpSize² elements
  template <alpaka::concepts::Acc TAcc, typename T>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void blockPrefixScan(
      const TAcc& acc, T const* ci, T* co, int32_t size, T* ws = nullptr) {
    if constexpr (!requires_single_thread_per_block_v<TAcc>) {
      const auto warpSize = alpaka::warp::getSize(acc);
      int32_t const blockDimension(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
      int32_t const blockThreadIdx(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
      ALPAKA_ASSERT_ACC(ws);
      ALPAKA_ASSERT_ACC(size <= warpSize * warpSize);
      ALPAKA_ASSERT_ACC(0 == blockDimension % warpSize);
      auto first = blockThreadIdx;
      ALPAKA_ASSERT_ACC(isPowerOf2(warpSize));
      auto laneId = blockThreadIdx & (warpSize - 1);
      auto warpUpRoundedSize = (size + warpSize - 1) / warpSize * warpSize;

      for (auto i = first; i < warpUpRoundedSize; i += blockDimension) {
        // When padding the warp, warpPrefixScan is a noop
        warpPrefixScan(acc, laneId, ci, co, i, i < size);
        if (i < size) {
          // Skipped in warp padding threads.
          auto warpId = i / warpSize;
          ALPAKA_ASSERT_ACC(warpId < warpSize);
          if ((warpSize - 1) == laneId)
            ws[warpId] = co[i];
        }
      }
      alpaka::syncBlockThreads(acc);
      if (size <= warpSize)
        return;
      if (blockThreadIdx < warpSize) {
        warpPrefixScan(acc, laneId, ws, blockThreadIdx);
      }
      alpaka::syncBlockThreads(acc);
      for (auto i = first + warpSize; i < size; i += blockDimension) {
        int32_t warpId = i / warpSize;
        co[i] += ws[warpId - 1];
      }
      alpaka::syncBlockThreads(acc);
    } else {
      co[0] = ci[0];
      for (int32_t i = 1; i < size; ++i)
        co[i] = ci[i] + co[i - 1];
    }
  }

  template <alpaka::concepts::Acc TAcc, typename T>
  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE void blockPrefixScan(const TAcc& acc,
                                                           T* __restrict__ c,
                                                           int32_t size,
                                                           T* __restrict__ ws = nullptr) {
    if constexpr (!requires_single_thread_per_block_v<TAcc>) {
      const auto warpSize = alpaka::warp::getSize(acc);
      int32_t const blockDimension(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
      int32_t const blockThreadIdx(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
      ALPAKA_ASSERT_ACC(ws);
      ALPAKA_ASSERT_ACC(size <= warpSize * warpSize);
      ALPAKA_ASSERT_ACC(0 == blockDimension % warpSize);
      auto first = blockThreadIdx;
      auto laneId = blockThreadIdx & (warpSize - 1);
      auto warpUpRoundedSize = (size + warpSize - 1) / warpSize * warpSize;

      for (auto i = first; i < warpUpRoundedSize; i += blockDimension) {
        // When padding the warp, warpPrefixScan is a noop
        warpPrefixScan(acc, laneId, c, i, i < size);
        if (i < size) {
          // Skipped in warp padding threads.
          auto warpId = i / warpSize;
          ALPAKA_ASSERT_ACC(warpId < warpSize);
          if ((warpSize - 1) == laneId)
            ws[warpId] = c[i];
        }
      }
      alpaka::syncBlockThreads(acc);
      if (size <= warpSize)
        return;
      if (blockThreadIdx < warpSize) {
        warpPrefixScan(acc, laneId, ws, blockThreadIdx);
      }
      alpaka::syncBlockThreads(acc);
      for (auto i = first + warpSize; i < size; i += blockDimension) {
        auto warpId = i / warpSize;
        c[i] += ws[warpId - 1];
      }
      alpaka::syncBlockThreads(acc);
    } else {
      for (int32_t i = 1; i < size; ++i)
        c[i] += c[i - 1];
    }
  }

  // Throws an exception with a message containing the shared memory requirements and limit.
  void throwSharedMemoryLimitExceeded(const size_t nElements,
                                      const uint32_t nBlocks,
                                      const size_t requiredSharedMem,
                                      const size_t sharedMemLimit);

  // Verify shared memory requirements
  template <alpaka::concepts::Acc TAcc, typename TSize>
  ALPAKA_FN_INLINE static void checkSharedMemoryPrefixScan(TSize nElements,
                                                           uint32_t nBlocks,
                                                           alpaka::Dev<TAcc> const& device) {
    auto requiredSharedMem = (nBlocks + alpaka::getPreferredWarpSize(device)) * sizeof(TSize);
    auto sharedMemLimit = alpaka::getAccDevProps<TAcc>(device).m_sharedMemSizeBytes;
    if (requiredSharedMem > sharedMemLimit) {
      throwSharedMemoryLimitExceeded(static_cast<size_t>(nElements), nBlocks, requiredSharedMem, sharedMemLimit);
    }
  }

  // in principle not limited.... in practice limited by shared memory size and occupancy.
  template <typename T>
  struct multiBlockPrefixScan {
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC void operator()(
        const TAcc& acc, T const* ci, T* co, uint32_t size, int32_t numBlocks, int32_t* pc, std::size_t warpSize) const {
      // Get shared variable. The workspace is needed only for multi-threaded accelerators.
      T* ws = nullptr;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        ws = alpaka::getDynSharedMem<T>(acc);
      }
      ALPAKA_ASSERT_ACC(warpSize == static_cast<std::size_t>(alpaka::warp::getSize(acc)));
      [[maybe_unused]] const auto elementsPerGrid = alpaka::getWorkDiv<alpaka::Grid, alpaka::Elems>(acc)[0u];
      const auto elementsPerBlock = alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[0u];
      const auto threadsPerBlock = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
      const auto blocksPerGrid = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      const auto blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
      ALPAKA_ASSERT_ACC(elementsPerGrid >= size);
      // first each block does a scan
      [[maybe_unused]] int off = elementsPerBlock * blockIdx;
      if (size - off > 0) {
        blockPrefixScan(acc, ci + off, co + off, alpaka::math::min(acc, elementsPerBlock, size - off), ws);
      }

      // count blocks that finished
      auto& isLastBlockDone = alpaka::declareSharedVar<bool, __COUNTER__>(acc);
      //__shared__ bool isLastBlockDone;
      if (0 == threadIdx) {
        alpaka::mem_fence(acc, alpaka::memory_scope::Device{});
        auto value = alpaka::atomicAdd(acc, pc, 1, alpaka::hierarchy::Blocks{});  // block counter
        isLastBlockDone = (value == (int(blocksPerGrid) - 1));
      }

      alpaka::syncBlockThreads(acc);

      if (!isLastBlockDone)
        return;

      ALPAKA_ASSERT_ACC(int(blocksPerGrid) == *pc);

      // good each block has done its work and now we are left in last block

      // let's get the partial sums from each block except the last, which receives 0.
      T* psum = nullptr;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        psum = ws + warpSize;
      } else {
        psum = alpaka::getDynSharedMem<T>(acc);
      }
      for (int32_t i = threadIdx, ni = blocksPerGrid; i < ni; i += threadsPerBlock) {
        auto j = elementsPerBlock * i + elementsPerBlock - 1;
        psum[i] = (j < size) ? co[j] : T(0);
      }
      alpaka::syncBlockThreads(acc);
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        if (blocksPerGrid <= warpSize * warpSize)
          blockPrefixScan(acc, psum, blocksPerGrid, ws);
        else {
          auto off = 0u;
          while (off + warpSize * warpSize < blocksPerGrid) {
            blockPrefixScan(acc, psum + off, warpSize * warpSize, ws);
            off = off + warpSize * warpSize - 1;
            // ^ this -1 is to keep the previous round total sum around
            alpaka::syncBlockThreads(acc);
          }
          blockPrefixScan(acc, psum + off, psum + off, blocksPerGrid - off, ws);
        }
      } else {
        blockPrefixScan(acc, psum, blocksPerGrid, ws);
      }
      // now it would have been handy to have the other blocks around...
      // Simplify the computation by having one version where threads per block = block size
      // and a second for the one thread per block accelerator.
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        //  Here threadsPerBlock == elementsPerBlock
        for (uint32_t i = threadIdx + threadsPerBlock, k = 0; i < size; i += threadsPerBlock, ++k) {
          co[i] += psum[k];
        }
      } else {
        // We are single threaded here, adding partial sums starting with the 2nd block.
        for (uint32_t i = elementsPerBlock; i < size; i++) {
          co[i] += psum[i / elementsPerBlock - 1];
        }
      }
    }
  };

  // Two kernel approach, not memory limited but more overhead due to kernel launch and global memory usage.
  // Kernel A: scan one level (tile per block) and emit one block sum per block.
  // Kernel B: add scanned block offsets to each block (except block 0).
  // Kernel A can be called recursively on smaller and smaller block-sums arrays orchestration happens on host

  // Kernel A
  template <typename T>
  struct scanTilesWriteBlockSums {
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC void operator()(
        const TAcc& acc, T const* ci, T* co, uint32_t size, T* blockSums, std::size_t warpSize) const {
      T* ws = nullptr;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        ws = alpaka::getDynSharedMem<T>(acc);
      }

      const auto elementsPerBlock = alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[0u];
      const auto blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];

      auto off = elementsPerBlock * blockIdx;
      if (off >= size)
        return;

      auto n = alpaka::math::min(acc, elementsPerBlock, size - off);
      blockPrefixScan(acc, ci + off, co + off, static_cast<int32_t>(n), ws);

      if (threadIdx == 0u) {
        blockSums[blockIdx] = co[off + n - 1];
      }
    }
  };

  // Fused Kernel A: scan tiles, then the last block to finish scans the block sums
  // and propagates offsets back into the output.
  // Removes the need for a separate single-block top level and its Kernel B pass.
  // Follows the same last-block-finishes pattern as multiBlockPrefixScan.
  template <typename T>
  struct scanTilesFused {
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                  T const* ci,
                                  T* co,
                                  uint32_t size,
                                  T* blockSums,
                                  int32_t numBlocks,
                                  int32_t* pc,
                                  std::size_t warpSize) const {
      T* ws = nullptr;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        ws = alpaka::getDynSharedMem<T>(acc);
      }
      ALPAKA_ASSERT_ACC(warpSize == static_cast<std::size_t>(alpaka::warp::getSize(acc)));

      const auto elementsPerBlock = alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[0u];
      const auto threadsPerBlock = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
      const auto blocksPerGrid = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      const auto blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];

      // Each block scans its tile (same as scanTilesWriteBlockSums)
      auto off = elementsPerBlock * blockIdx;
      if (off < size) {
        auto n = alpaka::math::min(acc, elementsPerBlock, size - off);
        blockPrefixScan(acc, ci + off, co + off, static_cast<int32_t>(n), ws);
        if (threadIdx == 0u) {
          blockSums[blockIdx] = co[off + n - 1];
        }
      }

      // Last block?
      auto& isLastBlockDone = alpaka::declareSharedVar<bool, __COUNTER__>(acc);
      if (0 == threadIdx) {
        alpaka::mem_fence(acc, alpaka::memory_scope::Device{});
        auto value = alpaka::atomicAdd(acc, pc, 1, alpaka::hierarchy::Blocks{});
        isLastBlockDone = (value == (int(blocksPerGrid) - 1));
      }
      alpaka::syncBlockThreads(acc);
      if (!isLastBlockDone)
        return;

      ALPAKA_ASSERT_ACC(int(blocksPerGrid) == *pc);

      // Last block scans all block sums and adds offsets
      T* psum = nullptr;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        psum = ws + warpSize;
      } else {
        psum = alpaka::getDynSharedMem<T>(acc);
      }

      for (uint32_t i = threadIdx; i < blocksPerGrid; i += threadsPerBlock) {
        psum[i] = blockSums[i];
      }
      alpaka::syncBlockThreads(acc);

      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        if (blocksPerGrid <= warpSize * warpSize)
          blockPrefixScan(acc, psum, static_cast<int32_t>(blocksPerGrid), ws);
        else {
          auto scanOff = 0u;
          while (scanOff + warpSize * warpSize < blocksPerGrid) {
            blockPrefixScan(acc, psum + scanOff, static_cast<int32_t>(warpSize * warpSize), ws);
            scanOff = scanOff + warpSize * warpSize - 1;
            alpaka::syncBlockThreads(acc);
          }
          blockPrefixScan(acc, psum + scanOff, psum + scanOff, static_cast<int32_t>(blocksPerGrid - scanOff), ws);
        }
      } else {
        blockPrefixScan(acc, psum, static_cast<int32_t>(blocksPerGrid), ws);
      }

      // Write scanned block sums back to global memory
      for (uint32_t i = threadIdx; i < blocksPerGrid; i += threadsPerBlock) {
        blockSums[i] = psum[i];
      }

      // Add scanned block offsets to output (skip first block, it needs no offset)
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        for (uint32_t i = threadIdx + threadsPerBlock, k = 0; i < size; i += threadsPerBlock, ++k) {
          co[i] += psum[k];
        }
      } else {
        for (uint32_t i = elementsPerBlock; i < size; i++) {
          co[i] += psum[i / elementsPerBlock - 1];
        }
      }
    }
  };

  // Kernel B
  template <typename T>
  struct addScannedBlockOffsets {
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC void operator()(const TAcc& acc, T* data, uint32_t size, T const* scannedBlockSums) const {
      const auto elementsPerBlock = alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[0u];
      const auto blockIdx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
      const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];

      if (blockIdx == 0u)
        return;

      T const add = scannedBlockSums[blockIdx - 1u];
      uint32_t begin = elementsPerBlock * blockIdx;
      uint32_t end = alpaka::math::min(acc, begin + elementsPerBlock, size);
      for (auto i = begin + threadIdx; i < end; i += elementsPerBlock) {
        data[i] += add;
      }
    }
  };

  // Decoupled lookback prefix scan (Merrill & Garland, 2016)
  // Single kernel launch with fully parallel offset propagation.  Each block scans
  // its tile, publishes a packed (status, value) descriptor via atomic 64-bit
  // exchange, and resolves its global prefix by looking back at predecessors.
  //
  // Key optimisations modelled on CUB's DeviceScan:
  //  1. ITEMS_PER_THREAD register-level scan to amortise shared-memory and
  //     lookback overhead (1024 threads x 4 items = 4096-element tiles).
  //  2. Warp-cooperative lookback: the first warp polls warpSize predecessors in
  //     parallel, combines them with a warp inclusive scan + ballot, and advances
  //     in window-sized steps until a PREFIX descriptor is found.
  //  3. Spin-loop backoff via a block-scope fence to prevent the compiler from
  //     collapsing the polling loop into an uncached tight spin.
  //
  // Requires sizeof(T) <= 4 for the packed 64-bit descriptor.
  // NOTE: the spin-loop relies on the GPU's forward progress guarantee for
  // resident blocks.  The grid is capped to device occupancy and each block
  // processes multiple tiles via a global counter, ensuring no deadlock.
  // The CPU serial backend takes the sequential fallback and never enters this.
  template <typename T, uint32_t ITEMS_PER_THREAD = 4u>
  struct scanDecoupledLookback {
    template <alpaka::concepts::Acc TAcc>
    ALPAKA_FN_ACC void operator()(const TAcc& acc,
                                  T const* ci,
                                  T* co,
                                  uint32_t size,
                                  uint64_t* tileDescriptors,
                                  int32_t* tileCounter,
                                  std::size_t warpSize) const {
      // Dynamic shared memory layout:
      //   ws[warpSize] - blockPrefixScan workspace
      //   threadAggs[threadsPerBlock] - per-thread aggregates for the block-wide scan
      T* ws = nullptr;
      T* threadAggs = nullptr;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        ws = alpaka::getDynSharedMem<T>(acc);
        threadAggs = ws + warpSize;
      }
      ALPAKA_ASSERT_ACC(warpSize == static_cast<std::size_t>(alpaka::warp::getSize(acc)));

      const auto threadsPerBlock = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
      const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
      uint32_t const elementsPerTile = threadsPerBlock * ITEMS_PER_THREAD;
      uint32_t const numTiles = (size + elementsPerTile - 1u) / elementsPerTile;
      using ShflT = std::conditional_t<std::is_floating_point_v<T>, T, std::int32_t>;

      // Each block processes tiles in order until all are consumed.
      // The grid is sized to the device occupancy, not to numTiles,
      // so blocks may loop multiple times for large inputs.
      while (true) {
        // Step 1: acquire a tile index
        auto& tileIdx = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
        if (threadIdx == 0u) {
          tileIdx = static_cast<uint32_t>(alpaka::atomicAdd(acc, tileCounter, 1, alpaka::hierarchy::Blocks{}));
        }
        alpaka::syncBlockThreads(acc);
        if (tileIdx >= numTiles)
          return;

        uint32_t tileOff = elementsPerTile * tileIdx;
        uint32_t tileEnd = alpaka::math::min(acc, tileOff + elementsPerTile, size);

        // Step 2: load ITEMS_PER_THREAD consecutive elements into registers
        T items[ITEMS_PER_THREAD];
        uint32_t threadOff = tileOff + threadIdx * ITEMS_PER_THREAD;
        CMS_UNROLL_LOOP
        for (uint32_t k = 0; k < ITEMS_PER_THREAD; ++k) {
          items[k] = (threadOff + k < tileEnd) ? ci[threadOff + k] : T(0);
        }

        // Step 3: sequential inclusive scan in registers (pure ALU, no sync)
        CMS_UNROLL_LOOP
        for (uint32_t k = 1; k < ITEMS_PER_THREAD; ++k) {
          items[k] += items[k - 1];
        }

        // Step 4: block-wide inclusive scan of per-thread aggregates
        // threadAggs[t] = items[ITEMS-1] for thread t
        threadAggs[threadIdx] = items[ITEMS_PER_THREAD - 1];
        alpaka::syncBlockThreads(acc);
        blockPrefixScan(acc, threadAggs, static_cast<int32_t>(threadsPerBlock), ws);
        // threadAggs[t] now holds the inclusive prefix sum of aggregates [0..t].

        // Apply the block-level exclusive prefix to the register items.
        T threadPrefix = (threadIdx > 0u) ? threadAggs[threadIdx - 1u] : T(0);
        CMS_UNROLL_LOOP
        for (uint32_t k = 0; k < ITEMS_PER_THREAD; ++k) {
          items[k] += threadPrefix;
        }
        // items[] now contains the tile-local inclusive prefix scan.

        // Step 5: publish descriptor and resolve global prefix
        T tileAggregate = threadAggs[threadsPerBlock - 1u];
        auto& exclusivePrefix = alpaka::declareSharedVar<T, __COUNTER__>(acc);

        if constexpr (!requires_single_thread_per_block_v<TAcc>) {
          // Only the first warp participates in the lookback.
          if (threadIdx < static_cast<uint32_t>(warpSize)) {
            int32_t laneId = static_cast<int32_t>(threadIdx);

            if (tileIdx == 0u) {
              // Tile 0 has no predecessors.
              if (laneId == 0) {
                alpaka::mem_fence(acc, alpaka::memory_scope::Device{});
                alpaka::atomicExch(acc,
                                   &tileDescriptors[0u],
                                   packTileDescriptor(kTileStatusPrefix, tileAggregate),
                                   alpaka::hierarchy::Blocks{});
                exclusivePrefix = T(0);
              }
            } else {
              // Publish local aggregate so successors can chain through us.
              if (laneId == 0) {
                alpaka::mem_fence(acc, alpaka::memory_scope::Device{});
                alpaka::atomicExch(acc,
                                   &tileDescriptors[tileIdx],
                                   packTileDescriptor(kTileStatusAggregate, tileAggregate),
                                   alpaka::hierarchy::Blocks{});
              }

              // Warp-cooperative lookback: poll warpSize predecessors per window.
              T runningPrefix = T(0);
              int32_t windowBase = static_cast<int32_t>(tileIdx) - 1;

              while (true) {
                int32_t predIdx = windowBase - laneId;

                // Every lane spins on its predecessor until valid.
                uint64_t pred = 0;
                bool valid = (predIdx < 0);
                while (!alpaka::warp::all(acc, static_cast<int>(valid))) {
                  if (!valid) {
                    pred = alpaka::atomicCas(
                        acc, &tileDescriptors[predIdx], uint64_t(0), uint64_t(0), alpaka::hierarchy::Blocks{});
                    valid = (unpackTileStatus(pred) != kTileStatusInvalid);
                  }
                  // Backoff: block-scope fence prevents the compiler from collapsing
                  // the spin into an uncached tight loop.
                  alpaka::mem_fence(acc, alpaka::memory_scope::Block{});
                }

                // Decode each lane's predecessor.
                T val = (predIdx >= 0) ? unpackTileValue<T>(pred) : T(0);
                uint32_t status = (predIdx >= 0) ? unpackTileStatus(pred) : kTileStatusPrefix;

                // Warp inclusive scan of values across the window.
                T scanned = val;
                CMS_UNROLL_LOOP
                for (int32_t offset = 1; offset < static_cast<int32_t>(warpSize); offset <<= 1) {
                  T y = static_cast<T>(alpaka::warp::shfl(acc, static_cast<ShflT>(scanned), laneId - offset));
                  if (laneId >= offset)
                    scanned += y;
                }

                // Ballot: identify lanes with PREFIX status.
                auto prefixBallot = alpaka::warp::ballot(acc, static_cast<int>(status == kTileStatusPrefix));

                if (prefixBallot) {
                  // Lowest set bit = most recent predecessor with a complete prefix.
                  // The warp scan result at that lane is the total contribution.
                  int firstPrefixLane = ctzBallot(prefixBallot);
                  T windowResult =
                      static_cast<T>(alpaka::warp::shfl(acc, static_cast<ShflT>(scanned), firstPrefixLane));
                  runningPrefix += windowResult;
                  break;
                }
                // No PREFIX in this window, accumulate total and slide back.
                T windowSum = static_cast<T>(
                    alpaka::warp::shfl(acc, static_cast<ShflT>(scanned), static_cast<int>(warpSize) - 1));
                runningPrefix += windowSum;
                windowBase -= static_cast<int32_t>(warpSize);
              }

              // Lane 0 publishes the result and promotes to PREFIX.
              if (laneId == 0) {
                exclusivePrefix = runningPrefix;
                alpaka::mem_fence(acc, alpaka::memory_scope::Device{});
                alpaka::atomicExch(acc,
                                   &tileDescriptors[tileIdx],
                                   packTileDescriptor(kTileStatusPrefix, exclusivePrefix + tileAggregate),
                                   alpaka::hierarchy::Blocks{});
              }
            }
          }
        } else {
          // Single-thread backend: trivial sequential lookback.
          if (tileIdx == 0u) {
            tileDescriptors[0] = packTileDescriptor(kTileStatusPrefix, tileAggregate);
            exclusivePrefix = T(0);
          } else {
            tileDescriptors[tileIdx] = packTileDescriptor(kTileStatusAggregate, tileAggregate);
            T running = T(0);
            for (int32_t lb = static_cast<int32_t>(tileIdx) - 1; lb >= 0; --lb) {
              running += unpackTileValue<T>(tileDescriptors[lb]);
              if (unpackTileStatus(tileDescriptors[lb]) == kTileStatusPrefix)
                break;
            }
            exclusivePrefix = running;
            tileDescriptors[tileIdx] = packTileDescriptor(kTileStatusPrefix, exclusivePrefix + tileAggregate);
          }
        }
        alpaka::syncBlockThreads(acc);

        // Step 6: add global exclusive prefix and store results
        if (exclusivePrefix != T(0)) {
          CMS_UNROLL_LOOP
          for (uint32_t k = 0; k < ITEMS_PER_THREAD; ++k) {
            items[k] += exclusivePrefix;
          }
        }
        CMS_UNROLL_LOOP
        for (uint32_t k = 0; k < ITEMS_PER_THREAD; ++k) {
          if (threadOff + k < tileEnd)
            co[threadOff + k] = items[k];
        }
        alpaka::syncBlockThreads(acc);
      }
    }
  };

  // Helper struct and function to compute the number of levels and block sizes for the two-kernel prefix scan.
  // 3 levels (0, 1, 2) can cover up to 1024^3 ~ 10^9 elements, enough for any realistic use case
  constexpr uint32_t iterativePrefixScanMaxLevels = 3;
  struct PrefixScanLevelPlan {
    uint32_t nLevels = 0;
    std::vector<uint32_t> levelSize;
    std::vector<uint32_t> levelBlocks;
  };

  // Throws an exception with a message containing the requested iterations and the set compile-time limit
  void throwIterativePrefixScanMaxLevelsExceeded(const size_t nElements, const uint32_t nLevels);

  ALPAKA_FN_INLINE static PrefixScanLevelPlan makePrefixScanLevelPlan(uint32_t nOnes, uint32_t nthreads) {
    PrefixScanLevelPlan p;
    if (nOnes == 0u) {
      return p;
    }

    p.levelSize.reserve(iterativePrefixScanMaxLevels);
    p.levelBlocks.reserve(iterativePrefixScanMaxLevels);

    p.nLevels = 1u;
    p.levelSize.emplace_back(nOnes);
    p.levelBlocks.emplace_back((nOnes + nthreads - 1u) / nthreads);

    while (p.levelBlocks[p.nLevels - 1u] > 1u) {
      if (p.nLevels >= iterativePrefixScanMaxLevels) {
        throwIterativePrefixScanMaxLevelsExceeded(nOnes, p.nLevels);
      }
      p.levelSize.emplace_back(p.levelBlocks[p.nLevels - 1u]);
      p.levelBlocks.emplace_back((p.levelSize[p.nLevels] + nthreads - 1u) / nthreads);
      ++p.nLevels;
    }
    return p;
  }

  template <alpaka::concepts::Acc TAcc, typename TQueue, typename T>
  ALPAKA_FN_INLINE static void iterativePrefixScan(T* input, T* output, uint32_t size, TQueue& queue) {
    if (size == 0u) {
      return;
    }

    if constexpr (!requires_single_thread_per_block_v<TAcc>) {
      constexpr uint32_t nthreads = 1024u;
      auto const plan = makePrefixScanLevelPlan(size, nthreads);
      assert(plan.nLevels > 0);

      std::vector<cms::alpakatools::device_buffer<alpaka::Dev<TQueue>, T[]>> blockSumsBuffers;
      blockSumsBuffers.reserve(plan.nLevels);
      for (uint32_t l = 0; l < plan.nLevels; ++l) {
        blockSumsBuffers.emplace_back(cms::alpakatools::make_device_buffer<T[]>(queue, plan.levelBlocks[l]));
      }

      auto const warpSize = alpaka::getPreferredWarpSize(alpaka::getDev(queue));

      //std::cout << "Running iterative prefix scan for " << size << " elements with " << plan.nLevels << " levels\n";
      //for (uint32_t l = 0; l < plan.nLevels; ++l) {
      //  std::cout << "  Level " << l << ": size = " << plan.levelSize[l] << ", blocks = " << plan.levelBlocks[l]
      //            << "\n";
      //}

      // Kernel A on level-0 input data
      auto workDiv = cms::alpakatools::make_workdiv<TAcc>(plan.levelBlocks[0], nthreads);
      alpaka::exec<TAcc>(queue,
                         workDiv,
                         scanTilesWriteBlockSums<T>{},
                         input,
                         output,
                         plan.levelSize[0],
                         blockSumsBuffers[0].data(),
                         warpSize);

      // Iterative use of kernel A on block-sum levels
      for (uint32_t l = 1; l < plan.nLevels; ++l) {
        auto workDiv = cms::alpakatools::make_workdiv<TAcc>(plan.levelBlocks[l], nthreads);
        alpaka::exec<TAcc>(queue,
                           workDiv,
                           scanTilesWriteBlockSums<T>{},
                           blockSumsBuffers[l - 1].data(),
                           blockSumsBuffers[l - 1].data(),
                           plan.levelSize[l],
                           blockSumsBuffers[l].data(),
                           warpSize);
      }

      // Kernel B from top-1 down to level 0
      for (int32_t l = static_cast<int32_t>(plan.nLevels) - 2; l >= 0; --l) {
        auto workDiv = cms::alpakatools::make_workdiv<TAcc>(plan.levelBlocks[l], nthreads);
        T* levelData = (l == 0) ? output : blockSumsBuffers[l - 1].data();
        alpaka::exec<TAcc>(
            queue, workDiv, addScannedBlockOffsets<T>{}, levelData, plan.levelSize[l], blockSumsBuffers[l].data());
      }
    } else {
      output[0] = input[0];
      for (uint32_t i = 1; i < size; ++i) {
        output[i] = input[i] + output[i - 1];
      }
    }
  }

  template <alpaka::concepts::Acc TAcc, typename TQueue, typename T>
  ALPAKA_FN_INLINE static void iterativePrefixScanFused(T* input, T* output, uint32_t size, TQueue& queue) {
    if (size == 0u) {
      return;
    }

    if constexpr (!requires_single_thread_per_block_v<TAcc>) {
      constexpr uint32_t nthreads = 1024u;
      auto const plan = makePrefixScanLevelPlan(size, nthreads);
      assert(plan.nLevels > 0);

      std::vector<cms::alpakatools::device_buffer<alpaka::Dev<TQueue>, T[]>> blockSumsBuffers;
      blockSumsBuffers.reserve(plan.nLevels);
      for (uint32_t l = 0; l < plan.nLevels; ++l) {
        blockSumsBuffers.emplace_back(cms::alpakatools::make_device_buffer<T[]>(queue, plan.levelBlocks[l]));
      }

      auto const warpSize = alpaka::getPreferredWarpSize(alpaka::getDev(queue));

      // The last level has 1 block by construction. Fuse it into the previous level using the
      // last-block-finishes pattern, eliminating one kernel launch and its Kernel B pass.

      // fusedLevel: the topmost multi-block level that will use scanTilesFused.
      // For nLevels >= 2: fusedLevel = nLevels - 2 (folds the single-block top level).
      // For nLevels == 1: fusedLevel = 0, but levelBlocks[0] == 1 so the fused path is skipped.
      uint32_t fusedLevel = (plan.nLevels >= 2u) ? plan.nLevels - 2u : 0u;

      // Kernel A on levels below the fused level
      for (uint32_t l = 0; l < fusedLevel; ++l) {
        T* ci = (l == 0) ? input : blockSumsBuffers[l - 1].data();
        T* co = (l == 0) ? output : blockSumsBuffers[l - 1].data();
        auto workDiv = cms::alpakatools::make_workdiv<TAcc>(plan.levelBlocks[l], nthreads);
        alpaka::exec<TAcc>(queue,
                           workDiv,
                           scanTilesWriteBlockSums<T>{},
                           ci,
                           co,
                           plan.levelSize[l],
                           blockSumsBuffers[l].data(),
                           warpSize);
      }

      // Fused kernel on the topmost multi-block level
      if (plan.levelBlocks[fusedLevel] > 1u) {
        T* ci = (fusedLevel == 0) ? input : blockSumsBuffers[fusedLevel - 1].data();
        T* co = (fusedLevel == 0) ? output : blockSumsBuffers[fusedLevel - 1].data();
        auto pcBuf = cms::alpakatools::make_device_buffer<int32_t>(queue);
        alpaka::memset(queue, pcBuf, 0);
        auto workDiv = cms::alpakatools::make_workdiv<TAcc>(plan.levelBlocks[fusedLevel], nthreads);
        alpaka::exec<TAcc>(queue,
                           workDiv,
                           scanTilesFused<T>{},
                           ci,
                           co,
                           plan.levelSize[fusedLevel],
                           blockSumsBuffers[fusedLevel].data(),
                           static_cast<int32_t>(plan.levelBlocks[fusedLevel]),
                           pcBuf.data(),
                           warpSize);
      } else {
        // Single block at fused level: normal Kernel A suffices (no cross-block work needed)
        T* ci = (fusedLevel == 0) ? input : blockSumsBuffers[fusedLevel - 1].data();
        T* co = (fusedLevel == 0) ? output : blockSumsBuffers[fusedLevel - 1].data();
        auto workDiv = cms::alpakatools::make_workdiv<TAcc>(1u, nthreads);
        alpaka::exec<TAcc>(queue,
                           workDiv,
                           scanTilesWriteBlockSums<T>{},
                           ci,
                           co,
                           plan.levelSize[fusedLevel],
                           blockSumsBuffers[fusedLevel].data(),
                           warpSize);
      }

      // Kernel B from fusedLevel-1 down to level 0
      // (the fused kernel already scanned block sums and added offsets for fusedLevel)
      for (int32_t l = static_cast<int32_t>(fusedLevel) - 1; l >= 0; --l) {
        auto workDiv = cms::alpakatools::make_workdiv<TAcc>(plan.levelBlocks[l], nthreads);
        T* levelData = (l == 0) ? output : blockSumsBuffers[l - 1].data();
        alpaka::exec<TAcc>(
            queue, workDiv, addScannedBlockOffsets<T>{}, levelData, plan.levelSize[l], blockSumsBuffers[l].data());
      }
    } else {
      output[0] = input[0];
      for (uint32_t i = 1; i < size; ++i) {
        output[i] = input[i] + output[i - 1];
      }
    }
  }

  // Decoupled lookback prefix scan host-side interface

  // Returns the number of tiles (== required tileDescriptor capacity) for a
  // decoupled lookback prefix scan over nElements items.
  // Call at initialisation time with the maximum expected nElements to size
  // the pre-allocated descriptor buffer.
  ALPAKA_FN_INLINE static uint32_t decoupledLookbackNumTiles(uint32_t nElements) {
    // Must stay in sync with the nthreads/itemsPerThread used below.
    constexpr uint32_t elementsPerTile = 1024u * 4u;
    return (nElements + elementsPerTile - 1u) / elementsPerTile;
  }

  // Pre-allocated workspace for the decoupled lookback prefix scan.
  // Bundles the tile descriptor array and the tile counter into a single object
  // that can be allocated once at construction time and reused across events.
  //
  // Usage:
  //   // At construction (once, sized for the worst case):
  //   auto ws = DecoupledLookbackWorkspace<TDev>::create(maxNOnes, queue);
  //
  //   // Per event (hot path, zero allocation):
  //   decoupledLookbackPrefixScan<TAcc>(poff, poff, nOnes, queue, ws);
  //
  template <typename TDev>
  struct DecoupledLookbackWorkspace {
    cms::alpakatools::device_buffer<TDev, uint64_t[]> tileDescriptors;
    cms::alpakatools::device_buffer<TDev, int32_t> tileCounter;
    uint32_t capacity;  // number of tiles this workspace can hold

    // Factory: create a workspace large enough for a prefix scan of up to maxElements items.
    template <typename TQueue>
    static DecoupledLookbackWorkspace create(uint32_t maxElements, TQueue& queue) {
      uint32_t numTiles = decoupledLookbackNumTiles(maxElements);
      return DecoupledLookbackWorkspace{cms::alpakatools::make_device_buffer<uint64_t[]>(queue, numTiles),
                                        cms::alpakatools::make_device_buffer<int32_t>(queue),
                                        numTiles};
    }

    // Zero-initialise before each use.  Enqueues two async memsets.
    template <typename TQueue>
    void reset(TQueue& queue) {
      alpaka::memset(queue, tileDescriptors, 0);
      alpaka::memset(queue, tileCounter, 0);
    }
  };

  // Decoupled lookback prefix scan pre-allocated overload.
  // The caller provides device-accessible buffers for the tile descriptors and
  // the tile counter, avoiding any per-call device memory allocation.
  //
  // tileDescriptors  must have capacity >= decoupledLookbackNumTiles(size).
  // tileCounter      must point to a single int32_t in device memory.
  //
  // IMPORTANT: both tileDescriptors (first numTiles entries) and *tileCounter
  // MUST be zero-initialised before each call.  Use alpaka::memset on the
  // owning device_buffer objects.
  template <alpaka::concepts::Acc TAcc, typename TQueue, typename T>
  ALPAKA_FN_INLINE static void decoupledLookbackPrefixScan(T* input,
                                                           T* output,
                                                           uint32_t size,
                                                           TQueue& queue,
                                                           uint64_t* tileDescriptors,
                                                           uint32_t descriptorCapacity,
                                                           int32_t* tileCounter) {
    if (size == 0u)
      return;

    if constexpr (!requires_single_thread_per_block_v<TAcc>) {
      constexpr uint32_t nthreads = 1024u;
      constexpr uint32_t itemsPerThread = 4u;
      constexpr uint32_t elementsPerTile = nthreads * itemsPerThread;
      uint32_t const numTiles = (size + elementsPerTile - 1u) / elementsPerTile;
      assert(descriptorCapacity >= numTiles);

      auto const& device = alpaka::getDev(queue);
      auto const warpSize = alpaka::getPreferredWarpSize(device);

      // Cap grid to device occupancy for the forward progress guarantee.
      auto const devProps = alpaka::getAccDevProps<TAcc>(device);
      uint32_t const maxResidentBlocks = static_cast<uint32_t>(devProps.m_multiProcessorCount) * 2u;
      uint32_t const nBlocks = std::min(numTiles, maxResidentBlocks);

      auto workDiv = cms::alpakatools::make_workdiv<TAcc>(nBlocks, nthreads);
      alpaka::exec<TAcc>(queue,
                         workDiv,
                         scanDecoupledLookback<T, itemsPerThread>{},
                         input,
                         output,
                         size,
                         tileDescriptors,
                         tileCounter,
                         warpSize);
    } else {
      output[0] = input[0];
      for (uint32_t i = 1; i < size; ++i)
        output[i] = input[i] + output[i - 1];
    }
  }

  // Decoupled lookback prefix scan — workspace overload.
  // Zero-initialises the workspace and delegates to the pre-allocated overload.
  // No per-call device memory allocation.
  template <alpaka::concepts::Acc TAcc, typename TQueue, typename T>
  ALPAKA_FN_INLINE static void decoupledLookbackPrefixScan(
      T* input, T* output, uint32_t size, TQueue& queue, DecoupledLookbackWorkspace<alpaka::Dev<TQueue>>& workspace) {
    if (size == 0u)
      return;

    if constexpr (!requires_single_thread_per_block_v<TAcc>) {
      workspace.reset(queue);
      decoupledLookbackPrefixScan<TAcc>(input,
                                        output,
                                        size,
                                        queue,
                                        workspace.tileDescriptors.data(),
                                        workspace.capacity,
                                        workspace.tileCounter.data());
    } else {
      output[0] = input[0];
      for (uint32_t i = 1; i < size; ++i)
        output[i] = input[i] + output[i - 1];
    }
  }

  // Decoupled lookback prefix scan self-allocating convenience overload.
  // Allocates tile descriptors and counter internally.  Suitable for one-off calls;
  // for production use with many concurrent streams, prefer the pre-allocated
  // overload to avoid allocator contention on the hot path.
  template <alpaka::concepts::Acc TAcc, typename TQueue, typename T>
  ALPAKA_FN_INLINE static void decoupledLookbackPrefixScan(T* input, T* output, uint32_t size, TQueue& queue) {
    if (size == 0u)
      return;

    if constexpr (!requires_single_thread_per_block_v<TAcc>) {
      uint32_t const numTiles = decoupledLookbackNumTiles(size);

      auto tileDescriptorsBuf = cms::alpakatools::make_device_buffer<uint64_t[]>(queue, numTiles);
      auto tileCounterBuf = cms::alpakatools::make_device_buffer<int32_t>(queue);
      alpaka::memset(queue, tileDescriptorsBuf, 0);
      alpaka::memset(queue, tileCounterBuf, 0);

      decoupledLookbackPrefixScan<TAcc>(
          input, output, size, queue, tileDescriptorsBuf.data(), numTiles, tileCounterBuf.data());
    } else {
      output[0] = input[0];
      for (uint32_t i = 1; i < size; ++i)
        output[i] = input[i] + output[i - 1];
    }
  }

}  // namespace cms::alpakatools

// declare the amount of block shared memory used by the multiBlockPrefixScan kernel
namespace alpaka::trait {
  // Variable size shared mem
  template <alpaka::concepts::Acc TAcc, typename T>
  struct BlockSharedMemDynSizeBytes<cms::alpakatools::multiBlockPrefixScan<T>, TAcc> {
    template <typename TVec>
    ALPAKA_FN_HOST_ACC static std::size_t getBlockSharedMemDynSizeBytes(
        cms::alpakatools::multiBlockPrefixScan<T> const& /* kernel */,
        TVec const& /* blockThreadExtent */,
        TVec const& /* threadElemExtent */,
        T const* /* ci */,
        T const* /* co */,
        int32_t /* size */,
        int32_t numBlocks,
        int32_t const* /* pc */,
        // This trait function does not receive the accelerator object to look up the warp size
        std::size_t warpSize) {
      // We need workspace (T[warpsize]) + partial sums (T[numblocks]).
      if constexpr (cms::alpakatools::requires_single_thread_per_block_v<TAcc>) {
        return sizeof(T) * numBlocks;
      } else {
        return sizeof(T) * (warpSize + numBlocks);
      }
    }
  };

  // Two-kernel approach requires only workspace for the first kernel, which is sized to the warp size.
  // Kernel A
  template <alpaka::concepts::Acc TAcc, typename T>
  struct BlockSharedMemDynSizeBytes<cms::alpakatools::scanTilesWriteBlockSums<T>, TAcc> {
    template <typename TVec>
    ALPAKA_FN_HOST_ACC static std::size_t getBlockSharedMemDynSizeBytes(
        cms::alpakatools::scanTilesWriteBlockSums<T> const& /* kernel */,
        TVec const& /* blockThreadExtent */,
        TVec const& /* threadElemExtent */,
        T const* /* ci */,
        T const* /* co */,
        uint32_t /* size */,
        T* /* blockSums */,
        // This trait function does not receive the accelerator object to look up the warp size
        std::size_t warpSize) {
      if constexpr (cms::alpakatools::requires_single_thread_per_block_v<TAcc>) {
        return 0;
      } else {
        return sizeof(T) * warpSize;
      }
    }
  };

  // Fused Kernel A: workspace (T[warpsize]) + partial sums (T[numblocks]), same as multiBlockPrefixScan.
  template <alpaka::concepts::Acc TAcc, typename T>
  struct BlockSharedMemDynSizeBytes<cms::alpakatools::scanTilesFused<T>, TAcc> {
    template <typename TVec>
    ALPAKA_FN_HOST_ACC static std::size_t getBlockSharedMemDynSizeBytes(
        cms::alpakatools::scanTilesFused<T> const& /* kernel */,
        TVec const& /* blockThreadExtent */,
        TVec const& /* threadElemExtent */,
        T const* /* ci */,
        T* /* co */,
        uint32_t /* size */,
        T* /* blockSums */,
        int32_t numBlocks,
        int32_t const* /* pc */,
        std::size_t warpSize) {
      if constexpr (cms::alpakatools::requires_single_thread_per_block_v<TAcc>) {
        return sizeof(T) * numBlocks;
      } else {
        return sizeof(T) * (warpSize + numBlocks);
      }
    }
  };

  // Decoupled lookback: blockPrefixScan workspace (T[warpSize]) + thread aggregates (T[threadsPerBlock]).
  template <alpaka::concepts::Acc TAcc, typename T, uint32_t ITEMS_PER_THREAD>
  struct BlockSharedMemDynSizeBytes<cms::alpakatools::scanDecoupledLookback<T, ITEMS_PER_THREAD>, TAcc> {
    template <typename TVec>
    ALPAKA_FN_HOST_ACC static std::size_t getBlockSharedMemDynSizeBytes(
        cms::alpakatools::scanDecoupledLookback<T, ITEMS_PER_THREAD> const& /* kernel */,
        TVec const& blockThreadExtent,
        TVec const& /* threadElemExtent */,
        T const* /* ci */,
        T* /* co */,
        uint32_t /* size */,
        uint64_t* /* tileDescriptors */,
        int32_t* /* tileCounter */,
        std::size_t warpSize) {
      if constexpr (cms::alpakatools::requires_single_thread_per_block_v<TAcc>) {
        return 0;
      } else {
        return sizeof(T) * (warpSize + static_cast<std::size_t>(blockThreadExtent[0u]));
      }
    }
  };

  // Kernel B does not require shared memory.
  template <alpaka::concepts::Acc TAcc, typename T>
  struct BlockSharedMemDynSizeBytes<cms::alpakatools::addScannedBlockOffsets<T>, TAcc> {
    template <typename TVec>
    ALPAKA_FN_HOST_ACC static std::size_t getBlockSharedMemDynSizeBytes(
        cms::alpakatools::addScannedBlockOffsets<T> const&, TVec const&, TVec const&, T const*, uint32_t, T const*) {
      return 0;
    }
  };

}  // namespace alpaka::trait

#endif  // HeterogeneousCore_AlpakaInterface_interface_prefixScan_h
