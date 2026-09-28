// Author: 'bformless'
// Location: Munich, Germany
// (c) 2019
// 
// ------------------------------------------------------------------------------------------------------------------------------------------------------------ 
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files(the "Software"),
// to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and /or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Date: 07-31-2019
// changed on: 09-28-2026

#include "sieve.hpp"

#include <algorithm>
#include <array>
#include <barrier>
#include <bit>
#include <chrono>
#include <condition_variable>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <mutex>
#include <ostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace fast_sieve {
namespace {

// Compact unsigned types match the byte-oriented wheel representation and 32-bit base primes.
using u8 = std::uint8_t;
using u32 = std::uint32_t;

// A wheel modulo 30 stores only residues coprime to 2, 3, and 5.
constexpr std::array<u8, 8> kResidues{1U, 7U, 11U, 13U, 17U, 19U, 23U, 29U};
// Gaps advance from one represented residue to the next within a wheel cycle.
constexpr std::array<u8, 8> kWheelGaps{6U, 4U, 2U, 4U, 2U, 4U, 6U, 2U};
// These five primes form a reusable periodic presieve and are handled separately as primes.
constexpr std::array<u32, 5> kPresievePrimes{7U, 11U, 13U, 17U, 19U};
// Primes removed by the wheel or presieve are inserted explicitly into range results.
constexpr std::array<u64, 8> kSmallPrimes{2U, 3U, 5U, 7U, 11U, 13U, 17U, 19U};
// The presieve repeats after the product of its prime moduli, measured in 30-number blocks.
constexpr u32 kPresievePeriodBlocks = 7U * 11U * 13U * 17U * 19U;
// Each segment holds one byte per wheel block; this size balances locality and task granularity.
constexpr std::size_t kSegmentBlocks = 262144U;
// Known canonical values validate the one-million Dave's Garage benchmark workload.
constexpr u64 kDaveBenchmarkPrimeCount = 78498U;
constexpr u64 kDaveBenchmarkLastPrime = 999983U;
static_assert(kSegmentBlocks <= static_cast<std::size_t>((std::numeric_limits<u32>::max)()));
static_assert(kPresievePeriodBlocks == 323323U);
// UINT64_MAX has an integer square root below 2^32, so u32 can represent every required base prime.
static_assert(kMaximumSupportedNumber / (static_cast<u64>((std::numeric_limits<u32>::max)()) + 1U) <
              static_cast<u64>((std::numeric_limits<u32>::max)()) + 1U);

// Maps each residue modulo 30 to its bit position, or 0xFF when the residue is not represented.
[[nodiscard]] constexpr std::array<u8, 30> make_residue_index() noexcept {
    std::array<u8, 30> table{};
    table.fill(0xFFU);
    for (u8 i = 0U; i < static_cast<u8>(kResidues.size()); ++i) {
        table[kResidues[i]] = i;
    }
    return table;
}

constexpr auto kResidueIndex = make_residue_index();

// Precomputes the destination bit mask for every product of two represented wheel residues.
[[nodiscard]] constexpr std::array<std::array<u8, 8>, 8> make_product_masks() noexcept {
    std::array<std::array<u8, 8>, 8> table{};
    for (u8 pi = 0U; pi < 8U; ++pi) {
        for (u8 qi = 0U; qi < 8U; ++qi) {
            const u8 residue = static_cast<u8>((static_cast<unsigned>(kResidues[pi]) *
                                                static_cast<unsigned>(kResidues[qi])) % 30U);
            const u8 ri = kResidueIndex[residue];
            table[pi][qi] = static_cast<u8>(u8{1U} << ri);
        }
    }
    return table;
}

constexpr auto kProductMasks = make_product_masks();

// Precomputes block carries caused when a wheel-state step crosses a multiple-of-30 boundary.
[[nodiscard]] constexpr std::array<std::array<u8, 8>, 8> make_block_carries() noexcept {
    std::array<std::array<u8, 8>, 8> table{};
    for (u8 pi = 0U; pi < 8U; ++pi) {
        const unsigned p_residue = kResidues[pi];
        for (u8 qi = 0U; qi < 8U; ++qi) {
            const unsigned product_residue =
                (p_residue * static_cast<unsigned>(kResidues[qi])) % 30U;
            const unsigned delta = p_residue * static_cast<unsigned>(kWheelGaps[qi]);
            table[pi][qi] = static_cast<u8>((product_residue + delta) / 30U);
        }
    }
    return table;
}

constexpr auto kBlockCarries = make_block_carries();

// Gives the number of prime-sized steps needed to reach a residue represented by the wheel.
[[nodiscard]] constexpr std::array<std::array<u8, 30>, 8> make_adjust_steps() noexcept {
    std::array<std::array<u8, 30>, 8> table{};
    for (u8 pi = 0U; pi < 8U; ++pi) {
        const unsigned p_residue = kResidues[pi];
        for (unsigned residue = 0U; residue < 30U; ++residue) {
            u8 steps = 0U;
            unsigned current = residue;
            while (kResidueIndex[current] == 0xFFU) {
                current = (current + p_residue) % 30U;
                ++steps;
            }
            table[pi][residue] = steps;
        }
    }
    return table;
}

constexpr auto kAdjustSteps = make_adjust_steps();

// Recovers the q wheel state from p's residue index and the residue of p*q.
[[nodiscard]] constexpr std::array<std::array<u8, 30>, 8> make_q_state_from_product() noexcept {
    std::array<std::array<u8, 30>, 8> table{};
    for (auto& row : table) {
        row.fill(0xFFU);
    }
    for (u8 pi = 0U; pi < 8U; ++pi) {
        for (u8 qi = 0U; qi < 8U; ++qi) {
            const u8 product_residue = static_cast<u8>(
                (static_cast<unsigned>(kResidues[pi]) * static_cast<unsigned>(kResidues[qi])) % 30U);
            table[pi][product_residue] = qi;
        }
    }
    return table;
}

constexpr auto kQStateFromProduct = make_q_state_from_product();

// Eight wheel cycles span 8*p blocks, enabling safe wide OR operations for small primes.
constexpr std::array<u32, 22> kPackedStrikePrimes{
    23U, 29U, 31U, 37U, 41U, 43U, 47U, 53U, 59U, 61U, 67U,
    71U, 73U, 79U, 83U, 89U, 97U, 101U, 103U, 107U, 109U, 113U};
constexpr std::size_t kPackedStrikeMaxPrime = 113U;
constexpr std::size_t kPackedStrikeTileBytes = 8U * kPackedStrikeMaxPrime;

// Stores all eight q-state packed strike tiles for one small base prime.
struct PackedStrikePattern final {
    std::array<std::array<u8, kPackedStrikeTileBytes>, 8> bytes{};
};

// Builds compile-time packed strike tiles used to mark several wheel cycles with dense OR operations.
[[nodiscard]] constexpr std::array<PackedStrikePattern, kPackedStrikePrimes.size()>
make_packed_strike_patterns() noexcept {
    std::array<PackedStrikePattern, kPackedStrikePrimes.size()> patterns{};
    for (std::size_t pi = 0U; pi < kPackedStrikePrimes.size(); ++pi) {
        const u32 prime = kPackedStrikePrimes[pi];
        const u8 p_index = kResidueIndex[static_cast<u8>(prime % 30U)];
        const u32 p_blocks = prime / 30U;
        for (u8 q_start = 0U; q_start < 8U; ++q_start) {
            std::array<u32, 8> offsets{};
            std::array<u8, 8> masks{};
            u32 cycle_offset = 0U;
            for (u8 i = 0U; i < 8U; ++i) {
                const u8 state = static_cast<u8>((q_start + i) & 7U);
                offsets[i] = cycle_offset;
                masks[i] = kProductMasks[p_index][state];
                cycle_offset += p_blocks * static_cast<u32>(kWheelGaps[state]) +
                                static_cast<u32>(kBlockCarries[p_index][state]);
            }
            for (u32 cycle = 0U; cycle < 8U; ++cycle) {
                const u32 cycle_base = cycle * prime;
                for (u8 i = 0U; i < 8U; ++i) {
                    const std::size_t position = static_cast<std::size_t>(cycle_base + offsets[i]);
                    patterns[pi].bytes[q_start][position] = static_cast<u8>(
                        patterns[pi].bytes[q_start][position] | masks[i]);
                }
            }
        }
    }
    return patterns;
}

constexpr auto kPackedStrikePatterns = make_packed_strike_patterns();

// Maps an eligible small prime directly to its packed-pattern table index.
[[nodiscard]] constexpr std::array<u8, kPackedStrikeMaxPrime + 1U>
make_packed_strike_index() noexcept {
    std::array<u8, kPackedStrikeMaxPrime + 1U> index{};
    index.fill(0xFFU);
    for (u8 i = 0U; i < static_cast<u8>(kPackedStrikePrimes.size()); ++i) {
        index[kPackedStrikePrimes[i]] = i;
    }
    return index;
}

constexpr auto kPackedStrikeIndex = make_packed_strike_index();

// Verifies compile-time wheel invariants before any generated table is accepted.
[[nodiscard]] constexpr bool wheel_tables_are_consistent() noexcept {
    unsigned gap_sum = 0U;
    for (const u8 gap : kWheelGaps) {
        gap_sum += static_cast<unsigned>(gap);
    }
    if (gap_sum != 30U) {
        return false;
    }

    for (u8 pi = 0U; pi < 8U; ++pi) {
        unsigned carry_sum = 0U;
        for (u8 qi = 0U; qi < 8U; ++qi) {
            carry_sum += static_cast<unsigned>(kBlockCarries[pi][qi]);
            if (kProductMasks[pi][qi] == 0U) {
                return false;
            }
        }
        if (carry_sum != static_cast<unsigned>(kResidues[pi])) {
            return false;
        }

        for (unsigned residue = 0U; residue < 30U; ++residue) {
            const unsigned adjusted =
                (residue + static_cast<unsigned>(kAdjustSteps[pi][residue]) *
                               static_cast<unsigned>(kResidues[pi])) % 30U;
            if (kResidueIndex[adjusted] == 0xFFU ||
                kQStateFromProduct[pi][adjusted] == 0xFFU) {
                return false;
            }
        }
    }
    return true;
}

static_assert(wheel_tables_are_consistent());

// Computes floor(sqrt(value)) and corrects any floating-point seed error using exact division tests.
[[nodiscard]] u64 isqrt_u64(const u64 value) noexcept {
    if (value < 2U) {
        return value;
    }

    u64 candidate = static_cast<u64>(std::sqrt(static_cast<long double>(value)));
    while (candidate > value / candidate) {
        --candidate;
    }
    while (candidate < (std::numeric_limits<u64>::max)()) {
        const u64 next = candidate + 1U;
        if (next > value / next) {
            break;
        }
        candidate = next;
    }
    return candidate;
}

// Builds all base primes from 23 through floor(sqrt(limit)) with an odd-only bit sieve.
// For the uint64_t maximum, floor(sqrt(limit)) is UINT32_MAX, so every base prime fits in u32.
[[nodiscard]] std::vector<u32> build_base_primes(const u64 limit) {
    const u64 root = isqrt_u64(limit);
    if (root < 23U) {
        return {};
    }

    const u64 odd_count = ((root - 3U) >> 1U) + 1U;
    const std::size_t word_count = static_cast<std::size_t>((odd_count + 63U) >> 6U);
    std::vector<u64> composite(word_count, 0U);

    const u64 sieve_root = isqrt_u64(root);
    for (u64 index = 0U;; ++index) {
        const u64 prime = 3U + (index << 1U);
        if (prime > sieve_root) {
            break;
        }
        const std::size_t word = static_cast<std::size_t>(index >> 6U);
        const unsigned bit = static_cast<unsigned>(index & 63U);
        if (((composite[word] >> bit) & 1U) != 0U) {
            continue;
        }

        const u64 start = ((prime * prime) - 3U) >> 1U;
        for (u64 mark = start; mark < odd_count; mark += prime) {
            const std::size_t mark_word = static_cast<std::size_t>(mark >> 6U);
            const unsigned mark_bit = static_cast<unsigned>(mark & 63U);
            composite[mark_word] |= (u64{1U} << mark_bit);
        }
    }

    std::vector<u32> primes;
    if (root >= 100U) {
        const long double x = static_cast<long double>(root);
        const long double estimate = x / (std::log(x) - 1.0L);
        if (estimate > 0.0L && estimate < static_cast<long double>(primes.max_size())) {
            primes.reserve(static_cast<std::size_t>(estimate));
        }
    }

    const u64 first_index = (23U - 3U) >> 1U;
    for (u64 index = first_index; index < odd_count; ++index) {
        const std::size_t word = static_cast<std::size_t>(index >> 6U);
        const unsigned bit = static_cast<unsigned>(index & 63U);
        if (((composite[word] >> bit) & 1U) == 0U) {
            primes.push_back(static_cast<u32>(3U + (index << 1U)));
        }
    }
    return primes;
}

// Builds one full periodic composite pattern for the presieve primes 7, 11, 13, 17, and 19.
[[nodiscard]] std::vector<u8> make_presieve_pattern() {
    std::vector<u8> period(1U, u8{0U});
    u32 current_period = 1U;
    for (const u32 prime : kPresievePrimes) {
        const u32 next_period = current_period * prime;
        std::vector<u8> expanded(static_cast<std::size_t>(next_period));
        for (u32 offset = 0U; offset < next_period; offset += current_period) {
            std::memcpy(expanded.data() + static_cast<std::size_t>(offset),
                        period.data(), period.size());
        }

        const u8 p_index = kResidueIndex[static_cast<u8>(prime % 30U)];
        u64 multiple = static_cast<u64>(prime);
        u8 q_state = 0U;
        const u64 end_value = static_cast<u64>(next_period) * 30U;
        while (multiple < end_value) {
            const std::size_t block = static_cast<std::size_t>(multiple / 30U);
            expanded[block] = static_cast<u8>(expanded[block] |
                                               kProductMasks[p_index][q_state]);
            multiple += static_cast<u64>(prime) * static_cast<u64>(kWheelGaps[q_state]);
            q_state = static_cast<u8>((q_state + 1U) & 7U);
        }
        period = std::move(expanded);
        current_period = next_period;
    }

    return period;
}

// Lazily constructs the immutable presieve once and reuses it for all later segments.
[[nodiscard]] const std::vector<u8>& presieve_pattern() {
    static const std::vector<u8> value = make_presieve_pattern();
    return value;
}


// Counts the explicitly handled small primes that lie inside an inclusive range.
[[nodiscard]] constexpr SieveResult small_prime_result_range(const u64 first_number,
                                                             const u64 last_number) noexcept {
    SieveResult result{};
    for (const u64 prime : kSmallPrimes) {
        if (prime < first_number) {
            continue;
        }
        if (prime > last_number) {
            break;
        }
        ++result.count;
        result.last_prime = prime;
    }
    return result;
}

// Returns the legacy [2, limit] small-prime contribution used by the benchmark path.
[[nodiscard]] constexpr SieveResult small_prime_result(const u64 limit) noexcept {
    return limit < 2U ? SieveResult{} : small_prime_result_range(2U, limit);
}

// Describes the wheel-30 block interval that must actually be sieved for a search range.
struct RangeGeometry final {
    u64 wheel_first_number = 0U;
    u64 first_block = 0U;
    u64 last_block = 0U;
    u64 block_count = 0U;
    u64 segment_count = 0U;
};

// Builds overflow-safe block geometry for numbers at or above the first wheel candidate, 23.
[[nodiscard]] constexpr RangeGeometry make_range_geometry(const u64 first_number,
                                                          const u64 last_number) noexcept {
    if (first_number > last_number || last_number < 23U) {
        return {};
    }

    RangeGeometry geometry{};
    geometry.wheel_first_number = first_number < 23U ? 23U : first_number;
    geometry.first_block = geometry.wheel_first_number / 30U;
    geometry.last_block = last_number / 30U;
    geometry.block_count = geometry.last_block - geometry.first_block + 1U;

    const u64 segment_blocks = static_cast<u64>(kSegmentBlocks);
    geometry.segment_count = geometry.block_count / segment_blocks;
    if ((geometry.block_count % segment_blocks) != 0U) {
        ++geometry.segment_count;
    }
    return geometry;
}

// Initializes one segment from the reusable presieve and masks numbers outside the requested range.
void seed_segment(std::vector<u8>& composite,
                  const std::size_t block_count,
                  const u64 first_block,
                  const u64 first_number,
                  const u64 last_number) noexcept {
    const auto& base = presieve_pattern();
    const std::size_t period = base.size();
    std::size_t phase = static_cast<std::size_t>(first_block % static_cast<u64>(period));
    std::size_t written = 0U;

    // Copy the periodic 7*11*13*17*19 presieve pattern into this independent segment.
    while (written < block_count) {
        const std::size_t available = period - phase;
        const std::size_t amount = std::min(available, block_count - written);
        std::memcpy(composite.data() + written, base.data() + phase, amount);
        written += amount;
        phase = 0U;
    }

    // Number 1 occupies residue 1 of block zero but is not prime.
    if (first_block == 0U) {
        composite[0] = static_cast<u8>(composite[0] | u8{1U});
    }

    // Mask wheel residues below the inclusive lower search bound in the first searched block.
    const u64 search_first_block = first_number / 30U;
    if (first_block == search_first_block) {
        const u64 block_base = first_block * 30U;
        const u64 minimum_residue = first_number - block_base;
        u8& first = composite[0];
        for (u8 ri = 0U; ri < static_cast<u8>(kResidues.size()); ++ri) {
            if (static_cast<u64>(kResidues[ri]) < minimum_residue) {
                first = static_cast<u8>(first | static_cast<u8>(u8{1U} << ri));
            }
        }
    }

    // Mask wheel residues above the inclusive upper search bound in the last searched block.
    const u64 search_last_block = last_number / 30U;
    const u64 segment_last_block = first_block + static_cast<u64>(block_count) - 1U;
    if (segment_last_block == search_last_block) {
        const u64 block_base = search_last_block * 30U;
        const u64 maximum_residue = last_number - block_base;
        u8& last = composite[block_count - 1U];
        for (u8 ri = 0U; ri < static_cast<u8>(kResidues.size()); ++ri) {
            if (static_cast<u64>(kResidues[ri]) > maximum_residue) {
                last = static_cast<u8>(last | static_cast<u8>(u8{1U} << ri));
            }
        }
    }
}

// Marks composites generated by one base prime inside a bounded wheel-30 segment.
void mark_prime(std::vector<u8>& composite,
                const std::size_t block_count,
                const u64 first_block,
                const u64 low,
                const u64 high,
                const u32 prime) noexcept {
    const u64 p = static_cast<u64>(prime);
    const u64 square = p * p;
    if (square > high) {
        return;
    }

    const u8 p_residue = static_cast<u8>(prime % 30U);
    const u8 p_index = kResidueIndex[p_residue];

    // Start at p*p, or advance to the first represented multiple in this segment.
    u64 multiple = square;
    u8 q_state = p_index;
    if (multiple < low) {
        const u64 remainder = low % p;
        const u64 delta = remainder == 0U ? 0U : (p - remainder);
        if (low > (std::numeric_limits<u64>::max)() - delta) {
            return;
        }
        multiple = low + delta;

        const u8 product_residue = static_cast<u8>(multiple % 30U);
        const u8 adjust = kAdjustSteps[p_index][product_residue];
        if (adjust != 0U) {
            const u64 add = p * static_cast<u64>(adjust);
            if (multiple > (std::numeric_limits<u64>::max)() - add) {
                return;
            }
            multiple += add;
        }
        if (multiple > high) {
            return;
        }
        q_state = kQStateFromProduct[p_index][static_cast<u8>(multiple % 30U)];
    }

    const u64 absolute_block = multiple / 30U;
    if (absolute_block < first_block) {
        return;
    }
    u32 block = static_cast<u32>(absolute_block - first_block);
    const u32 block_limit = static_cast<u32>(block_count);
    if (block >= block_limit) {
        return;
    }

    // Precompute the eight wheel offsets and masks for one complete q-state cycle.
    const u32 p_blocks = prime / 30U;
    std::array<u32, 8> offsets{};
    std::array<u8, 8> masks{};
    u32 cycle_offset = 0U;
    for (u8 offset = 0U; offset < 8U; ++offset) {
        const u8 state = static_cast<u8>((q_state + offset) & 7U);
        offsets[offset] = cycle_offset;
        masks[offset] = kProductMasks[p_index][state];
        cycle_offset += p_blocks * static_cast<u32>(kWheelGaps[state]) +
                        static_cast<u32>(kBlockCarries[p_index][state]);
    }

    // Small strike primes use a precomputed packed pattern to reduce per-composite work.
    if (prime <= static_cast<u32>(kPackedStrikeMaxPrime)) {
        const u8 packed_index = kPackedStrikeIndex[prime];
        if (packed_index != 0xFFU) {
            const std::size_t tile_bytes = static_cast<std::size_t>(prime) * 8U;
            const u8* const pattern = kPackedStrikePatterns[packed_index].bytes[q_state].data();
            const u32 tile_blocks = prime * 8U;
            while (tile_blocks <= block_limit - block) {
                u8* const destination = composite.data() + block;
                if (prime <= 61U) {
                    for (std::size_t byte = 0U; byte < tile_bytes; byte += sizeof(u64)) {
                        u64 data_word = 0U;
                        u64 mask_word = 0U;
                        std::memcpy(&data_word, destination + byte, sizeof(data_word));
                        std::memcpy(&mask_word, pattern + byte, sizeof(mask_word));
                        data_word |= mask_word;
                        std::memcpy(destination + byte, &data_word, sizeof(data_word));
                    }
                } else {
                    for (std::size_t byte = 0U; byte < tile_bytes; ++byte) {
                        destination[byte] = static_cast<u8>(destination[byte] | pattern[byte]);
                    }
                }
                block += tile_blocks;
            }
        }
    }

    // Apply complete eight-state wheel cycles while the whole cycle fits in the segment.
    while (prime <= block_limit - block) {
        u8* const base = composite.data() + block;
        base[offsets[0]] = static_cast<u8>(base[offsets[0]] | masks[0]);
        base[offsets[1]] = static_cast<u8>(base[offsets[1]] | masks[1]);
        base[offsets[2]] = static_cast<u8>(base[offsets[2]] | masks[2]);
        base[offsets[3]] = static_cast<u8>(base[offsets[3]] | masks[3]);
        base[offsets[4]] = static_cast<u8>(base[offsets[4]] | masks[4]);
        base[offsets[5]] = static_cast<u8>(base[offsets[5]] | masks[5]);
        base[offsets[6]] = static_cast<u8>(base[offsets[6]] | masks[6]);
        base[offsets[7]] = static_cast<u8>(base[offsets[7]] | masks[7]);
        block += prime;
    }

    // Finish the partial wheel cycle at the right edge without accessing beyond the segment.
    const u32 remaining = block_limit - block;
    for (u8 offset = 0U; offset < 8U && offsets[offset] < remaining; ++offset) {
        const u32 index = block + offsets[offset];
        composite[index] = static_cast<u8>(composite[index] | masks[offset]);
    }
}

// Counts clear candidate bits; every clear bit represents one prime after sieving.
[[nodiscard]] u64 count_clear_bits(const std::vector<u8>& composite,
                                   const std::size_t block_count) noexcept {
    u64 set_bits = 0U;
    std::size_t index = 0U;
    for (; index + sizeof(u64) <= block_count; index += sizeof(u64)) {
        u64 word = 0U;
        std::memcpy(&word, composite.data() + index, sizeof(word));
        set_bits += static_cast<u64>(std::popcount(word));
    }
    for (; index < block_count; ++index) {
        set_bits += static_cast<u64>(std::popcount(static_cast<unsigned>(composite[index])));
    }
    return static_cast<u64>(block_count) * 8U - set_bits;
}

// Finds the numerically greatest clear wheel candidate in a processed segment.
[[nodiscard]] u64 find_last_prime(const std::vector<u8>& composite,
                                  const std::size_t block_count,
                                  const u64 first_block) noexcept {
    constexpr u64 maximum = (std::numeric_limits<u64>::max)();
    for (std::size_t index = block_count; index > 0U; --index) {
        const u8 available = static_cast<u8>(~composite[index - 1U]);
        if (available == 0U) {
            continue;
        }
        const u64 absolute_block = first_block + static_cast<u64>(index - 1U);
        const u64 base = absolute_block * 30U;
        for (int ri = 7; ri >= 0; --ri) {
            const u8 mask = static_cast<u8>(u8{1U} << static_cast<unsigned>(ri));
            const u64 residue = static_cast<u64>(kResidues[static_cast<std::size_t>(ri)]);
            if ((available & mask) != 0U && residue <= maximum - base) {
                return base + residue;
            }
        }
    }
    return 0U;
}

// Processes an arbitrary contiguous block span of at most kSegmentBlocks wheel blocks.
[[nodiscard]] SieveResult process_block_span(std::vector<u8>& composite,
                                             const u64 first_block,
                                             const std::size_t block_count,
                                             const u64 first_number,
                                             const u64 last_number,
                                             const std::vector<u32>& base_primes) noexcept {
    seed_segment(composite, block_count, first_block, first_number, last_number);

    const u64 block_low = first_block * 30U;
    const u64 low = std::max(first_number, block_low);
    const u64 end_block = first_block + static_cast<u64>(block_count);
    const u64 max_block_plus_one = ((std::numeric_limits<u64>::max)() / 30U) + 1U;
    const u64 natural_high = end_block >= max_block_plus_one
                                 ? (std::numeric_limits<u64>::max)()
                                 : end_block * 30U - 1U;
    const u64 high = std::min(last_number, natural_high);

    for (const u32 prime : base_primes) {
        const u64 p = static_cast<u64>(prime);
        if (p * p > high) {
            break;
        }
        mark_prime(composite, block_count, first_block, low, high, prime);
    }

    SieveResult result{};
    result.count = count_clear_bits(composite, block_count);
    if (result.count != 0U) {
        result.last_prime = find_last_prime(composite, block_count, first_block);
    }
    return result;
}

// Implements the allocation behavior expected by the canonical Dave's Garage benchmark.
class DaveBenchmarkSieve final {
public:
    explicit DaveBenchmarkSieve(const u64 limit)
        : limit_(limit),
          composite_(static_cast<std::size_t>((limit / 30U) + 1U)) {}

    // Runs the single-buffer benchmark sieve and returns its validated summary values.
    [[nodiscard]] SieveResult run() noexcept {
        const std::size_t block_count = composite_.size();
        seed_segment(composite_, block_count, 0U, 0U, limit_);

        const u64 root = isqrt_u64(limit_);
        for (u64 block = 0U; block < static_cast<u64>(block_count); ++block) {
            const u64 base = block * 30U;
            for (u8 ri = 0U; ri < static_cast<u8>(kResidues.size()); ++ri) {
                const u64 prime = base + static_cast<u64>(kResidues[ri]);
                if (prime < 23U) {
                    continue;
                }
                if (prime > root) {
                    return finish();
                }
                const u8 mask = static_cast<u8>(u8{1U} << ri);
                if ((composite_[static_cast<std::size_t>(block)] & mask) == 0U) {
                    mark_prime(composite_, block_count, 0U, 0U, limit_,
                               static_cast<u32>(prime));
                }
            }
        }
        return finish();
    }

private:
    // Adds the five explicitly handled primes that are excluded by the presieve pattern.
    [[nodiscard]] SieveResult finish() const noexcept {
        SieveResult result = small_prime_result(limit_);
        result.count += count_clear_bits(composite_, composite_.size());
        result.last_prime = find_last_prime(composite_, composite_.size(), 0U);
        return result;
    }

    u64 limit_;
    std::vector<u8> composite_;
};

// Stores one independent benchmark worker's pass count and final sieve result.
struct DaveWorkerState final {
    u64 passes = 0U;
    SieveResult last{};
    std::exception_ptr error{};
};

// Runs benchmark passes until the common deadline after synchronizing worker start time.
void run_dave_worker(DaveWorkerState& state,
                     std::barrier<>& start_barrier,
                     const std::chrono::steady_clock::time_point& deadline) noexcept {
    try {
        start_barrier.arrive_and_wait();
        do {
            DaveBenchmarkSieve sieve(kDaveBenchmarkLimit);
            state.last = sieve.run();
            ++state.passes;
        } while (std::chrono::steady_clock::now() < deadline);
    } catch (...) {
        state.error = std::current_exception();
    }
}

// Owns reusable segment storage for one independent normal-search worker.
struct WorkerState final {
    explicit WorkerState(const std::size_t block_capacity) : composite(block_capacity) {}

    std::vector<u8> composite;
    SieveResult result{};
    u64 first_block = 0U;
    std::size_t block_count = 0U;
};

// Processes one range-relative segment and records enough geometry for ordered output.
void process_relative_segment(WorkerState& state,
                              const u64 relative_segment,
                              const RangeGeometry& geometry,
                              const u64 last_number,
                              const std::vector<u32>& base_primes) noexcept {
    const u64 segment_blocks = static_cast<u64>(kSegmentBlocks);
    state.first_block = geometry.first_block + relative_segment * segment_blocks;
    const u64 remaining_blocks = geometry.last_block - state.first_block + 1U;
    state.block_count = static_cast<std::size_t>(std::min(segment_blocks, remaining_blocks));
    state.result = process_block_span(state.composite,
                                      state.first_block,
                                      state.block_count,
                                      geometry.wheel_first_number,
                                      last_number,
                                      base_primes);
}

// Reuses a fixed set of worker threads for ordered output batches without per-segment thread creation.
class SegmentBatchPool final {
public:
    SegmentBatchPool(const unsigned worker_count,
                     std::vector<WorkerState>& states,
                     const RangeGeometry& geometry,
                     const u64 last_number,
                     const std::vector<u32>& base_primes)
        : states_(states),
          geometry_(geometry),
          last_number_(last_number),
          base_primes_(base_primes) {
        threads_.reserve(worker_count);
        try {
            for (unsigned worker = 0U; worker < worker_count; ++worker) {
                threads_.emplace_back([this, worker]() noexcept { worker_loop(worker); });
            }
        } catch (...) {
            stop();
            throw;
        }
    }

    SegmentBatchPool(const SegmentBatchPool&) = delete;
    SegmentBatchPool& operator=(const SegmentBatchPool&) = delete;

    ~SegmentBatchPool() {
        stop();
    }

    // Starts one generation of at most worker_count independent segments and waits for completion.
    void process_batch(const u64 batch_start, const unsigned active_workers) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            batch_start_ = batch_start;
            active_workers_ = active_workers;
            completed_workers_ = 0U;
            ++generation_;
        }
        start_cv_.notify_all();

        std::unique_lock<std::mutex> lock(mutex_);
        done_cv_.wait(lock, [this, active_workers] {
            return completed_workers_ == active_workers;
        });
    }

private:
    // Waits for new work generations and processes only the slot assigned to this worker index.
    void worker_loop(const unsigned worker) noexcept {
        u64 observed_generation = 0U;
        std::unique_lock<std::mutex> lock(mutex_);
        for (;;) {
            start_cv_.wait(lock, [this, observed_generation] {
                return stopping_ || generation_ != observed_generation;
            });
            if (stopping_) {
                return;
            }

            observed_generation = generation_;
            const u64 batch_start = batch_start_;
            const unsigned active_workers = active_workers_;
            if (worker >= active_workers) {
                continue;
            }

            lock.unlock();
            process_relative_segment(states_[worker],
                                     batch_start + static_cast<u64>(worker),
                                     geometry_,
                                     last_number_,
                                     base_primes_);
            lock.lock();

            ++completed_workers_;
            if (completed_workers_ == active_workers_) {
                done_cv_.notify_one();
            }
        }
    }

    // Wakes all sleeping workers and joins them exactly once, including constructor-failure cleanup.
    void stop() noexcept {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (stopping_) {
                return;
            }
            stopping_ = true;
        }
        start_cv_.notify_all();
        for (auto& thread : threads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }

    std::vector<WorkerState>& states_;
    const RangeGeometry& geometry_;
    u64 last_number_;
    const std::vector<u32>& base_primes_;
    std::vector<std::jthread> threads_;
    std::mutex mutex_;
    std::condition_variable start_cv_;
    std::condition_variable done_cv_;
    u64 generation_ = 0U;
    u64 batch_start_ = 0U;
    unsigned active_workers_ = 0U;
    unsigned completed_workers_ = 0U;
    bool stopping_ = false;
};

// Assigns disjoint segments to one counting worker with a fixed cyclic stride.
void worker_run(WorkerState& state,
                const u64 first_relative_segment,
                const u64 stride,
                const RangeGeometry& geometry,
                const u64 last_number,
                const std::vector<u32>& base_primes) noexcept {
    SieveResult local{};
    for (u64 segment = first_relative_segment;
         segment < geometry.segment_count;
         segment += stride) {
        process_relative_segment(state, segment, geometry, last_number, base_primes);
        local.count += state.result.count;
        local.last_prime = std::max(local.last_prime, state.result.last_prime);
    }
    state.result = local;
}

// Appends one decimal prime followed by a newline to an output buffer.
void append_prime_line(std::string& buffer, const u64 value) {
    std::array<char, 32> digits{};
    const auto converted = std::to_chars(digits.data(), digits.data() + digits.size(), value);
    if (converted.ec != std::errc{}) {
        throw std::runtime_error("internal decimal conversion failed");
    }
    buffer.append(digits.data(), converted.ptr);
    buffer.push_back('\n');
}

// Appends all primes from one processed segment that also lie in the selected output range.
void append_segment_primes(std::string& output,
                           const std::vector<u8>& composite,
                           const std::size_t block_count,
                           const u64 first_block,
                           const u64 output_first_number,
                           const u64 output_last_number) {
    constexpr u64 maximum = (std::numeric_limits<u64>::max)();
    for (std::size_t block = 0U; block < block_count; ++block) {
        const u8 available = static_cast<u8>(~composite[block]);
        if (available == 0U) {
            continue;
        }

        const u64 absolute_block = first_block + static_cast<u64>(block);
        const u64 base = absolute_block * 30U;
        for (u8 ri = 0U; ri < static_cast<u8>(kResidues.size()); ++ri) {
            const u8 mask = static_cast<u8>(u8{1U} << ri);
            const u64 residue = static_cast<u64>(kResidues[ri]);
            if ((available & mask) == 0U || residue > maximum - base) {
                continue;
            }
            const u64 prime = base + residue;
            if (prime >= output_first_number && prime <= output_last_number) {
                append_prime_line(output, prime);
            }
        }
    }
}

// Flushes a bounded text buffer to the selected stream and verifies stream state.
void flush_output_buffer(std::ostream& output, std::string& buffer) {
    if (buffer.empty()) {
        return;
    }
    output.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    if (!output) {
        throw std::runtime_error("output stream write failed");
    }
    buffer.clear();
}

} // namespace

// Returns a conservative useful upper bound: no more workers than independent segments or hardware threads.
unsigned max_useful_threads(const u64 first_number, const u64 last_number) noexcept {
    const RangeGeometry geometry = make_range_geometry(first_number, last_number);
    if (geometry.segment_count <= 1U) {
        return 1U;
    }

    unsigned hardware = std::thread::hardware_concurrency();
    if (hardware == 0U) {
        hardware = 1U;
    }
    const u64 bounded = std::min<u64>(static_cast<u64>(hardware), geometry.segment_count);
    return static_cast<unsigned>(bounded);
}

// Searches an inclusive range without producing the prime list.
SieveResult run_range(const u64 first_number,
                      const u64 last_number,
                      unsigned thread_count) {
    if (first_number < 2U || first_number > last_number) {
        throw std::invalid_argument("invalid prime search range");
    }

    SieveResult total = small_prime_result_range(first_number, last_number);
    const RangeGeometry geometry = make_range_geometry(first_number, last_number);
    if (geometry.segment_count == 0U) {
        return total;
    }

    // Materialize the presieve before entering noexcept worker code so allocation failures propagate.
    static_cast<void>(presieve_pattern());

    if (thread_count == 0U) {
        thread_count = 1U;
    }
    thread_count = std::min(thread_count, max_useful_threads(first_number, last_number));

    // Base primes only need to reach sqrt(last_number), independent of first_number.
    const std::vector<u32> base_primes = build_base_primes(last_number);
    const unsigned workers = static_cast<unsigned>(
        std::min<u64>(static_cast<u64>(thread_count), geometry.segment_count));
    const std::size_t worker_block_capacity = static_cast<std::size_t>(
        std::min<u64>(static_cast<u64>(kSegmentBlocks), geometry.block_count));

    std::vector<WorkerState> states;
    states.reserve(workers);
    for (unsigned worker = 0U; worker < workers; ++worker) {
        states.emplace_back(worker_block_capacity);
    }

    if (workers == 1U) {
        worker_run(states[0], 0U, 1U, geometry, last_number, base_primes);
    } else {
        std::vector<std::jthread> threads;
        threads.reserve(workers);
        for (unsigned worker = 0U; worker < workers; ++worker) {
            threads.emplace_back([&, worker]() noexcept {
                worker_run(states[worker],
                           static_cast<u64>(worker),
                           static_cast<u64>(workers),
                           geometry,
                           last_number,
                           base_primes);
            });
        }
        for (auto& thread : threads) {
            thread.join();
        }
    }

    for (const WorkerState& state : states) {
        total.count += state.result.count;
        total.last_prime = std::max(total.last_prime, state.result.last_prime);
    }
    return total;
}

// Searches the range and emits only the selected output subrange in strictly ascending order.
SieveResult run_range_and_write(const u64 first_number,
                                const u64 last_number,
                                const u64 output_first_number,
                                const u64 output_last_number,
                                unsigned thread_count,
                                std::ostream& output) {
    if (first_number < 2U || first_number > last_number ||
        output_first_number < first_number || output_last_number > last_number ||
        output_first_number > output_last_number) {
        throw std::invalid_argument("invalid search or output range");
    }

    SieveResult total = small_prime_result_range(first_number, last_number);
    std::string buffer;
    buffer.reserve(1U << 20U);

    // The five presieved small primes are written explicitly when selected for output.
    for (const u64 prime : kSmallPrimes) {
        if (prime >= output_first_number && prime <= output_last_number) {
            append_prime_line(buffer, prime);
        }
    }

    const RangeGeometry geometry = make_range_geometry(first_number, last_number);
    if (geometry.segment_count != 0U) {
        // Materialize the presieve before entering noexcept worker code so allocation failures propagate.
        static_cast<void>(presieve_pattern());

        if (thread_count == 0U) {
            thread_count = 1U;
        }
        thread_count = std::min(thread_count, max_useful_threads(first_number, last_number));

        const std::vector<u32> base_primes = build_base_primes(last_number);
        const unsigned workers = static_cast<unsigned>(
            std::min<u64>(static_cast<u64>(thread_count), geometry.segment_count));
        const std::size_t worker_block_capacity = static_cast<std::size_t>(
            std::min<u64>(static_cast<u64>(kSegmentBlocks), geometry.block_count));

        std::vector<WorkerState> states;
        states.reserve(workers);
        for (unsigned worker = 0U; worker < workers; ++worker) {
            states.emplace_back(worker_block_capacity);
        }

        // Process at most one segment per worker at a time, then serialize completed segments in order.
        // A persistent pool bounds memory and avoids repeatedly creating threads for very large ranges.
        if (workers == 1U) {
            for (u64 segment = 0U; segment < geometry.segment_count; ++segment) {
                process_relative_segment(states[0], segment, geometry, last_number, base_primes);
                total.count += states[0].result.count;
                total.last_prime = std::max(total.last_prime, states[0].result.last_prime);
                append_segment_primes(buffer,
                                      states[0].composite,
                                      states[0].block_count,
                                      states[0].first_block,
                                      output_first_number,
                                      output_last_number);
                if (buffer.size() >= (1U << 20U)) {
                    flush_output_buffer(output, buffer);
                }
            }
        } else {
            SegmentBatchPool pool(workers, states, geometry, last_number, base_primes);
            for (u64 batch_start = 0U; batch_start < geometry.segment_count;
                 batch_start += static_cast<u64>(workers)) {
                const unsigned active = static_cast<unsigned>(std::min<u64>(
                    static_cast<u64>(workers), geometry.segment_count - batch_start));
                pool.process_batch(batch_start, active);

                for (unsigned worker = 0U; worker < active; ++worker) {
                    total.count += states[worker].result.count;
                    total.last_prime = std::max(total.last_prime, states[worker].result.last_prime);
                    append_segment_primes(buffer,
                                          states[worker].composite,
                                          states[worker].block_count,
                                          states[worker].first_block,
                                          output_first_number,
                                          output_last_number);
                    if (buffer.size() >= (1U << 20U)) {
                        flush_output_buffer(output, buffer);
                    }
                }
            }
        }
    }

    flush_output_buffer(output, buffer);
    return total;
}

// Preserves the previous [2, limit] counting API.
SieveResult run(const u64 limit, const unsigned thread_count) {
    if (limit < 2U) {
        return {};
    }
    return run_range(2U, limit, thread_count);
}

// Preserves the previous [2, limit] output API while using the new one-prime-per-line format.
SieveResult run_and_print(const u64 limit, std::ostream& output) {
    if (limit < 2U) {
        return {};
    }
    return run_range_and_write(2U, limit, 2U, limit, 1U, output);
}

// Returns the runtime-reported hardware-thread ceiling for the canonical benchmark.
unsigned max_benchmark_threads() noexcept {
    const unsigned hardware = std::thread::hardware_concurrency();
    return hardware == 0U ? 1U : hardware;
}

// Executes the canonical 5-second benchmark and validates the known one-million result.
DaveBenchmarkResult run_daves_benchmark(unsigned thread_count, double minimum_seconds) {
    if (!(minimum_seconds > 0.0)) {
        minimum_seconds = 5.0;
    }
    if (thread_count == 0U) {
        thread_count = 1U;
    }
    thread_count = std::min(thread_count, max_benchmark_threads());

    // Materialize the reusable wheel seed before timing benchmark passes.
    static_cast<void>(presieve_pattern());

    DaveBenchmarkResult result{};
    result.threads = thread_count;

    const auto duration = std::chrono::duration<double>(minimum_seconds);
    if (thread_count == 1U) {
        const auto start = std::chrono::steady_clock::now();
        const auto deadline = start +
            std::chrono::duration_cast<std::chrono::steady_clock::duration>(duration);
        do {
            DaveBenchmarkSieve sieve(kDaveBenchmarkLimit);
            result.sieve = sieve.run();
            ++result.passes;
        } while (std::chrono::steady_clock::now() < deadline);
        const auto end = std::chrono::steady_clock::now();
        result.elapsed_seconds = std::chrono::duration<double>(end - start).count();
    } else {
        std::vector<DaveWorkerState> states(static_cast<std::size_t>(thread_count));
        std::barrier start_barrier(static_cast<std::ptrdiff_t>(thread_count) + 1);
        auto deadline = std::chrono::steady_clock::time_point{};
        std::vector<std::jthread> workers;
        workers.reserve(static_cast<std::size_t>(thread_count));
        for (unsigned worker = 0U; worker < thread_count; ++worker) {
            workers.emplace_back([&, worker]() noexcept {
                run_dave_worker(states[static_cast<std::size_t>(worker)],
                                start_barrier,
                                deadline);
            });
        }

        const auto start = std::chrono::steady_clock::now();
        deadline = start +
            std::chrono::duration_cast<std::chrono::steady_clock::duration>(duration);
        start_barrier.arrive_and_wait();
        for (auto& worker : workers) {
            worker.join();
        }
        const auto end = std::chrono::steady_clock::now();
        result.elapsed_seconds = std::chrono::duration<double>(end - start).count();

        bool first = true;
        for (const DaveWorkerState& state : states) {
            if (state.error) {
                std::rethrow_exception(state.error);
            }
            result.passes += state.passes;
            if (first) {
                result.sieve = state.last;
                first = false;
            } else if (state.last.count != result.sieve.count ||
                       state.last.last_prime != result.sieve.last_prime) {
                result.valid = false;
                return result;
            }
        }
    }

    result.valid = result.passes != 0U &&
                   result.sieve.count == kDaveBenchmarkPrimeCount &&
                   result.sieve.last_prime == kDaveBenchmarkLastPrime;
    return result;
}

} // namespace fast_sieve
