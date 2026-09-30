// Author: 'bformless'
// Location: Munich, Germany
// (c) 2019
// 
// ------------------------------------------------------------------------------------------------------------------------------------------------------------ 
// License: AGPLv3-only
// ------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// Date: 07-31-2019
// changed on: 09-28-2026

#ifndef FAST_SIEVE_SIEVE_HPP
#define FAST_SIEVE_SIEVE_HPP

#include <cstdint>
#include <iosfwd>
#include <limits>

namespace fast_sieve {

// Unsigned 64-bit values are used for the complete supported search domain.
using u64 = std::uint64_t;

// The largest number accepted by the public range API and command-line interface.
inline constexpr u64 kMaximumSupportedNumber = (std::numeric_limits<u64>::max)();

// Dave's Garage benchmark always uses the canonical one-million upper bound.
inline constexpr u64 kDaveBenchmarkLimit = 1000000U;

// Summarizes a prime search over an inclusive numeric range.
struct SieveResult final {
    u64 count = 0U;
    u64 last_prime = 0U;
};

// Stores the canonical Dave's Garage benchmark result and validation state.
struct DaveBenchmarkResult final {
    u64 passes = 0U;
    double elapsed_seconds = 0.0;
    unsigned threads = 1U;
    SieveResult sieve{};
    bool valid = false;
};

// Searches all primes in the inclusive range [first_number, last_number].
[[nodiscard]] SieveResult run_range(u64 first_number,
                                    u64 last_number,
                                    unsigned thread_count = 1U);

// Searches the inclusive range and writes selected primes one per line in ascending order.
// The output interval must be a non-empty subrange of the searched interval.
[[nodiscard]] SieveResult run_range_and_write(u64 first_number,
                                              u64 last_number,
                                              u64 output_first_number,
                                              u64 output_last_number,
                                              unsigned thread_count,
                                              std::ostream& output);

// Returns the largest thread count that can perform independent range-segment work.
// The result never exceeds std::thread::hardware_concurrency() when that value is known.
[[nodiscard]] unsigned max_useful_threads(u64 first_number, u64 last_number) noexcept;

// Compatibility wrapper that searches [2, limit].
[[nodiscard]] SieveResult run(u64 limit, unsigned thread_count = 1U);

// Compatibility wrapper that searches [2, limit] and writes every prime one per line.
[[nodiscard]] SieveResult run_and_print(u64 limit, std::ostream& output);

// Returns the maximum benchmark worker count reported by the C++ runtime.
[[nodiscard]] unsigned max_benchmark_threads() noexcept;

// Executes and validates the canonical Dave's Garage benchmark workload.
[[nodiscard]] DaveBenchmarkResult run_daves_benchmark(unsigned thread_count = 1U,
                                                       double minimum_seconds = 5.0);

} // namespace fast_sieve

#endif
