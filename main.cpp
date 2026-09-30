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

#include "sieve.hpp"

#include <charconv>
#include <chrono>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <string>
#include <string_view>

namespace {

using fast_sieve::u64;

constexpr std::string_view kDefaultDaveLabel = "bformless_fast_sieve_cpp20_wheel30_p19_pack113";

// Holds the fully parsed command-line configuration for a normal range search.
struct ProgramOptions final {
    u64 first_number = 0U;
    u64 last_number = 0U;
    u64 output_first_number = 0U;
    u64 output_last_number = 0U;
    unsigned threads = 1U;
    std::string output_file;
    bool have_first_number = false;
    bool have_last_number = false;
    bool have_output_first_number = false;
    bool have_output_last_number = false;
    bool use_output_file = false;
    bool output_console_explicit = false;
    bool have_threads = false;
};

// Converts a complete decimal token to uint64_t without locale-dependent parsing.
[[nodiscard]] bool parse_u64(const std::string_view text, u64& value) noexcept {
    if (text.empty()) {
        return false;
    }
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value, 10);
    return parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size();
}

// Converts a complete decimal token to unsigned for validated thread counts.
[[nodiscard]] bool parse_unsigned(const std::string_view text, unsigned& value) noexcept {
    if (text.empty()) {
        return false;
    }
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value, 10);
    return parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size();
}

// Rejects benchmark labels that would corrupt the semicolon-separated benchmark record.
[[nodiscard]] bool valid_benchmark_label(const std::string_view label) noexcept {
    if (label.empty()) {
        return false;
    }
    for (const char ch : label) {
        if (ch == ';' || ch == '\r' || ch == '\n') {
            return false;
        }
    }
    return true;
}

// Prints the complete non-interactive command-line interface.
void print_help(std::ostream& output, const char* executable) {
    output
        << "SIEVE OF ERATOSTHENES - segmented wheel-30 - (p)(c)2019 by bformless\n\n"
        << "Search mode:\n"
        << "  " << executable << " --firstnumber <N> --lastnumber <N> [options]\n\n"
        << "Required search range:\n"
        << "  --firstnumber <N>        Inclusive first number. Minimum: 2.\n"
        << "  --lastnumber <N>         Inclusive last number. Maximum: "
        << fast_sieve::kMaximumSupportedNumber << ".\n"
        << "                           This is the numeric uint64_t limit; practical\n"
        << "                           runtime and memory depend on the requested range.\n\n"
        << "Prime output:\n"
        << "  By default, no individual prime numbers are written. Only summary metadata\n"
        << "  is reported on stderr after the search completes.\n"
        << "  --outputconsole          Write prime numbers to stdout, one prime per line.\n"
        << "  --outputfile <PATH>      Write prime numbers to PATH, one prime per line.\n"
        << "  --output_startnumber <N> First prime value eligible for output.\n"
        << "                           Default when output is enabled: --firstnumber.\n"
        << "  --output_endnumber <N>   Last prime value eligible for output.\n"
        << "                           Default when output is enabled: --lastnumber.\n"
        << "                           Output bounds require --outputconsole or --outputfile\n"
        << "                           and must stay inside the search range.\n\n"
        << "Parallelism:\n"
        << "  --threads <N>            Requested sieve worker count. Default: 1.\n"
        << "                           The effective count is automatically capped at the\n"
        << "                           smaller of available hardware threads and independent\n"
        << "                           range segments. Requests above that cap are accepted\n"
        << "                           and reduced instead of being rejected. This avoids\n"
        << "                           useless oversubscription but cannot guarantee that\n"
        << "                           every additional thread is faster on every system.\n\n"
        << "Benchmark mode:\n"
        << "  " << executable << " --DavesBenchmark [--threads <N>]"
        << " [--benchmark-label <TEXT>]\n"
        << "  --DavesBenchmark         Run the canonical 5-second, 1,000,000 limit benchmark.\n"
        << "  --benchmark-label <TEXT> Override the benchmark record label. The label must\n"
        << "                           be non-empty and contain no semicolon, CR, or LF.\n\n"
        << "General:\n"
        << "  --help, -h               Print this complete help text.\n\n"
        << "Examples:\n"
        << "  " << executable
        << " --firstnumber 2 --lastnumber 1000000\n"
        << "                           Search only; do not list individual primes.\n"
        << "  " << executable
        << " --firstnumber 2 --lastnumber 1000000000 --threads 8"
        << " --outputfile primes.txt\n"
        << "  " << executable
        << " --firstnumber 2 --lastnumber 1000000 --outputconsole"
        << " --output_startnumber 1000 --output_endnumber 2000\n";
}

// Reports a command-line error together with the complete help text.
[[nodiscard]] int command_line_error(const char* executable, const std::string_view message) {
    std::cerr << "error: " << message << "\n\n";
    print_help(std::cerr, executable);
    return 2;
}

// Returns the following argv token or reports a missing option value.
[[nodiscard]] bool take_value(const int argc,
                              char** argv,
                              int& index,
                              std::string_view& value) noexcept {
    if (index + 1 >= argc) {
        return false;
    }
    ++index;
    value = argv[index];
    return true;
}

// Parses and validates all normal search-mode switches without interactive input.
[[nodiscard]] int parse_search_options(const int argc,
                                       char** argv,
                                       ProgramOptions& options) {
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument = argv[index];
        std::string_view value;

        if (argument == "--firstnumber") {
            if (options.have_first_number || !take_value(argc, argv, index, value) ||
                !parse_u64(value, options.first_number)) {
                return command_line_error(argv[0], "invalid or duplicate --firstnumber");
            }
            options.have_first_number = true;
        } else if (argument == "--lastnumber") {
            if (options.have_last_number || !take_value(argc, argv, index, value) ||
                !parse_u64(value, options.last_number)) {
                return command_line_error(argv[0], "invalid or duplicate --lastnumber");
            }
            options.have_last_number = true;
        } else if (argument == "--threads") {
            if (options.have_threads || !take_value(argc, argv, index, value) ||
                !parse_unsigned(value, options.threads) || options.threads == 0U) {
                return command_line_error(argv[0],
                                          "--threads requires one integer greater than zero and may appear only once");
            }
            options.have_threads = true;
        } else if (argument == "--outputfile") {
            if (options.use_output_file || options.output_console_explicit ||
                !take_value(argc, argv, index, value) || value.empty()) {
                return command_line_error(argv[0],
                                          "--outputfile requires one non-empty path and is mutually exclusive with --outputconsole");
            }
            options.output_file.assign(value);
            options.use_output_file = true;
        } else if (argument == "--outputconsole") {
            if (options.output_console_explicit || options.use_output_file) {
                return command_line_error(argv[0],
                                          "--outputconsole is duplicate or conflicts with --outputfile");
            }
            options.output_console_explicit = true;
        } else if (argument == "--output_startnumber") {
            if (options.have_output_first_number || !take_value(argc, argv, index, value) ||
                !parse_u64(value, options.output_first_number)) {
                return command_line_error(argv[0],
                                          "invalid or duplicate --output_startnumber");
            }
            options.have_output_first_number = true;
        } else if (argument == "--output_endnumber") {
            if (options.have_output_last_number || !take_value(argc, argv, index, value) ||
                !parse_u64(value, options.output_last_number)) {
                return command_line_error(argv[0],
                                          "invalid or duplicate --output_endnumber");
            }
            options.have_output_last_number = true;
        } else if (argument == "--help" || argument == "-h") {
            print_help(std::cout, argv[0]);
            return 1;
        } else {
            return command_line_error(argv[0], std::string("unknown option: ") + std::string(argument));
        }
    }

    if (!options.have_first_number || !options.have_last_number) {
        return command_line_error(argv[0],
                                  "both --firstnumber and --lastnumber are required in search mode");
    }
    if (options.first_number < 2U) {
        return command_line_error(argv[0], "--firstnumber must be at least 2");
    }
    if (options.first_number > options.last_number) {
        return command_line_error(argv[0], "--firstnumber must not exceed --lastnumber");
    }

    const bool output_enabled = options.use_output_file || options.output_console_explicit;
    if (!output_enabled &&
        (options.have_output_first_number || options.have_output_last_number)) {
        return command_line_error(
            argv[0],
            "--output_startnumber and --output_endnumber require --outputconsole or --outputfile");
    }

    if (output_enabled) {
        if (!options.have_output_first_number) {
            options.output_first_number = options.first_number;
        }
        if (!options.have_output_last_number) {
            options.output_last_number = options.last_number;
        }
        if (options.output_first_number < options.first_number ||
            options.output_last_number > options.last_number ||
            options.output_first_number > options.output_last_number) {
            return command_line_error(
                argv[0],
                "output bounds must form a non-empty subrange of [--firstnumber, --lastnumber]");
        }
    }

    // Accept an oversized request and reduce it to the amount of independent work.
    const unsigned maximum_threads =
        fast_sieve::max_useful_threads(options.first_number, options.last_number);
    if (options.threads > maximum_threads) {
        std::cerr << "note: requested --threads " << options.threads
                  << " reduced to " << maximum_threads
                  << " for this range and runtime\n";
        options.threads = maximum_threads;
    }

    return 0;
}

// Runs the preserved Dave's Garage benchmark through the same switch-based thread syntax.
[[nodiscard]] int run_daves_benchmark_mode(const int argc, char** argv) {
    unsigned threads = 1U;
    std::string_view label = kDefaultDaveLabel;
    bool have_threads = false;
    bool have_label = false;

    for (int index = 2; index < argc; ++index) {
        const std::string_view argument = argv[index];
        std::string_view value;
        if (argument == "--threads") {
            if (have_threads || !take_value(argc, argv, index, value) ||
                !parse_unsigned(value, threads) || threads == 0U) {
                return command_line_error(argv[0], "invalid benchmark --threads value");
            }
            have_threads = true;
        } else if (argument == "--benchmark-label") {
            if (have_label || !take_value(argc, argv, index, value) ||
                !valid_benchmark_label(value)) {
                return command_line_error(argv[0], "invalid benchmark label");
            }
            label = value;
            have_label = true;
        } else if (argument == "--help" || argument == "-h") {
            print_help(std::cout, argv[0]);
            return 0;
        } else {
            return command_line_error(argv[0], "invalid benchmark option");
        }
    }

    // Apply the same non-failing cap policy to benchmark worker requests.
    const unsigned maximum_threads = fast_sieve::max_benchmark_threads();
    if (threads > maximum_threads) {
        std::cerr << "note: requested benchmark --threads " << threads
                  << " reduced to " << maximum_threads
                  << " for this runtime\n";
        threads = maximum_threads;
    }

    const fast_sieve::DaveBenchmarkResult result =
        fast_sieve::run_daves_benchmark(threads, 5.0);
    const double average = result.passes == 0U
        ? 0.0
        : result.elapsed_seconds / static_cast<double>(result.passes);

    std::cerr.imbue(std::locale::classic());
    std::cout.imbue(std::locale::classic());
    std::cerr << std::setprecision(9)
              << "Passes: " << result.passes
              << ", Time: " << result.elapsed_seconds
              << ", Avg: " << average
              << ", Limit: " << fast_sieve::kDaveBenchmarkLimit
              << ", Count: " << result.sieve.count
              << ", Valid: " << (result.valid ? "true" : "false") << '\n';

    if (!result.valid) {
        return 1;
    }

    std::cout << label << ';'
              << result.passes << ';'
              << std::fixed << std::setprecision(6) << result.elapsed_seconds << ';'
              << result.threads
              << ";algorithm=wheel,faithful=yes,bits=1\n";
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    // Disable synchronized C stdio because all program I/O uses C++ streams.
    std::ios::sync_with_stdio(false);

    try {
        // With no arguments the program performs no work and prints the complete help text.
        if (argc == 1) {
            print_help(std::cout, argv[0]);
            return 0;
        }

        // Keep the benchmark as a separate deterministic command mode.
        if (std::string_view(argv[1]) == "--DavesBenchmark") {
            return run_daves_benchmark_mode(argc, argv);
        }

        ProgramOptions options{};
        const int parse_status = parse_search_options(argc, argv, options);
        if (parse_status != 0) {
            return parse_status == 1 ? 0 : parse_status;
        }

        // Prime listing is opt-in. Without an output switch, run the faster count-only path.
        const bool output_enabled = options.use_output_file || options.output_console_explicit;
        std::ofstream file_output;
        std::ostream* prime_output = nullptr;
        if (options.use_output_file) {
            file_output.open(options.output_file, std::ios::out | std::ios::trunc);
            if (!file_output.is_open()) {
                std::cerr << "fatal: cannot open output file: " << options.output_file << '\n';
                return 1;
            }
            prime_output = &file_output;
        } else if (options.output_console_explicit) {
            prime_output = &std::cout;
        }

        const auto start = std::chrono::steady_clock::now();
        fast_sieve::SieveResult result{};
        if (output_enabled) {
            result = fast_sieve::run_range_and_write(
                options.first_number,
                options.last_number,
                options.output_first_number,
                options.output_last_number,
                options.threads,
                *prime_output);
            prime_output->flush();
            if (!*prime_output) {
                std::cerr << "fatal: output stream write failed\n";
                return 1;
            }
        } else {
            result = fast_sieve::run_range(
                options.first_number, options.last_number, options.threads);
        }
        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end - start;

        // Keep metadata on stderr so stdout and output files contain prime values only.
        std::cerr.imbue(std::locale::classic());
        std::cerr << std::setprecision(9)
                  << "elapsed_seconds=" << elapsed.count() << '\n'
                  << "search_firstnumber=" << options.first_number << '\n'
                  << "search_lastnumber=" << options.last_number << '\n'
                  << "output_mode="
                  << (options.use_output_file ? "file"
                      : options.output_console_explicit ? "console" : "none") << '\n';
        if (output_enabled) {
            std::cerr << "output_startnumber=" << options.output_first_number << '\n'
                      << "output_endnumber=" << options.output_last_number << '\n';
        }
        std::cerr << "threads=" << options.threads << '\n'
                  << "prime_count=" << result.count << '\n'
                  << "last_prime=" << result.last_prime << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "fatal: " << error.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "fatal: unknown exception\n";
        return 1;
    }
}
