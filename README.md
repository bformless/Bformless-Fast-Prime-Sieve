# Sieve of Eratosthenes - Segmented Wheel-30 C++20

A portable C++20 prime sieve for inclusive 64-bit search ranges. The implementation combines a segmented wheel modulo 30, a reusable presieve, packed strike patterns for small base primes, optional multithreading, ordered prime output, and a built-in Dave's Garage benchmark mode.

The source uses the C++ standard library only. No compiler intrinsics, inline assembly, or platform-specific source paths are required.

## Requirements

- C++20 compiler and standard library
- `std::barrier`, `std::jthread`, and `<bit>` support
- Exact-width `std::uint32_t` and `std::uint64_t` integer types
- Hosted C++ implementation with thread support for multithreaded operation

The numeric command-line domain extends through `UINT64_MAX` (`18446744073709551615`). This is a numeric interface limit, not a practical runtime or memory guarantee. Base-prime generation grows with `sqrt(--lastnumber)`, so very large upper bounds can require substantial time and memory.

## Build

### GCC

```sh
g++ -std=c++20 -O3 -DNDEBUG -pthread -Wall -Wextra -Wpedantic main.cpp sieve.cpp -o bformless_fast_sieve
```

### Clang

```sh
clang++ -std=c++20 -O3 -DNDEBUG -pthread -Wall -Wextra -Wpedantic main.cpp sieve.cpp -o bformless_fast_sieve
```

### MSVC

```bat
cl /std:c++20 /O2 /EHsc /W4 main.cpp sieve.cpp /Fe:bformless_fast_sieve.exe
```

The exact compiler, flags, runtime library, operating system, and hardware affect benchmark results. Architecture-specific flags can be added for local experiments, but they are not required by the source code.

## Search mode

A normal search requires both inclusive range bounds:

```sh
bformless_fast_sieve --firstnumber 2 --lastnumber 1000000
```

By default, individual primes are not printed. The program reports summary metadata on `stderr`, including elapsed time, effective thread count, prime count, and the last prime found.

To write primes to standard output:

```sh
bformless_fast_sieve --firstnumber 2 --lastnumber 1000000 --outputconsole
```

To write primes to a file:

```sh
bformless_fast_sieve --firstnumber 2 --lastnumber 1000000000 --threads 8 --outputfile primes.txt
```

Output can be restricted to a subrange while the full search range is still sieved and counted:

```sh
bformless_fast_sieve --firstnumber 2 --lastnumber 1000000 --outputconsole --output_startnumber 1000 --output_endnumber 2000
```

Prime output is one decimal prime per line and remains in ascending numeric order. `stdout` is reserved for prime output or the benchmark record; summary and diagnostic metadata are written to `stderr`.

## Dave's Garage benchmark mode

Run the built-in five-second benchmark with the canonical limit of `1,000,000`:

```sh
bformless_fast_sieve --DavesBenchmark --threads 1
```

A successful benchmark validates a prime count of `78498` and last prime `999983`. The final semicolon-separated record uses the default label `bformless_fast_sieve_cpp20_wheel30_p19_pack113` unless `--benchmark-label` is supplied.

Example with a custom label:

```sh
bformless_fast_sieve --DavesBenchmark --threads 14 --benchmark-label my_test
```

## Exit status

- `0`: successful search, successful benchmark, or help display
- `1`: runtime failure or failed benchmark validation
- `2`: invalid command line

## Command-line options

| Option | Meaning |
|---|---|
| `--firstnumber <N>` | Inclusive first search value. Required in search mode. Minimum: `2`. |
| `--lastnumber <N>` | Inclusive last search value. Required in search mode. Maximum: `18446744073709551615`. |
| `--outputconsole` | Write selected primes to `stdout`, one per line. Mutually exclusive with `--outputfile`. |
| `--outputfile <PATH>` | Write selected primes to `PATH`, one per line. Mutually exclusive with `--outputconsole`. |
| `--output_startnumber <N>` | First prime value eligible for output. Requires an output mode. Default: `--firstnumber`. |
| `--output_endnumber <N>` | Last prime value eligible for output. Requires an output mode. Default: `--lastnumber`. |
| `--threads <N>` | Requested worker count. Must be greater than zero. In search mode the effective count is capped by available hardware threads and independent range segments; benchmark mode is capped by available hardware threads. Default: `1`. |
| `--DavesBenchmark` | Run the canonical five-second benchmark with limit `1,000,000`. This selects benchmark mode instead of normal search mode. |
| `--benchmark-label <TEXT>` | Override the first field of the benchmark record. Benchmark mode only. The label must be non-empty and contain no semicolon, carriage return, or line feed. |
| `--help`, `-h` | Print the complete help text. |

Running the program without arguments also prints the help text and performs no sieve work.

## How the sieve works

The sieve uses a wheel modulo 30, so each 30-number block stores only the eight residues coprime to `2`, `3`, and `5`: `1, 7, 11, 13, 17, 19, 23, 29`. A reusable periodic presieve removes multiples of `7`, `11`, `13`, `17`, and `19`. Base primes from `23` through `floor(sqrt(last_number))` then strike composite candidates inside independent segments of up to `262144` wheel blocks. Small strike primes through `113` use precomputed packed patterns; the remaining base primes use the general wheel strike path. After sieving, every clear candidate bit represents a prime, while the small primes excluded by the wheel and presieve are handled explicitly.

Independent segments can be processed by multiple workers. Count-only searches aggregate segment results directly. When prime values are requested, processed segments are emitted in range order so the output remains strictly ascending.

## Benchmarks

Reference measurements supplied for the release:

- CPU: i9-9940X
- Clock: 3.3 GHz
- Hyper-Threading: on
- Turbo: off
- SpeedStep: off
- Compiler: MSVC
- Build flags: GS GL W4 Gy Zc:wchar_t Zi Gm- O2 sdl Zc:inline fp:precise D "NDEBUG" "_CONSOLE" "_UNICODE" GT WX- Zc:forScope Gd Oi MD std:c++20 FC EHsc nologo Ot
- Operating System: Windows 10

These are single recorded runs under the stated machine settings. No repeated-run dispersion data were supplied, so the values document those measurements only and are not a universal performance claim.

### 14 threads

```text
bformless_fast_sieve.exe --DavesBenchmark --threads 14
Passes: 724315, Time: 5.0004656, Avg: 6.90371675e-06, Limit: 1000000, Count: 78498, Valid: true
bformless_fast_sieve_cpp20_wheel30_p19_pack113;724315;5.000466;14;algorithm=wheel,faithful=yes,bits=1
```

### 1 thread

```text
bformless_fast_sieve.exe --DavesBenchmark --threads 1
Passes: 44571, Time: 5.0000346, Avg: 0.000112181342, Limit: 1000000, Count: 78498, Valid: true
bformless_fast_sieve_cpp20_wheel30_p19_pack113;44571;5.000035;1;algorithm=wheel,faithful=yes,bits=1
```

### 1 Trillion Number Search with 14 threads

```text
bformless_fast_sieve.exe --firstnumber 2 --lastnumber 1000000000000 --threads 14
elapsed_seconds=41.0583033
search_firstnumber=2
search_lastnumber=1000000000000
output_mode=none
threads=14
prime_count=37607912018
last_prime=999999999989
```

### 1 Trillion Number Search with 1 thread

```text
bformless_fast_sieve.exe --firstnumber 2 --lastnumber 1000000000000 --threads 1
elapsed_seconds=561.441015
search_firstnumber=2
search_lastnumber=1000000000000
output_mode=none
threads=1
prime_count=37607912018
last_prime=999999999989
```