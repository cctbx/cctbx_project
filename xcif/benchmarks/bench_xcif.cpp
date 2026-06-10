// cctbx_project/xcif/benchmarks/bench_xcif.cpp
//
// Level 1 (pure C++) benchmark: measures wall-clock time for
// xcif::parse() on a memory-mapped CIF file.
//
// Usage:  bench_xcif <file.cif> [repeat] [warmup]
//         repeat  — timed iterations (default 5)
//         warmup  — untimed warmup iterations (default 3)
//
// Output: per-iteration times and a machine-parseable RESULT line.

#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>

#include "xcif/mapped_file.h"
#include "xcif/data_model.h"

static double median(std::vector<double>& v) {
  std::size_t n = v.size();
  std::sort(v.begin(), v.end());
  if (n % 2 == 1) return v[n / 2];
  return (v[n / 2 - 1] + v[n / 2]) / 2.0;
}

static double stddev(const std::vector<double>& v, double mean) {
  double sum = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    double d = v[i] - mean;
    sum += d * d;
  }
  return std::sqrt(sum / v.size());
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr,
      "Usage: %s <file.cif> [repeat] [warmup]\n"
      "  repeat  — timed iterations (default 5)\n"
      "  warmup  — untimed warmup iterations (default 3)\n",
      argv[0]);
    return 1;
  }

  const char* path = argv[1];
  int repeat = (argc >= 3) ? std::atoi(argv[2]) : 5;
  int warmup = (argc >= 4) ? std::atoi(argv[3]) : 3;

  if (repeat < 1) repeat = 1;
  if (warmup < 0) warmup = 0;

  // Check file exists and report size.
  xcif::MappedFile mf_check(path);
  std::fprintf(stderr, "File: %s (%.1f MB)\n",
               path, mf_check.size() / (1024.0 * 1024.0));
  std::fprintf(stderr, "Warmup: %d  Repeat: %d\n\n", warmup, repeat);

  // Warmup: parse and discard (uses zero-copy parse_file).
  for (int i = 0; i < warmup; ++i) {
    xcif::Document doc = xcif::parse_file(path);
    (void)doc;
  }

  // Timed runs.
  typedef std::chrono::high_resolution_clock Clock;
  std::vector<double> times_ms;
  times_ms.reserve(repeat);

  for (int i = 0; i < repeat; ++i) {
    Clock::time_point t0 = Clock::now();
    xcif::Document doc = xcif::parse_file(path);
    Clock::time_point t1 = Clock::now();
    (void)doc;

    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    times_ms.push_back(ms);
    std::fprintf(stderr, "  [%d] %.1f ms\n", i + 1, ms);
  }

  // Statistics.
  double sum = 0.0;
  for (std::size_t i = 0; i < times_ms.size(); ++i) sum += times_ms[i];
  double mean = sum / times_ms.size();
  double med = median(times_ms);
  double sd = stddev(times_ms, mean);
  double mn = *std::min_element(times_ms.begin(), times_ms.end());
  double mx = *std::max_element(times_ms.begin(), times_ms.end());

  std::fprintf(stderr,
    "\n  median=%.1f ms  min=%.1f ms  max=%.1f ms  stddev=%.1f ms\n",
    med, mn, mx, sd);

  // Machine-parseable result line on stdout.
  std::printf("RESULT: median_ms=%.2f min_ms=%.2f max_ms=%.2f stddev_ms=%.2f\n",
              med, mn, mx, sd);

  return 0;
}
