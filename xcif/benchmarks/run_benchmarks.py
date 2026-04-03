from __future__ import absolute_import, division, print_function
# xcif/benchmarks/run_benchmarks.py
#
# Wall-clock and peak-RSS benchmark: xcif vs gemmi vs ucif.
#
# Usage:
#   libtbx.python xcif/benchmarks/run_benchmarks.py <file.cif> [file2.cif ...] \
#       [--repeat N] [--warmup N] [--level1]
#
# Each parser runs in a subprocess so ru_maxrss is isolated per parser.

import sys
import os
import time
import math
import subprocess
import json
import resource
import platform

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _median(v):
  s = sorted(v)
  n = len(s)
  if n % 2 == 1:
    return s[n // 2]
  return (s[n // 2 - 1] + s[n // 2]) / 2.0

def _stddev(v, mean):
  return math.sqrt(sum((x - mean) ** 2 for x in v) / len(v))

def _rss_kb():
  """Peak RSS of current process in KB (platform-normalized)."""
  ru = resource.getrusage(resource.RUSAGE_SELF)
  rss = ru.ru_maxrss
  if platform.system() == "Darwin":
    rss = rss // 1024  # macOS reports bytes
  return rss

def _format_ms(ms):
  if ms < 1000:
    return "%.1f ms" % ms
  return "%.2f s" % (ms / 1000.0)

# ---------------------------------------------------------------------------
# Subprocess worker — runs a single parser in isolation
# ---------------------------------------------------------------------------

_WORKER_SCRIPT = r'''
from __future__ import absolute_import, division, print_function
import sys, os, time, json, resource, platform

def rss_kb():
  ru = resource.getrusage(resource.RUSAGE_SELF)
  rss = ru.ru_maxrss
  if platform.system() == "Darwin":
    rss = rss // 1024
  return rss

parser_name = sys.argv[1]
cif_path    = sys.argv[2]
repeat      = int(sys.argv[3])
warmup      = int(sys.argv[4])

parse_fn = None
if parser_name == "ucif":
  from iotbx.cif import reader as cif_reader
  parse_fn = lambda: cif_reader(file_path=cif_path)
elif parser_name == "gemmi":
  import gemmi
  parse_fn = lambda: gemmi.cif.read(cif_path)
elif parser_name == "xcif":
  import xcif
  parse_fn = lambda: xcif.reader(file_path=cif_path)
else:
  sys.exit("unknown parser: " + parser_name)

# Warmup
for _ in range(warmup):
  parse_fn()

# Timed runs
times = []
for _ in range(repeat):
  t0 = time.perf_counter()
  parse_fn()
  t1 = time.perf_counter()
  times.append((t1 - t0) * 1000.0)  # ms

peak_rss = rss_kb()

result = {"times_ms": times, "peak_rss_kb": peak_rss}
sys.stdout.write(json.dumps(result) + "\n")
'''

def _run_parser_subprocess(parser_name, cif_path, repeat, warmup, python_exe):
  """Run a parser in a child process; return (times_ms, peak_rss_kb) or None."""
  cmd = [python_exe, "-c", _WORKER_SCRIPT,
         parser_name, cif_path, str(repeat), str(warmup)]
  try:
    proc = subprocess.Popen(
      cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate(timeout=600)
  except Exception as e:
    return None, str(e)

  if proc.returncode != 0:
    err = stderr.decode("utf-8", errors="replace").strip()
    return None, err

  try:
    data = json.loads(stdout.decode("utf-8").strip())
    return data, None
  except Exception as e:
    return None, "JSON parse error: %s" % e

# ---------------------------------------------------------------------------
# Level 1 (pure C++) benchmark
# ---------------------------------------------------------------------------

def _run_level1(cif_path, repeat, warmup):
  """Run bench_xcif C++ binary; return parsed result dict or None."""
  try:
    import libtbx.load_env
    build_dir = libtbx.env.under_build("xcif")
    bench_bin = os.path.join(build_dir, "benchmarks", "bench_xcif")
    if not os.path.isfile(bench_bin):
      return None, "bench_xcif binary not found at %s" % bench_bin
  except Exception:
    return None, "libtbx not available; cannot locate bench_xcif"

  cmd = [bench_bin, cif_path, str(repeat), str(warmup)]
  try:
    proc = subprocess.Popen(
      cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate(timeout=600)
  except Exception as e:
    return None, str(e)

  if proc.returncode != 0:
    err = stderr.decode("utf-8", errors="replace").strip()
    return None, err

  # Parse the RESULT line from stdout
  for line in stdout.decode("utf-8").splitlines():
    if line.startswith("RESULT:"):
      parts = {}
      for token in line[len("RESULT:"):].split():
        k, v = token.split("=")
        parts[k] = float(v)
      return parts, None

  return None, "no RESULT line in bench_xcif output"

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
  args = sys.argv[1:]

  # Parse arguments
  cif_paths = []
  repeat = 5
  warmup = 3
  run_level1 = False
  i = 0
  while i < len(args):
    if args[i] == "--repeat":
      repeat = int(args[i + 1])
      i += 2
    elif args[i] == "--warmup":
      warmup = int(args[i + 1])
      i += 2
    elif args[i] == "--level1":
      run_level1 = True
      i += 1
    else:
      cif_paths.append(args[i])
      i += 1

  if not cif_paths:
    print("Usage: libtbx.python run_benchmarks.py <file.cif> [file2.cif ...]")
    print("       [--repeat N] [--warmup N] [--level1]")
    sys.exit(1)

  # Validate files
  for p in cif_paths:
    if not os.path.isfile(p):
      print("File not found: %s" % p)
      sys.exit(1)

  # Determine Python executable
  python_exe = sys.executable

  parsers = ["ucif", "gemmi", "xcif"]

  for cif_path in cif_paths:
    size_mb = os.path.getsize(cif_path) / (1024.0 * 1024.0)
    print("=" * 70)
    print("File: %s (%.1f MB)" % (os.path.basename(cif_path), size_mb))
    print("Warmup: %d  Repeat: %d" % (warmup, repeat))
    print("=" * 70)

    # ── Level 2 (Python boundary) ──────────────────────────────────
    print("\n--- Level 2: Python parse (subprocess-isolated) ---\n")
    results = {}

    for name in parsers:
      sys.stdout.write("  %-12s " % name)
      sys.stdout.flush()
      data, err = _run_parser_subprocess(
        name, cif_path, repeat, warmup, python_exe)
      if data is None:
        print("SKIP (%s)" % err.split("\n")[-1][:60])
        continue

      times = data["times_ms"]
      rss = data["peak_rss_kb"]
      med = _median(times)
      mn = min(times)
      mx = max(times)
      avg = sum(times) / len(times)
      sd = _stddev(times, avg)

      results[name] = {
        "median_ms": med, "min_ms": mn, "max_ms": mx,
        "stddev_ms": sd, "peak_rss_kb": rss,
      }
      print("median=%s  min=%s  max=%s  sd=%s  RSS=%.0f MB" % (
        _format_ms(med), _format_ms(mn), _format_ms(mx),
        _format_ms(sd), rss / 1024.0))

    # Summary table
    if len(results) >= 2:
      gemmi_med = results.get("gemmi", {}).get("median_ms")
      gemmi_rss = results.get("gemmi", {}).get("peak_rss_kb")
      print("\n  %-12s %10s %10s %10s" % ("Parser", "Median", "vs gemmi", "RSS MB"))
      print("  " + "-" * 46)
      for name in sorted(results, key=lambda k: results[k]["median_ms"]):
        r = results[name]
        ratio_str = ""
        if gemmi_med:
          ratio = r["median_ms"] / gemmi_med
          ratio_str = "%.2fx" % ratio
        print("  %-12s %10s %10s %10.0f" % (
          name, _format_ms(r["median_ms"]), ratio_str,
          r["peak_rss_kb"] / 1024.0))

    # ── Level 1 (pure C++) ─────────────────────────────────────────
    if run_level1:
      print("\n--- Level 1: Pure C++ parse (xcif only) ---\n")
      l1, err = _run_level1(cif_path, repeat, warmup)
      if l1 is None:
        print("  SKIP (%s)" % err)
      else:
        print("  xcif C++:  median=%s  min=%s  max=%s  sd=%s" % (
          _format_ms(l1["median_ms"]),
          _format_ms(l1["min_ms"]),
          _format_ms(l1["max_ms"]),
          _format_ms(l1["stddev_ms"])))
        # Compare Level 1 vs Level 2
        xcif_l2 = results.get("xcif", {}).get("median_ms")
        if xcif_l2:
          overhead = xcif_l2 - l1["median_ms"]
          print("  Python boundary overhead: %s (%.0f%%)" % (
            _format_ms(overhead),
            (overhead / l1["median_ms"]) * 100 if l1["median_ms"] > 0 else 0))

    print()

if __name__ == "__main__":
  main()
