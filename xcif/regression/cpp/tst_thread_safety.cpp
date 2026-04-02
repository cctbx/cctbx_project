// cctbx_project/xcif/regression/cpp/tst_thread_safety.cpp
//
// Step 14: Thread-safety tests.
// Verifies that xcif can safely parse distinct documents on multiple
// threads simultaneously and that concurrent read access to a shared
// parsed Document is safe.

#include <xcif/data_model.h>
#include <xcif/numeric.h>
#include "test_utils.h"

#include <atomic>
#include <string>
#include <thread>
#include <vector>

using namespace xcif;

// ─── Helper: generate a unique CIF string for each thread ──────────

static std::string make_cif(int id, int n_rows) {
  std::string s = "data_block_" + std::to_string(id) + "\n";
  s += "_entry.id " + std::to_string(id) + "\n";
  s += "_cell.length_a " + std::to_string(10.0 + id) + "\n";
  s += "loop_\n_atom.id\n_atom.x\n_atom.y\n_atom.z\n";
  for (int r = 0; r < n_rows; ++r) {
    s += std::to_string(r + 1) + " "
       + std::to_string(1.0 * r + 0.1 * id) + " "
       + std::to_string(2.0 * r) + " "
       + std::to_string(3.0 * r) + "\n";
  }
  return s;
}

// ─── Test 1: Concurrent independent parses ─────────────────────────

static std::atomic<int> concurrent_parse_failures(0);

void test_concurrent_parse() {
  const int N_THREADS = 8;
  const int N_ROWS = 200;
  std::vector<std::thread> threads;

  for (int t = 0; t < N_THREADS; ++t) {
    threads.emplace_back([t, N_ROWS]() {
      std::string cif = make_cif(t, N_ROWS);
      Document doc = parse(cif.c_str(), cif.size(), "<thread>");

      // Verify block name
      if (doc.size() != 1) {
        ++concurrent_parse_failures;
        return;
      }
      std::string expected_name = "block_" + std::to_string(t);
      if (doc[0].name() != expected_name) {
        ++concurrent_parse_failures;
        return;
      }
      // Verify entry id
      string_view eid = doc[0].find_value("_entry.id");
      if (eid != std::to_string(t)) {
        ++concurrent_parse_failures;
        return;
      }
      // Verify loop dimensions
      const Loop* lp = doc[0].find_loop("_atom");
      if (lp == NULL || lp->width() != 4 || lp->length() != (std::size_t)N_ROWS) {
        ++concurrent_parse_failures;
        return;
      }
      // Spot-check a value
      std::vector<string_view> ids = lp->column("_atom.id");
      if (ids.size() != (std::size_t)N_ROWS || ids[0] != "1") {
        ++concurrent_parse_failures;
        return;
      }
    });
  }
  for (auto& th : threads) th.join();
  CHECK_EQ(concurrent_parse_failures.load(), 0);
}

// ─── Test 2: Concurrent reads on a shared Document ─────────────────

static std::atomic<int> concurrent_read_failures(0);

void test_concurrent_reads() {
  // Parse once, then read from multiple threads
  std::string cif = make_cif(42, 500);
  Document doc = parse(cif.c_str(), cif.size(), "<shared>");

  const int N_THREADS = 8;
  std::vector<std::thread> threads;

  for (int t = 0; t < N_THREADS; ++t) {
    threads.emplace_back([&doc, t]() {
      // Each thread reads the same document
      if (doc.size() != 1) { ++concurrent_read_failures; return; }
      const Block& blk = doc[0];
      if (blk.name() != "block_42") { ++concurrent_read_failures; return; }

      string_view cell_a = blk.find_value("_cell.length_a");
      if (cell_a.empty()) { ++concurrent_read_failures; return; }

      const Loop* lp = blk.find_loop("_atom");
      if (lp == NULL) { ++concurrent_read_failures; return; }
      if (lp->length() != 500) { ++concurrent_read_failures; return; }

      // Read all columns
      std::vector<string_view> ids = lp->column("_atom.id");
      std::vector<string_view> xs  = lp->column("_atom.x");
      std::vector<string_view> ys  = lp->column("_atom.y");
      std::vector<string_view> zs  = lp->column("_atom.z");
      if (ids.size() != 500 || xs.size() != 500) {
        ++concurrent_read_failures;
        return;
      }

      // Numeric conversions on shared data
      for (std::size_t i = 0; i < ids.size(); ++i) {
        as_int(ids[i]);
        as_double(xs[i]);
      }
    });
  }
  for (auto& th : threads) th.join();
  CHECK_EQ(concurrent_read_failures.load(), 0);
}

// ─── Test 3: Concurrent parse + numeric on same thread ─────────────

static std::atomic<int> mixed_failures(0);

void test_concurrent_parse_and_numeric() {
  const int N_THREADS = 8;
  std::vector<std::thread> threads;

  for (int t = 0; t < N_THREADS; ++t) {
    threads.emplace_back([t]() {
      // Parse a small CIF
      std::string cif = make_cif(t, 50);
      Document doc = parse(cif.c_str(), cif.size(), "<mixed>");
      const Loop* lp = doc[0].find_loop("_atom");
      if (lp == NULL) { ++mixed_failures; return; }

      // Extract and convert numerics
      std::vector<string_view> xs = lp->column("_atom.x");
      double sum = 0.0;
      for (std::size_t i = 0; i < xs.size(); ++i) {
        sum += as_double(xs[i]);
      }
      // Just verify we got a number (not NaN, not crash)
      if (sum != sum) { // NaN check
        ++mixed_failures;
      }
    });
  }
  for (auto& th : threads) th.join();
  CHECK_EQ(mixed_failures.load(), 0);
}

// ─── Test 4: Parse errors on multiple threads don't corrupt state ──

static std::atomic<int> error_thread_failures(0);

void test_concurrent_parse_errors() {
  const int N_THREADS = 8;
  std::vector<std::thread> threads;

  for (int t = 0; t < N_THREADS; ++t) {
    threads.emplace_back([t]() {
      // Alternate between valid and invalid CIF
      bool expect_error = (t % 2 == 0);
      std::string cif;
      if (expect_error) {
        cif = "_tag value\n"; // missing data_ header
      } else {
        cif = make_cif(t, 10);
      }
      try {
        Document doc = parse(cif.c_str(), cif.size(), "<error_test>");
        if (expect_error) {
          ++error_thread_failures; // should have thrown
        } else {
          if (doc.size() != 1) ++error_thread_failures;
        }
      } catch (const CifError&) {
        if (!expect_error) {
          ++error_thread_failures; // should not have thrown
        }
      } catch (...) {
        ++error_thread_failures; // unexpected exception type
      }
    });
  }
  for (auto& th : threads) th.join();
  CHECK_EQ(error_thread_failures.load(), 0);
}

// ─── Test 5: Many small parses (stress test) ───────────────────────

static std::atomic<int> stress_failures(0);

void test_stress_many_parses() {
  const int N_THREADS = 4;
  const int N_PARSES = 100;
  std::vector<std::thread> threads;

  for (int t = 0; t < N_THREADS; ++t) {
    threads.emplace_back([t, N_PARSES]() {
      for (int i = 0; i < N_PARSES; ++i) {
        int id = t * N_PARSES + i;
        std::string cif = make_cif(id, 5);
        Document doc = parse(cif.c_str(), cif.size(), "<stress>");
        if (doc.size() != 1) { ++stress_failures; return; }
        std::string expected = "block_" + std::to_string(id);
        if (doc[0].name() != expected) { ++stress_failures; return; }
      }
    });
  }
  for (auto& th : threads) th.join();
  CHECK_EQ(stress_failures.load(), 0);
}

// ─── Main ──────────────────────────────────────────────────────────

void run_all_tests() {
  test_concurrent_parse();
  test_concurrent_reads();
  test_concurrent_parse_and_numeric();
  test_concurrent_parse_errors();
  test_stress_many_parses();
}

XCIF_TEST_MAIN()
