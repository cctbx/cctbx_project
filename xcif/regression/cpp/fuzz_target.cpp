// cctbx_project/xcif/regression/cpp/fuzz_target.cpp
//
// Step 14: Fuzz target for the xcif parser.
//
// When built with a fuzzing engine (e.g., libFuzzer via -fsanitize=fuzzer),
// the engine calls LLVMFuzzerTestOneInput with random byte buffers.
//
// When built without a fuzzing engine (normal build), the standalone main()
// runs a built-in corpus of adversarial inputs to verify the parser never
// crashes on malformed data.  This mode is what runs in the test suite.

#include <xcif/data_model.h>
#include <xcif/numeric.h>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <stdexcept>

// ── Core fuzz function ─────────────────────────────────────────────
// Contract: MUST NOT crash or trigger UB.  Exceptions are expected
// and swallowed.  Returns 0 always (libFuzzer convention).

static int fuzz_one(const uint8_t* data, size_t size) {
  // 1. Fuzz the parser
  try {
    xcif::Document doc = xcif::parse(
      reinterpret_cast<const char*>(data), size, "<fuzz>");
    // If parse succeeds, exercise the AST to catch use-after-free etc.
    for (std::size_t i = 0; i < doc.size(); ++i) {
      const xcif::Block& blk = doc[i];
      (void)blk.name();
      const std::vector<std::pair<xcif::string_view, xcif::string_view>>&
        pairs = blk.pairs();
      for (std::size_t p = 0; p < pairs.size(); ++p) {
        (void)pairs[p].first.size();
        (void)pairs[p].second.size();
      }
      for (std::size_t li = 0; li < blk.loops().size(); ++li) {
        const xcif::Loop& lp = blk.loops()[li];
        for (std::size_t r = 0; r < lp.length(); ++r) {
          for (std::size_t c = 0; c < lp.width(); ++c) {
            (void)lp.value(r, c);
          }
        }
      }
    }
  } catch (const xcif::CifError&) {
    // expected for malformed input
  } catch (const std::exception&) {
    // other C++ exceptions are acceptable
  }

  // 2. Fuzz numeric conversions on small slices
  if (size > 0 && size <= 64) {
    std::string s(reinterpret_cast<const char*>(data), size);
    xcif::string_view sv(s.data(), s.size());
    try { (void)xcif::as_double(sv); } catch (...) {}
    try { (void)xcif::as_int(sv); } catch (...) {}
    try { (void)xcif::as_double_with_su(sv); } catch (...) {}
    (void)xcif::is_null(sv);
    (void)xcif::is_unknown(sv);
    (void)xcif::is_inapplicable(sv);
  }

  return 0;
}

// ── libFuzzer entry point ──────────────────────────────────────────
#ifdef XCIF_FUZZ_ENGINE

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  return fuzz_one(data, size);
}

#else
// ── Standalone adversarial corpus ──────────────────────────────────
// Runs as a normal test in the suite.  Each input exercises an edge
// case that has historically caused parser crashes in CIF libraries.

static int failures = 0;

static void run_case(const char* name, const char* input, size_t len) {
  int rc = fuzz_one(reinterpret_cast<const uint8_t*>(input), len);
  if (rc != 0) {
    std::fprintf(stderr, "FAIL: fuzz case '%s' returned %d\n", name, rc);
    ++failures;
  }
}

static void run_case(const char* name, const char* input) {
  run_case(name, input, std::strlen(input));
}

int main() {
  // Empty and minimal inputs
  run_case("empty", "", 0);
  run_case("null_byte", "\0", 1);
  run_case("single_newline", "\n");
  run_case("just_whitespace", "   \t\n\n  \t\n");

  // Missing data_ header
  run_case("tag_without_data", "_tag value\n");
  run_case("loop_without_data", "loop_\n_tag\nval\n");
  run_case("save_without_data", "save_foo\nsave_\n");

  // Truncated keywords
  run_case("truncated_data", "data_");
  run_case("truncated_loop", "data_x\nloop_");
  run_case("truncated_save", "data_x\nsave_");
  run_case("truncated_tag", "data_x\n_");

  // Unterminated strings
  run_case("unterminated_single_quote", "data_x\n_t 'unterminated");
  run_case("unterminated_double_quote", "data_x\n_t \"unterminated");
  run_case("unterminated_semicolon", "data_x\n_t\n;unterminated text");

  // Loop edge cases
  run_case("loop_no_tags", "data_x\nloop_\n1 2 3\n");
  run_case("loop_no_values", "data_x\nloop_\n_a\n_b\n");
  run_case("loop_uneven_values", "data_x\nloop_\n_a\n_b\n1 2 3\n");
  run_case("loop_single_value", "data_x\nloop_\n_a\nv\n");

  // Deeply nested save frames
  run_case("save_frame", "data_d\nsave_s1\n_t v\nsave_\n");
  run_case("save_no_end", "data_d\nsave_s1\n_t v\n");

  // Very long tag name
  {
    std::string long_tag = "data_x\n_";
    long_tag.append(10000, 'a');
    long_tag += " value\n";
    run_case("long_tag", long_tag.c_str(), long_tag.size());
  }

  // Very long value
  {
    std::string long_val = "data_x\n_t ";
    long_val.append(10000, 'z');
    long_val += "\n";
    run_case("long_value", long_val.c_str(), long_val.size());
  }

  // Very long semicolon text field
  {
    std::string long_semi = "data_x\n_t\n;";
    long_semi.append(10000, 'q');
    long_semi += "\n;\n";
    run_case("long_semicolon", long_semi.c_str(), long_semi.size());
  }

  // Binary garbage
  {
    const char garbage[] = "\x01\x02\xff\xfe\x80\x00\x7f\xc0\xc1";
    run_case("binary_garbage", garbage, sizeof(garbage) - 1);
  }

  // Binary garbage after valid header
  {
    std::string mixed = "data_x\n_t ";
    mixed += std::string("\xff\xfe\x00\x01\x80", 5);
    mixed += "\n";
    run_case("mixed_valid_binary", mixed.c_str(), mixed.size());
  }

  // Numeric edge cases
  run_case("numeric_empty", "");
  run_case("numeric_dot", ".");
  run_case("numeric_question", "?");
  run_case("numeric_plus", "+");
  run_case("numeric_minus", "-");
  run_case("numeric_e", "1e");
  run_case("numeric_su_open", "1.0(");
  run_case("numeric_su_empty", "1.0()");
  run_case("numeric_su_nested", "1.0((3))");
  run_case("numeric_overflow", "1e999");
  run_case("numeric_neg_overflow", "-1e999");
  run_case("numeric_many_digits",
    "123456789012345678901234567890.123456789");

  // data_ with special characters in block name
  run_case("block_name_special", "data_a-b.c_d+e\n");
  run_case("block_name_numeric", "data_123\n");

  // Multiple data_ blocks with same name
  run_case("duplicate_blocks", "data_x\n_a 1\ndata_x\n_b 2\n");

  // Comment-only file
  run_case("comments_only", "# comment 1\n# comment 2\n");

  // Comment after data
  run_case("comment_after", "data_x\n_t v # inline comment\n");

  if (failures == 0) {
    std::printf("All fuzz corpus cases passed.\n");
  } else {
    std::fprintf(stderr, "%d fuzz corpus case(s) FAILED.\n", failures);
  }
  return failures > 0 ? 1 : 0;
}

#endif // XCIF_FUZZ_ENGINE
