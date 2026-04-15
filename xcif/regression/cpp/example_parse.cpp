// cctbx_project/xcif/regression/cpp/example_parse.cpp
//
// Human-readable CIF parser demo.
// Parses a CIF file into a Document, walks its blocks, tag-value pairs,
// loops, and save frames, and prints a structured summary.
//
// This file is intentionally verbose with comments — it is the canonical
// "getting started" example for the xcif parser and numeric conversions.
//
// Build (from the phenix build directory):
//   libtbx.scons -j 4
//
// Run on the bundled example:
//   build/xcif/regression/cpp/example_parse \
//       modules/cctbx_project/xcif/regression/example.cif

#include <cstdio>
#include <cmath>
#include <string>

// ---------------------------------------------------------------------------
// xcif headers.
// ---------------------------------------------------------------------------

// MappedFile: zero-copy read-only memory-mapped I/O (POSIX + Windows).
#include "xcif/mapped_file.h"

// data_model.h: Document / Block / Loop / parse() — the parsed CIF AST.
// All tag names and values are zero-copy string_views into the source buffer.
// The Document owns the buffer (MappedFile or string copy), so the views
// are valid for the Document's lifetime.
#include "xcif/data_model.h"

// numeric.h: as_double(), as_int(), as_double_with_su(), is_null(), etc.
// Free functions that operate on string_view values.
#include "xcif/numeric.h"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Print a string_view safely (it is NOT NUL-terminated).
static void print_sv(const xcif::string_view& sv) {
  std::printf("%.*s", (int)sv.size(), sv.data());
}

// Print a string_view, truncating multi-line content for readability.
static void print_sv_oneline(const xcif::string_view& sv, int max_len = 60) {
  const char* p = sv.data();
  std::size_t len = sv.size();
  // Find first newline.
  std::size_t nl = 0;
  while (nl < len && p[nl] != '\n' && p[nl] != '\r') ++nl;
  std::size_t show = nl;
  if (show > static_cast<std::size_t>(max_len)) show = max_len;
  std::printf("%.*s", (int)show, p);
  if (show < len) std::printf(" ...");
}

// ---------------------------------------------------------------------------
// print_block — recursively prints a Block (used for both data blocks
//               and save frames).
// ---------------------------------------------------------------------------
static void print_block(const xcif::Block& blk, int indent) {
  char pad[32];
  int n = indent * 2;
  if (n > 30) n = 30;
  for (int i = 0; i < n; ++i) pad[i] = ' ';
  pad[n] = '\0';

  // ---- Tag-value pairs ----
  // Block::pairs_ is private, but we can probe individual tags.
  // For a demo, we list a few well-known tags from example.cif.
  // In real code you would iterate via the Python bindings (Step 10).

  // ---- Loops ----
  const std::vector<xcif::Loop>& loops = blk.loops();
  std::printf("%s  Loops: %zu\n", pad, loops.size());

  for (std::size_t li = 0; li < loops.size(); ++li) {
    const xcif::Loop& lp = loops[li];
    std::printf("%s  Loop %zu: %zu tags x %zu rows\n",
                pad, li, lp.width(), lp.length());

    // Print tag names.
    std::printf("%s    Tags:", pad);
    const std::vector<xcif::string_view>& tags = lp.tags();
    for (std::size_t ti = 0; ti < tags.size(); ++ti) {
      std::printf(" ");
      print_sv(tags[ti]);
    }
    std::printf("\n");

    // Print first 3 rows (if any) as a table.
    std::size_t show_rows = lp.length();
    if (show_rows > 3) show_rows = 3;
    for (std::size_t r = 0; r < show_rows; ++r) {
      std::printf("%s    [%zu]", pad, r);
      for (std::size_t c = 0; c < lp.width(); ++c) {
        std::printf("  ");
        print_sv_oneline(lp.value(r, c), 20);
      }
      std::printf("\n");
    }
    if (lp.length() > 3)
      std::printf("%s    ... (%zu more rows)\n", pad, lp.length() - 3);
  }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::fprintf(stderr,
      "Usage: %s <file.cif>\n\n"
      "Parses a CIF file and prints a structured summary.\n\n"
      "  %s modules/cctbx_project/xcif/regression/example.cif\n\n",
      argv[0], argv[0]);
    return 1;
  }

  const char* path = argv[1];

  // ------------------------------------------------------------------
  // Step 1 — Memory-map the file.
  // ------------------------------------------------------------------
  xcif::MappedFile mf(path);
  std::printf("=== xcif parser demo ===\n");
  std::printf("File: %s (%zu bytes)\n\n", mf.path(), mf.size());

  // ------------------------------------------------------------------
  // Step 2 — Parse into a Document.
  //
  // parse() takes a buffer + length.  The Document copies the buffer
  // into an internal vector<char> so all string_views remain valid.
  // For file I/O the MappedFile version (Step 12) will avoid the copy.
  // ------------------------------------------------------------------
  xcif::Document doc = xcif::parse(mf.data(), mf.size(), mf.path());
  std::printf("Blocks: %zu\n", doc.size());

  // ------------------------------------------------------------------
  // Step 3 — Walk the data blocks.
  // ------------------------------------------------------------------
  for (std::size_t bi = 0; bi < doc.size(); ++bi) {
    const xcif::Block& blk = doc[bi];
    std::printf("\n--- Block: ");
    print_sv(blk.name());
    std::printf(" ---\n");

    // ---- Demonstrate tag-value lookup ----
    // find_value() is case-insensitive and returns an empty string_view
    // if the tag is not present.
    std::printf("\n  Tag-value lookups (case-insensitive):\n");

    xcif::string_view cell_a = blk.find_value("_cell.length_a");
    if (!cell_a.empty()) {
      std::printf("    _cell.length_a   = ");
      print_sv(cell_a);
      std::printf("  (as_double: %.3f)\n", xcif::as_double(cell_a));

      // Demonstrate as_double_with_su() — extracts value and standard
      // uncertainty from "50.840(10)" → value=50.840, su=0.010
      std::pair<double, double> val_su = xcif::as_double_with_su(cell_a);
      std::printf("    as_double_with_su: value=%.3f  su=%.3f\n",
                  val_su.first, val_su.second);
    }

    // Show all three cell lengths with SU extraction.
    xcif::string_view cell_b = blk.find_value("_cell.length_b");
    xcif::string_view cell_c = blk.find_value("_cell.length_c");
    if (!cell_b.empty() && !cell_c.empty()) {
      std::pair<double, double> b_su = xcif::as_double_with_su(cell_b);
      std::pair<double, double> c_su = xcif::as_double_with_su(cell_c);
      std::printf("    _cell.length_b   = ");
      print_sv(cell_b);
      std::printf("  (value=%.3f  su=%.3f)\n", b_su.first, b_su.second);
      std::printf("    _cell.length_c   = ");
      print_sv(cell_c);
      std::printf("  (value=%.3f  su=%.3f)\n", c_su.first, c_su.second);
    }

    xcif::string_view sg = blk.find_value("_symmetry.space_group_name_H-M");
    if (!sg.empty()) {
      std::printf("    _symmetry.space_group_name_H-M = ");
      print_sv(sg);
      std::printf("\n");
    }

    xcif::string_view rfree = blk.find_value("_refine.ls_R_factor_R_free");
    if (!rfree.empty()) {
      std::printf("    _refine.ls_R_factor_R_free = ");
      print_sv(rfree);
      double val = xcif::as_double(rfree);
      std::printf("  (as_double: %.4f)\n", val);
    }

    // ---- Demonstrate unknown/inapplicable detection ----
    xcif::string_view reflns = blk.find_value("_refine.ls_number_reflns_obs");
    if (!reflns.empty()) {
      std::printf("    _refine.ls_number_reflns_obs = ");
      print_sv(reflns);
      if (xcif::is_inapplicable(reflns))
        std::printf("  [inapplicable]");
      else if (xcif::is_unknown(reflns))
        std::printf("  [unknown]");
      std::printf("\n");
    }

    xcif::string_view biso = blk.find_value("_refine.B_iso_mean");
    if (!biso.empty()) {
      std::printf("    _refine.B_iso_mean = ");
      print_sv(biso);
      if (xcif::is_null(biso))
        std::printf("  [null → NaN: %s]",
                    std::isnan(xcif::as_double(biso)) ? "yes" : "no");
      std::printf("\n");
    }

    // Demonstrate has_tag() — checking existence without reading value.
    std::printf("    has_tag(_cell.volume)  = %s\n",
                blk.has_tag("_cell.volume") ? "true" : "false");
    std::printf("    has_tag(_nonexistent)  = %s\n",
                blk.has_tag("_nonexistent") ? "true" : "false");

    // ---- Demonstrate loop access ----
    std::printf("\n");
    print_block(blk, 0);

    // ---- Demonstrate find_loop() by tag name ----
    // Use column_index() to look up columns by tag name — never hardcode
    // column positions, since CIF files vary in column order and count.
    // find_loop() accepts a full tag name OR a category prefix.
    const xcif::Loop* atom_loop = blk.find_loop("_atom_site");
    if (atom_loop) {
      std::size_t ci_id   = atom_loop->column_index("_atom_site.id");
      std::size_t ci_atom = atom_loop->column_index("_atom_site.label_atom_id");
      std::size_t ci_x    = atom_loop->column_index("_atom_site.Cartn_x");
      std::size_t ci_y    = atom_loop->column_index("_atom_site.Cartn_y");
      std::size_t ci_z    = atom_loop->column_index("_atom_site.Cartn_z");
      std::size_t ci_b    = atom_loop->column_index("_atom_site.B_iso_or_equiv");

      std::printf("\n  Atom coordinates (from find_loop + column_index):\n");
      std::printf("    %-4s  %-4s  %8s  %8s  %8s  %6s\n",
                  "ID", "Atom", "X", "Y", "Z", "B");
      std::size_t n = atom_loop->length();
      if (n > 5) n = 5;
      for (std::size_t r = 0; r < n; ++r) {
        xcif::string_view id   = atom_loop->value(r, ci_id);
        xcif::string_view atom = atom_loop->value(r, ci_atom);
        xcif::string_view x    = atom_loop->value(r, ci_x);
        xcif::string_view y    = atom_loop->value(r, ci_y);
        xcif::string_view z    = atom_loop->value(r, ci_z);
        xcif::string_view b    = atom_loop->value(r, ci_b);
        std::printf("    %-4.*s  %-4.*s  %8.3f  %8.3f  %8.3f  %6.2f\n",
                    (int)id.size(), id.data(),
                    (int)atom.size(), atom.data(),
                    xcif::as_double(x),
                    xcif::as_double(y),
                    xcif::as_double(z),
                    xcif::as_double(b));
      }
      if (atom_loop->length() > 5)
        std::printf("    ... (%zu more atoms)\n",
                    atom_loop->length() - 5);
    }

    // ---- Demonstrate column extraction ----
    const xcif::Loop* reflns_loop = blk.find_loop("_reflns.d_resolution_high");
    if (reflns_loop) {
      std::printf("\n  Column extraction (_reflns.d_resolution_high):\n");
      std::vector<xcif::string_view> col =
        reflns_loop->column("_reflns.d_resolution_high");
      for (std::size_t i = 0; i < col.size(); ++i) {
        std::printf("    [%zu] ", i);
        print_sv(col[i]);
        std::printf("  → as_double = %.2f\n", xcif::as_double(col[i]));
      }
    }

    // ---- Demonstrate as_int ----
    xcif::string_view sg_num = blk.find_value("_symmetry.Int_Tables_number");
    if (!sg_num.empty()) {
      std::printf("\n  Integer conversion:\n");
      std::printf("    _symmetry.Int_Tables_number = ");
      print_sv(sg_num);
      std::printf("  (as_int: %d)\n", xcif::as_int(sg_num));
    }
  }

  std::printf("\nDone.\n");
  return 0;
}
