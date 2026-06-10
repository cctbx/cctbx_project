// cctbx_project/xcif/regression/cpp/example_tokenize.cpp
//
// Human-readable CIF tokenizer demo.
// Reads a CIF file via MappedFile, drives the Tokenizer, and prints every
// token with its type, source location, and text.
//
// This file is intentionally verbose with comments — it is the canonical
// "getting started" example for xcif.
//
// Build (from the phenix build directory):
//   libtbx.scons -j 4
//
// Run on the bundled example:
//   build/xcif/regression/cpp/example_tokenize \
//       modules/cctbx_project/xcif/regression/example.cif

#include <cstdio>
#include <cstring>

// ---------------------------------------------------------------------------
// xcif headers.
//
// The SCons build system adds  xcif/include/  to the compiler's include path,
// so we write "xcif/foo.h" rather than a relative path.
// ---------------------------------------------------------------------------

// MappedFile: RAII wrapper for read-only memory-mapped I/O.
//   POSIX   -> mmap(2) + MADV_SEQUENTIAL
//   Windows -> CreateFileMapping / MapViewOfFile
// The file is zero-copied: no bytes are read into heap memory.
#include "xcif/mapped_file.h"

// Tokenizer: walks the mapped buffer and returns one Token per call to
// next().  Tokens are zero-copy views (ptr + len) into the mapping.
// The mapping MUST outlive the Tokenizer and any tokens it produces.
#include "xcif/tokenizer.h"

// ---------------------------------------------------------------------------
// token_type_name — printable label for each TokenType enum value
// ---------------------------------------------------------------------------
static const char* token_type_name(xcif::TokenType t) {
  switch (t) {
    case xcif::TOKEN_EOF:          return "EOF";
    case xcif::TOKEN_TAG:          return "TAG";
    case xcif::TOKEN_VALUE:        return "VALUE";
    case xcif::TOKEN_BLOCK_HEADER: return "BLOCK_HEADER";
    case xcif::TOKEN_LOOP:         return "LOOP";
    case xcif::TOKEN_SAVE_HEADER:  return "SAVE_HEADER";
    case xcif::TOKEN_SAVE_END:     return "SAVE_END";
    default:                       return "???";
  }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // ------------------------------------------------------------------
  // Argument handling
  // ------------------------------------------------------------------
  if (argc < 2) {
    std::fprintf(stderr,
      "Usage: %s <file.cif>\n"
      "\n"
      "Tokenizes a CIF file and prints every token.\n"
      "Run on the bundled example file:\n"
      "\n"
      "  %s modules/cctbx_project/xcif/regression/example.cif\n"
      "\n",
      argv[0], argv[0]);
    return 1;
  }

  const char* path = argv[1];

  // ------------------------------------------------------------------
  // Step 1 — Open the file with MappedFile.
  //
  // This call maps the file read-only into the process address space.
  // No bytes are copied to heap.  On failure it throws std::runtime_error
  // with a descriptive message (OS error string included).
  //
  // mf is RAII: the mapping is released when mf goes out of scope.
  // ------------------------------------------------------------------
  xcif::MappedFile mf(path);

  std::printf("=== xcif tokenizer demo ===\n");
  std::printf("File : %s\n", mf.path());
  std::printf("Size : %zu bytes\n", mf.size());

  if (mf.empty()) {
    std::printf("(empty file — nothing to tokenize)\n");
    return 0;
  }

  // ------------------------------------------------------------------
  // Step 2 — Construct the Tokenizer.
  //
  // Tokenizer(const MappedFile&) is a convenience overload that
  // forwards mf.data(), mf.size(), and mf.path() to the primary
  // Tokenizer(const char* data, std::size_t len, const char* source)
  // constructor.
  //
  // The tokenizer keeps a raw pointer into the mapped region.
  // mf must therefore outlive tok (and all Tokens tok produces).
  // ------------------------------------------------------------------
  xcif::Tokenizer tok(mf);

  // ------------------------------------------------------------------
  // Step 3 — Consume tokens.
  //
  // tok.next() advances through the buffer and returns the next Token.
  // It returns TOKEN_EOF once the file is exhausted, and continues to
  // return TOKEN_EOF on every subsequent call (safe to call in a loop).
  //
  // Each Token carries:
  //   .type   TokenType enum (TAG, VALUE, BLOCK_HEADER, LOOP, …)
  //   .ptr    pointer into the mapped buffer — NOT NUL-terminated
  //   .len    byte length of the token text
  //   .line   1-based line number in the source file
  //   .col    1-based column number
  //
  // t.as_str() returns a std::string copy — convenient but allocates.
  // For display we use printf's "%.*s" width specifier to avoid copying.
  // ------------------------------------------------------------------
  std::printf("\n");
  std::printf("%-15s  %5s  %4s  %s\n", "TYPE", "LINE", "COL", "TEXT");
  std::printf("%-15s  %5s  %4s  %s\n",
              "---------------", "-----", "----",
              "-----------------------------------------------");

  // Count tokens by type for the summary.
  int counts[7] = {};   // TOKEN_EOF=0 .. TOKEN_SAVE_END=6
  int total = 0;

  for (;;) {
    xcif::Token t = tok.next();

    if (t.type == xcif::TOKEN_EOF)
      break;

    // Multi-line values (semicolon text fields) can be many lines long.
    // Display only the first line and append " ..." if truncated so the
    // table stays readable.
    const char* text_end = t.ptr + t.len;
    const char* first_nl = t.ptr;
    while (first_nl < text_end && *first_nl != '\n' && *first_nl != '\r')
      ++first_nl;

    std::size_t display_len = static_cast<std::size_t>(first_nl - t.ptr);
    bool multiline = (display_len < t.len);

    std::printf("%-15s  %5d  %4d  %.*s%s\n",
                token_type_name(t.type),
                t.line,
                t.col,
                (int)display_len, t.ptr,
                multiline ? " ..." : "");

    if (t.type >= 0 && t.type <= 6)
      counts[t.type]++;
    ++total;
  }

  // ------------------------------------------------------------------
  // Summary
  // ------------------------------------------------------------------
  std::printf("\n--- %d tokens total ---\n", total);
  for (int i = 1; i <= 6; ++i) {   // skip EOF (i=0)
    if (counts[i] > 0) {
      std::printf("  %-15s : %d\n",
                  token_type_name(static_cast<xcif::TokenType>(i)),
                  counts[i]);
    }
  }

  return 0;
}
