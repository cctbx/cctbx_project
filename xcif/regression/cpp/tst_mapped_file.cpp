// cctbx_project/xcif/regression/cpp/tst_mapped_file.cpp
//
// Tests for xcif::MappedFile (mapped_file.h).
// Each test creates a temporary file, maps it, verifies behaviour, then
// removes the file.  Helpers use only standard C/C++ facilities.

#include "test_utils.h"
#include "xcif/mapped_file.h"
#include "xcif/tokenizer.h"

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Write `len` bytes to a new temporary file; return the path.
#ifdef _WIN32
static std::string write_temp(const char* content, std::size_t len) {
  char buf[L_tmpnam];
  if (!std::tmpnam(buf)) {
    throw std::runtime_error("tmpnam failed");
  }
  std::string path(buf);
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) {
    throw std::runtime_error(std::string("cannot create temp file: ") + path);
  }
  if (len > 0) std::fwrite(content, 1, len, f);
  std::fclose(f);
  return path;
}
#else
#include <unistd.h>
static std::string write_temp(const char* content, std::size_t len) {
  char tmpl[] = "/tmp/xcif_test_XXXXXX";
  int fd = mkstemp(tmpl);
  if (fd < 0) {
    throw std::runtime_error("mkstemp failed");
  }
  if (len > 0) {
    ssize_t n = write(fd, content, len);
    if (n < 0 || static_cast<std::size_t>(n) != len) {
      close(fd);
      throw std::runtime_error("write failed");
    }
  }
  close(fd);
  return std::string(tmpl);
}
#endif

static void rm_temp(const std::string& path) {
  std::remove(path.c_str());
}

// Succeeds if `expr` throws any std::exception; fails otherwise.
#define CHECK_THROWS(expr) \
  do { \
    bool _threw = false; \
    try { (void)(expr); } catch (const std::exception&) { _threw = true; } \
    if (!_threw) { \
      std::fprintf(stderr, "FAIL: %s:%d: expected exception from: %s\n", \
                   __FILE__, __LINE__, #expr); \
      ++xcif_test_failures; \
    } \
  } while (0)

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

static void test_map_basic_content() {
  const char* content = "data_test _a 42\n";
  std::string path = write_temp(content, std::strlen(content));
  {
    xcif::MappedFile mf(path.c_str());
    CHECK_EQ(mf.size(), std::strlen(content));
    CHECK(mf.data() != nullptr);
    CHECK(!mf.empty());
    CHECK_EQ(std::string(mf.data(), mf.size()), std::string(content));
  }
  rm_temp(path);
}

static void test_map_empty_file() {
  std::string path = write_temp("", 0);
  {
    xcif::MappedFile mf(path.c_str());
    CHECK_EQ(mf.size(), (std::size_t)0);
    CHECK(mf.data() == nullptr);
    CHECK(mf.empty());
  }
  rm_temp(path);
}

static void test_path_accessor() {
  const char* content = "data_x\n";
  std::string path = write_temp(content, std::strlen(content));
  {
    xcif::MappedFile mf(path.c_str());
    CHECK_EQ(std::string(mf.path()), path);
  }
  rm_temp(path);
}

static void test_string_path_ctor() {
  // Verify the std::string overload delegates correctly.
  const char* content = "data_x _b 99\n";
  std::string path = write_temp(content, std::strlen(content));
  {
    xcif::MappedFile mf(path);   // std::string overload
    CHECK_EQ(mf.size(), std::strlen(content));
    CHECK(mf.data() != nullptr);
  }
  rm_temp(path);
}

static void test_nonexistent_file_throws() {
  CHECK_THROWS(xcif::MappedFile("/xcif_no_such_file_xyz_12345.cif"));
}

static void test_null_path_throws() {
  CHECK_THROWS(xcif::MappedFile(static_cast<const char*>(nullptr)));
}

static void test_empty_path_throws() {
  CHECK_THROWS(xcif::MappedFile(""));
}

static void test_tokenizer_integration() {
  // Map a real CIF snippet and tokenize through the MappedFile layer.
  // Tokens must point into the mapped region (zero-copy).
  const char* content = "data_test\n_cell.length_a 5.0\n_cell.length_b 6.0\n";
  std::string path = write_temp(content, std::strlen(content));
  {
    xcif::MappedFile mf(path.c_str());
    xcif::Tokenizer  tok(mf.data(), mf.size(), mf.path());

    xcif::Token t0 = tok.next();
    CHECK_EQ(t0.type, xcif::TOKEN_BLOCK_HEADER);
    CHECK_EQ(t0.as_str(), std::string("data_test"));

    xcif::Token t1 = tok.next();
    CHECK_EQ(t1.type, xcif::TOKEN_TAG);
    CHECK_EQ(t1.as_str(), std::string("_cell.length_a"));

    xcif::Token t2 = tok.next();
    CHECK_EQ(t2.type, xcif::TOKEN_VALUE);
    CHECK_EQ(t2.as_str(), std::string("5.0"));

    // Token pointers fall inside the mapped region.
    CHECK_GE(t0.ptr, mf.data());
    CHECK_LT(t0.ptr, mf.data() + mf.size());
  }
  rm_temp(path);
}

static void test_move_constructor() {
  const char* content = "data_moved _x 1\n";
  std::string path = write_temp(content, std::strlen(content));
  {
    xcif::MappedFile mf1(path.c_str());
    const char*  ptr = mf1.data();
    std::size_t  sz  = mf1.size();

    xcif::MappedFile mf2(std::move(mf1));

    CHECK_EQ(mf2.data(), ptr);
    CHECK_EQ(mf2.size(), sz);
    CHECK(!mf2.empty());

    // Moved-from must be hollow.
    CHECK(mf1.data() == nullptr);
    CHECK_EQ(mf1.size(), (std::size_t)0);
    CHECK(mf1.empty());
  }
  rm_temp(path);
}

static void test_move_assignment() {
  const char* content = "data_mv2 _y 2\n";
  std::string path = write_temp(content, std::strlen(content));
  {
    xcif::MappedFile mf1(path.c_str());
    const char*  ptr = mf1.data();
    std::size_t  sz  = mf1.size();

    // mf2 holds its own mapping; overwrite it via move assignment.
    xcif::MappedFile mf2(path.c_str());
    mf2 = std::move(mf1);

    CHECK_EQ(mf2.data(), ptr);
    CHECK_EQ(mf2.size(), sz);
    CHECK(mf1.data() == nullptr);
    CHECK_EQ(mf1.size(), (std::size_t)0);
  }
  rm_temp(path);
}

static void test_self_move_assignment_safe() {
  // `if (this != &other)` in operator= must make self-move a no-op.
  // Use a pointer indirection to suppress compiler -Wself-move warnings.
  const char* content = "data_self _z 3\n";
  std::string path = write_temp(content, std::strlen(content));
  {
    xcif::MappedFile  mf(path.c_str());
    std::size_t       expected_size = mf.size();
    const char*       expected_data = mf.data();

    xcif::MappedFile* p = &mf;
    *p = std::move(*p);   // self-move via pointer — must be a no-op

    CHECK_EQ(mf.size(), expected_size);
    CHECK_EQ(mf.data(), expected_data);
  }
  rm_temp(path);
}

static void test_bytes_preserved_exactly() {
  // Verify the mapping is a faithful byte-for-byte copy of the file.
  const char content[] = "data_bin\n;\nline\xC3\xA9\n;\n";  // UTF-8 content
  std::size_t len = sizeof(content) - 1;  // exclude trailing NUL
  std::string path = write_temp(content, len);
  {
    xcif::MappedFile mf(path.c_str());
    CHECK_EQ(mf.size(), len);
    CHECK(std::memcmp(mf.data(), content, len) == 0);
  }
  rm_temp(path);
}

static void test_large_file() {
  // Map a file larger than a typical OS page (4 KB) to exercise multi-page
  // mapping and the MADV_SEQUENTIAL readahead path.
  const std::size_t N = 8192;
  std::string content;
  content.reserve(N);
  for (std::size_t i = 0; i < N; ++i) {
    content.push_back(static_cast<char>('A' + (i % 26)));
  }
  std::string path = write_temp(content.c_str(), content.size());
  {
    xcif::MappedFile mf(path.c_str());
    CHECK_EQ(mf.size(), N);
    CHECK(std::memcmp(mf.data(), content.c_str(), N) == 0);
  }
  rm_temp(path);
}

// ---------------------------------------------------------------------------
// Test registry
// ---------------------------------------------------------------------------

static void run_all_tests() {
  test_map_basic_content();
  test_map_empty_file();
  test_path_accessor();
  test_string_path_ctor();
  test_nonexistent_file_throws();
  test_null_path_throws();
  test_empty_path_throws();
  test_tokenizer_integration();
  test_move_constructor();
  test_move_assignment();
  test_self_move_assignment_safe();
  test_bytes_preserved_exactly();
  test_large_file();
}

XCIF_TEST_MAIN()
