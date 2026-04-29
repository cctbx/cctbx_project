// cctbx_project/xcif/include/xcif/mapped_file.h
#ifndef XCIF_MAPPED_FILE_H
#define XCIF_MAPPED_FILE_H

#include <cstddef>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// Platform detection
// ---------------------------------------------------------------------------
#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#else
#  include <sys/types.h>
#  include <sys/stat.h>
#  include <sys/mman.h>
#  include <fcntl.h>
#  include <unistd.h>
#  include <cerrno>
#  include <cstring>
#endif

namespace xcif {

// ---------------------------------------------------------------------------
// MappedFile — RAII wrapper for memory-mapped read-only file access.
//
// Usage:
//   xcif::MappedFile mf("model.cif");
//   xcif::Tokenizer  tok(mf.data(), mf.size(), mf.path());
//
// The mapping is held for the lifetime of the MappedFile object.  The
// Tokenizer (and any Tokens it produces) point into this mapping, so the
// MappedFile must outlive them.
//
// Move-only (non-copyable).  C++14-compatible.
// ---------------------------------------------------------------------------
class MappedFile {
 public:
  // Default constructor: creates an empty (unmapped) MappedFile.
  MappedFile() noexcept
    : data_(nullptr),
      size_(0)
#ifdef _WIN32
      , file_handle_(INVALID_HANDLE_VALUE)
      , mapping_handle_(nullptr)
#endif
  {}

  // Map the file at `path` read-only into the process address space.
  // Throws std::runtime_error on failure.
  explicit MappedFile(const char* path)
    : path_(path ? path : ""),
      data_(nullptr),
      size_(0)
#ifdef _WIN32
      , file_handle_(INVALID_HANDLE_VALUE)
      , mapping_handle_(nullptr)
#endif
  {
    if (!path || path[0] == '\0') {
      throw std::runtime_error("xcif::MappedFile: empty file path");
    }
    open_and_map(path);
  }

  explicit MappedFile(const std::string& path)
    : MappedFile(path.c_str())
  {}

  ~MappedFile() { close_mapping(); }

  // Move semantics
  MappedFile(MappedFile&& other) noexcept
    : path_(std::move(other.path_)),
      data_(other.data_),
      size_(other.size_)
#ifdef _WIN32
      , file_handle_(other.file_handle_)
      , mapping_handle_(other.mapping_handle_)
#endif
  {
    other.data_ = nullptr;
    other.size_ = 0;
#ifdef _WIN32
    other.file_handle_    = INVALID_HANDLE_VALUE;
    other.mapping_handle_ = nullptr;
#endif
  }

  MappedFile& operator=(MappedFile&& other) noexcept {
    if (this != &other) {
      close_mapping();
      path_ = std::move(other.path_);
      data_ = other.data_;
      size_ = other.size_;
#ifdef _WIN32
      file_handle_         = other.file_handle_;
      mapping_handle_      = other.mapping_handle_;
      other.file_handle_   = INVALID_HANDLE_VALUE;
      other.mapping_handle_= nullptr;
#endif
      other.data_ = nullptr;
      other.size_ = 0;
    }
    return *this;
  }

  // Non-copyable
  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;

  // Accessors
  const char* data() const noexcept { return data_; }
  std::size_t size() const noexcept { return size_; }
  const char* path() const noexcept { return path_.c_str(); }
  bool        empty() const noexcept { return size_ == 0; }

 private:
  std::string path_;
  const char* data_;
  std::size_t size_;

#ifdef _WIN32
  HANDLE file_handle_;
  HANDLE mapping_handle_;

  void open_and_map(const char* path) {
    file_handle_ = CreateFileA(
        path,
        GENERIC_READ,
        FILE_SHARE_READ,
        nullptr,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
        nullptr);

    if (file_handle_ == INVALID_HANDLE_VALUE) {
      throw std::runtime_error(
          std::string("xcif::MappedFile: cannot open file: ") + path);
    }

    LARGE_INTEGER li;
    if (!GetFileSizeEx(file_handle_, &li)) {
      CloseHandle(file_handle_);
      file_handle_ = INVALID_HANDLE_VALUE;
      throw std::runtime_error(
          std::string("xcif::MappedFile: cannot get file size: ") + path);
    }

    size_ = static_cast<std::size_t>(li.QuadPart);

    // An empty file cannot be memory-mapped on Windows.
    // Return a valid object with data_==nullptr and size_==0.
    if (size_ == 0) {
      CloseHandle(file_handle_);
      file_handle_ = INVALID_HANDLE_VALUE;
      return;
    }

    mapping_handle_ = CreateFileMappingA(
        file_handle_,
        nullptr,
        PAGE_READONLY,
        0, 0,       // map entire file
        nullptr);

    if (mapping_handle_ == nullptr) {
      CloseHandle(file_handle_);
      file_handle_ = INVALID_HANDLE_VALUE;
      throw std::runtime_error(
          std::string("xcif::MappedFile: CreateFileMapping failed: ") + path);
    }

    void* base = MapViewOfFile(mapping_handle_, FILE_MAP_READ, 0, 0, 0);
    if (base == nullptr) {
      CloseHandle(mapping_handle_);
      CloseHandle(file_handle_);
      mapping_handle_ = nullptr;
      file_handle_    = INVALID_HANDLE_VALUE;
      throw std::runtime_error(
          std::string("xcif::MappedFile: MapViewOfFile failed: ") + path);
    }

    data_ = static_cast<const char*>(base);
  }

  void close_mapping() {
    if (data_) {
      UnmapViewOfFile(static_cast<const void*>(data_));
      data_ = nullptr;
    }
    if (mapping_handle_) {
      CloseHandle(mapping_handle_);
      mapping_handle_ = nullptr;
    }
    if (file_handle_ != INVALID_HANDLE_VALUE) {
      CloseHandle(file_handle_);
      file_handle_ = INVALID_HANDLE_VALUE;
    }
    size_ = 0;
  }

#else // POSIX

  void open_and_map(const char* path) {
    int fd = ::open(path, O_RDONLY);
    if (fd == -1) {
      throw std::runtime_error(
          std::string("xcif::MappedFile: cannot open file: ") + path +
          " (" + std::strerror(errno) + ")");
    }

    struct stat st;
    if (::fstat(fd, &st) == -1) {
      int err = errno;
      ::close(fd);
      throw std::runtime_error(
          std::string("xcif::MappedFile: fstat failed: ") + path +
          " (" + std::strerror(err) + ")");
    }

    if (st.st_size < 0) {
      ::close(fd);
      throw std::runtime_error(
          std::string("xcif::MappedFile: file reports negative size: ") + path);
    }
    size_ = static_cast<std::size_t>(st.st_size);

    if (size_ == 0) {
      // mmap(2) with length 0 is undefined; just return empty mapping.
      ::close(fd);
      return;
    }

    void* base = ::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd, 0);
    int mmap_errno = errno;

    // The fd can be closed immediately after a successful mmap.
    ::close(fd);

    if (base == MAP_FAILED) {
      size_ = 0;
      throw std::runtime_error(
          std::string("xcif::MappedFile: mmap failed: ") + path +
          " (" + std::strerror(mmap_errno) + ")");
    }

    // Advise the kernel we will read sequentially — this triggers
    // aggressive readahead which substantially improves throughput on
    // large mmCIF files (often 50-200 MB).
    // MADV_SEQUENTIAL is POSIX.1-2001 but not universally defined
    // (some BSDs, older kernels), so guard against its absence.
#ifdef MADV_SEQUENTIAL
    ::madvise(base, size_, MADV_SEQUENTIAL);
#endif

    data_ = static_cast<const char*>(base);
  }

  void close_mapping() {
    if (data_ && size_ > 0) {
      ::munmap(const_cast<void*>(static_cast<const void*>(data_)), size_);
    }
    data_ = nullptr;
    size_ = 0;
  }

#endif // _WIN32 / POSIX
};

} // namespace xcif

#endif // XCIF_MAPPED_FILE_H
