// cctbx_project/xcif/include/xcif/string_view.h
//
// Minimal string_view for C++14.  When the project moves to C++17 the
// preprocessor guard below collapses this to a type alias for
// std::string_view — zero code changes required at call sites.
#ifndef XCIF_STRING_VIEW_H
#define XCIF_STRING_VIEW_H

#if __cplusplus >= 201703L

#include <string_view>
namespace xcif { using string_view = std::string_view; }

#else // C++14 implementation

#include <algorithm>
#include <cstring>
#include <ostream>
#include <string>

namespace xcif {

class string_view {
public:
  typedef const char*  const_iterator;
  typedef std::size_t  size_type;
  static const size_type npos = static_cast<size_type>(-1);

  // Constructors
  string_view() : ptr_(0), len_(0) {}
  string_view(const char* s, size_type n) : ptr_(s), len_(n) {}
  string_view(const char* s) : ptr_(s), len_(s ? std::strlen(s) : 0) {}
  string_view(const std::string& s) : ptr_(s.data()), len_(s.size()) {}

  // Element access
  const char* data()   const { return ptr_; }
  size_type   size()   const { return len_; }
  size_type   length() const { return len_; }
  bool        empty()  const { return len_ == 0; }
  char operator[](size_type i) const { return ptr_[i]; }
  char front() const { return ptr_[0]; }
  char back()  const { return ptr_[len_ - 1]; }

  // Iterators
  const_iterator begin() const { return ptr_; }
  const_iterator end()   const { return ptr_ + len_; }

  // Modifiers
  void remove_prefix(size_type n) { ptr_ += n; len_ -= n; }
  void remove_suffix(size_type n) { len_ -= n; }

  // Operations
  string_view substr(size_type pos, size_type count = npos) const {
    if (pos > len_) pos = len_;
    if (count > len_ - pos) count = len_ - pos;
    return string_view(ptr_ + pos, count);
  }

  size_type find(char c, size_type pos = 0) const {
    for (size_type i = pos; i < len_; ++i) {
      if (ptr_[i] == c) return i;
    }
    return npos;
  }

  int compare(const string_view& o) const {
    size_type n = std::min(len_, o.len_);
    int r = (n == 0) ? 0 : std::memcmp(ptr_, o.ptr_, n);
    if (r != 0) return r;
    return (len_ < o.len_) ? -1 : (len_ > o.len_) ? 1 : 0;
  }

  // Conversion — matches std::string(std::string_view) in C++17
  explicit operator std::string() const { return std::string(ptr_, len_); }

private:
  const char* ptr_;
  size_type   len_;
};

// ── Comparison operators ───────────────────────────────────────────

inline bool operator==(const string_view& a, const string_view& b) {
  return a.size() == b.size() &&
         (a.size() == 0 || std::memcmp(a.data(), b.data(), a.size()) == 0);
}
inline bool operator!=(const string_view& a, const string_view& b) {
  return !(a == b);
}
inline bool operator<(const string_view& a, const string_view& b) {
  return a.compare(b) < 0;
}

// vs const char*
inline bool operator==(const string_view& a, const char* b) {
  return a == string_view(b);
}
inline bool operator==(const char* a, const string_view& b) {
  return string_view(a) == b;
}
inline bool operator!=(const string_view& a, const char* b) {
  return !(a == b);
}
inline bool operator!=(const char* a, const string_view& b) {
  return !(a == b);
}

// vs std::string
inline bool operator==(const string_view& a, const std::string& b) {
  return a == string_view(b);
}
inline bool operator==(const std::string& a, const string_view& b) {
  return string_view(a) == b;
}
inline bool operator!=(const string_view& a, const std::string& b) {
  return !(a == b);
}
inline bool operator!=(const std::string& a, const string_view& b) {
  return !(a == b);
}

// Stream output
inline std::ostream& operator<<(std::ostream& os, const string_view& sv) {
  if (sv.size() > 0) os.write(sv.data(), sv.size());
  return os;
}

} // namespace xcif

#endif // __cplusplus check
#endif // XCIF_STRING_VIEW_H
