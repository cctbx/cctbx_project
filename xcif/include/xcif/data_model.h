// cctbx_project/xcif/include/xcif/data_model.h
#ifndef XCIF_DATA_MODEL_H
#define XCIF_DATA_MODEL_H

#include <xcif/string_view.h>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace xcif {

namespace detail { class Parser; }

// ── Case-insensitive hash / equal for string_view ──────────────────

struct CIHash {
  std::size_t operator()(const string_view& sv) const {
    std::size_t h = 14695981039346656037ULL;
    for (std::size_t i = 0; i < sv.size(); ++i) {
      unsigned char c = static_cast<unsigned char>(sv[i]);
      if (c >= 'A' && c <= 'Z') c += 32;
      h ^= c;
      h *= 1099511628211ULL;
    }
    return h;
  }
};

struct CIEqual {
  bool operator()(const string_view& a, const string_view& b) const {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
      unsigned char ca = static_cast<unsigned char>(a[i]);
      unsigned char cb = static_cast<unsigned char>(b[i]);
      if (ca >= 'A' && ca <= 'Z') ca += 32;
      if (cb >= 'A' && cb <= 'Z') cb += 32;
      if (ca != cb) return false;
    }
    return true;
  }
};

template <typename T>
using ci_map = std::unordered_map<string_view, T, CIHash, CIEqual>;

// ── CifError ───────────────────────────────────────────────────────

class CifError : public std::runtime_error {
public:
  CifError(const std::string& msg, int line, int col,
           const std::string& source);
  int line() const { return line_; }
  int col()  const { return col_; }
private:
  int line_;
  int col_;
};

// ── Loop ───────────────────────────────────────────────────────────

class Loop {
public:
  std::size_t width()  const { return tags_.size(); }
  std::size_t length() const {
    return tags_.empty() ? 0 : values_.size() / tags_.size();
  }
  const std::vector<string_view>& tags() const { return tags_; }
  bool has_tag(const string_view& tag) const {
    return tag_index_.find(tag) != tag_index_.end();
  }
  // Returns the column index for a tag, or (size_t)-1 if not found.
  std::size_t column_index(const string_view& tag) const {
    ci_map<std::size_t>::const_iterator it = tag_index_.find(tag);
    return it != tag_index_.end() ? it->second : std::size_t(-1);
  }
  string_view value(std::size_t row, std::size_t col) const {
    return values_[row * tags_.size() + col];
  }
  std::vector<string_view> column(const string_view& tag) const;

private:
  friend class detail::Parser;
  std::vector<string_view> tags_;
  std::vector<string_view> values_;
  ci_map<std::size_t> tag_index_;
};

// ── Block ──────────────────────────────────────────────────────────

class Block {
public:
  string_view name() const { return name_; }
  bool has_tag(const string_view& tag) const {
    return pair_index_.find(tag) != pair_index_.end();
  }
  string_view find_value(const string_view& tag) const;
  const Loop* find_loop(const string_view& tag) const;
  const std::vector<Loop>& loops() const { return loops_; }
  const Block* find_save_frame(const string_view& name) const;
  const std::vector<std::pair<string_view, string_view>>& pairs() const {
    return pairs_;
  }

private:
  friend class detail::Parser;
  string_view name_;
  std::vector<std::pair<string_view, string_view>> pairs_;
  ci_map<std::size_t> pair_index_;
  std::vector<Loop> loops_;
  ci_map<std::size_t> loop_tag_index_;
  std::vector<Block> save_frames_;
  ci_map<std::size_t> save_frame_index_;
};

// ── Document ───────────────────────────────────────────────────────

class Document {
public:
  std::size_t size() const { return blocks_.size(); }
  const Block& operator[](std::size_t i) const { return blocks_[i]; }
  const Block* find_block(const string_view& name) const;

private:
  friend class detail::Parser;
  friend Document parse(const char*, std::size_t, const char*);
  std::vector<char> owned_buf_;
  std::vector<Block> blocks_;
  ci_map<std::size_t> block_index_;
};

// ── Free functions ─────────────────────────────────────────────────

Document parse(const char* data, std::size_t length,
               const char* source = "<input>");

inline Document parse(const char* data,
                      const char* source = "<input>") {
  return parse(data, std::strlen(data), source);
}

} // namespace xcif
#endif // XCIF_DATA_MODEL_H
