// cctbx_project/xcif/data_model.cpp
#include <xcif/data_model.h>
#include <xcif/tokenizer.h>
#include <xcif/mapped_file.h>
#include <sstream>

namespace xcif {

// ── CifError ───────────────────────────────────────────────────────

static std::string format_error(const std::string& msg, int line, int col,
                                const std::string& source) {
  std::ostringstream os;
  os << source << ":" << line << ":" << col << ": " << msg;
  return os.str();
}

CifError::CifError(const std::string& msg, int line, int col,
                   const std::string& source)
  : std::runtime_error(format_error(msg, line, col, source))
  , line_(line), col_(col) {}

// ── Loop::column ───────────────────────────────────────────────────

std::vector<string_view> Loop::column(const string_view& tag) const {
  ci_map<std::size_t>::const_iterator it = tag_index_.find(tag);
  if (it == tag_index_.end()) return std::vector<string_view>();
  std::size_t ci = it->second;
  std::size_t w = tags_.size();
  std::size_t n = length();
  std::vector<string_view> result;
  result.reserve(n);
  for (std::size_t r = 0; r < n; ++r)
    result.push_back(values_[r * w + ci]);
  return result;
}

// ── Block accessors ────────────────────────────────────────────────

string_view Block::find_value(const string_view& tag) const {
  ci_map<std::size_t>::const_iterator it = pair_index_.find(tag);
  if (it == pair_index_.end()) return string_view();
  return pairs_[it->second].second;
}

const Loop* Block::find_loop(const string_view& tag) const {
  // Fast path: exact tag match (e.g. "_atom_site.Cartn_x").
  ci_map<std::size_t>::const_iterator it = loop_tag_index_.find(tag);
  if (it != loop_tag_index_.end()) return &loops_[it->second];

  // Category prefix match (e.g. "_atom_site" or "_atom_site.").
  // A query is treated as a prefix if it contains no '.' (after the
  // leading '_'), or ends with '.'.  We scan loop_tag_index_ for the
  // first tag whose category matches.
  bool has_dot = false;
  for (std::size_t i = 1; i < tag.size(); ++i) {
    if (tag.data()[i] == '.') { has_dot = true; break; }
  }
  if (has_dot && tag.data()[tag.size() - 1] != '.') return 0;

  // Build the prefix to match: ensure it ends with '.'.
  // e.g. "_atom_site" → "_atom_site.", "_atom_site." → "_atom_site."
  std::string prefix;
  if (!has_dot) {
    prefix.assign(tag.data(), tag.size());
    prefix += '.';
  } else {
    prefix.assign(tag.data(), tag.size());
  }

  CIEqual eq;
  string_view prefix_sv(prefix.data(), prefix.size());
  for (it = loop_tag_index_.begin(); it != loop_tag_index_.end(); ++it) {
    const string_view& key = it->first;
    if (key.size() > prefix_sv.size()) {
      string_view key_prefix(key.data(), prefix_sv.size());
      if (eq(key_prefix, prefix_sv))
        return &loops_[it->second];
    }
  }
  return 0;
}

const Block* Block::find_save_frame(const string_view& nm) const {
  ci_map<std::size_t>::const_iterator it = save_frame_index_.find(nm);
  if (it == save_frame_index_.end()) return 0;
  return &save_frames_[it->second];
}

// ── Document accessor ──────────────────────────────────────────────

const Block* Document::find_block(const string_view& nm) const {
  ci_map<std::size_t>::const_iterator it = block_index_.find(nm);
  if (it == block_index_.end()) return 0;
  return &blocks_[it->second];
}

// ── Parser ─────────────────────────────────────────────────────────

namespace detail {

class Parser {
public:
  Parser(const char* data, std::size_t len, const char* src,
         bool strict = true)
    : tok_(data, len, src), source_(src), strict_(strict) { advance(); }

  void run(Document& doc) {
    while (cur_.type != TOKEN_EOF) {
      if (cur_.type == TOKEN_BLOCK_HEADER) {
        doc.blocks_.push_back(Block());
        Block& blk = doc.blocks_.back();
        // Two kinds of TOKEN_BLOCK_HEADER tokens come from the tokenizer:
        //   "data_X"  (any length >=5) — block name is X
        //   "global_" (exactly 7 chars) — block name is the full "global_"
        // For "global_", point the string_view at a static literal so the
        // name survives Document copies/moves and matches the non-strict
        // synthesized block's storage strategy.
        if (cur_.len == 7 &&
            (cur_.ptr[0] == 'g' || cur_.ptr[0] == 'G')) {
          static const char kGlobalName[] = "global_";
          blk.name_ = string_view(kGlobalName, 7);
        } else {
          blk.name_ = string_view(cur_.ptr + 5, cur_.len - 5);
        }
        doc.block_index_[blk.name_] = doc.blocks_.size() - 1;
        advance();
        parse_content(blk, false);
      } else if (!strict_) {
        // Non-strict: accumulate pre-block-header content into an
        // implicit block named "global_". Created on demand; reused if
        // content appears both before any data_ block and again between
        // two data_ blocks (rare but legal under ucif-compat semantics).
        Block& global_blk = find_or_create_global_block(doc);
        parse_content(global_blk, false);
      } else {
        error("data outside of a data block");
      }
    }
  }

  Block& find_or_create_global_block(Document& doc) {
    // The name "global_" lives in static storage so the string_view has
    // stable backing regardless of how the Document is copied or moved.
    // (Storing it in a std::string member of Document triggers SSO: the
    // string's buffer moves with the object, dangling any string_view
    // into it.)
    static const char kGlobalName[] = "global_";
    static const std::size_t kGlobalLen = sizeof(kGlobalName) - 1;
    string_view name(kGlobalName, kGlobalLen);
    ci_map<std::size_t>::iterator it = doc.block_index_.find(name);
    if (it != doc.block_index_.end()) {
      return doc.blocks_[it->second];
    }
    doc.blocks_.push_back(Block());
    Block& blk = doc.blocks_.back();
    blk.name_ = name;
    doc.block_index_[name] = doc.blocks_.size() - 1;
    return blk;
  }

private:
  Tokenizer tok_;
  Token cur_;
  std::string source_;
  bool strict_;

  void advance() { cur_ = tok_.next(); }

  void error(const std::string& msg) {
    throw CifError(msg, cur_.line, cur_.col, source_);
  }

  void error_at(const std::string& msg, int line, int col) {
    throw CifError(msg, line, col, source_);
  }

  string_view sv() const { return string_view(cur_.ptr, cur_.len); }

  bool tag_exists(const Block& blk, const string_view& tag) {
    return blk.pair_index_.find(tag) != blk.pair_index_.end() ||
           blk.loop_tag_index_.find(tag) != blk.loop_tag_index_.end();
  }

  // Returns true if ended by SAVE_END
  bool parse_content(Block& blk, bool in_save) {
    while (cur_.type != TOKEN_EOF && cur_.type != TOKEN_BLOCK_HEADER) {
      switch (cur_.type) {
        case TOKEN_TAG:
          parse_pair(blk);
          break;
        case TOKEN_LOOP:
          parse_loop(blk);
          break;
        case TOKEN_SAVE_HEADER: {
          if (in_save) error("nested save frames are not allowed");
          string_view sname(cur_.ptr + 5, cur_.len - 5);
          blk.save_frames_.push_back(Block());
          Block& sf = blk.save_frames_.back();
          sf.name_ = sname;
          blk.save_frame_index_[sname] = blk.save_frames_.size() - 1;
          advance();
          if (!parse_content(sf, true))
            error("unclosed save frame");
          break;
        }
        case TOKEN_SAVE_END:
          if (!in_save) error("save_ without matching save frame");
          advance();
          return true;
        case TOKEN_VALUE:
          error("data outside of a data block");
          break;
        default:
          error("unexpected token");
          break;
      }
    }
    return false;
  }

  void parse_pair(Block& blk) {
    string_view tag = sv();
    int tline = cur_.line, tcol = cur_.col;
    if (tag_exists(blk, tag))
      error("duplicate tag '" + std::string(tag) + "'");
    advance();
    if (cur_.type != TOKEN_VALUE)
      error_at("missing value for tag '" + std::string(tag) + "'",
               tline, tcol);
    string_view val = sv();
    blk.pairs_.push_back(std::make_pair(tag, val));
    blk.pair_index_[tag] = blk.pairs_.size() - 1;
    advance();
  }

  void parse_loop(Block& blk) {
    advance(); // consume loop_
    Loop lp;
    // Read tags
    while (cur_.type == TOKEN_TAG) {
      string_view tag = sv();
      if (tag_exists(blk, tag) ||
          lp.tag_index_.find(tag) != lp.tag_index_.end())
        error("duplicate tag '" + std::string(tag) + "'");
      lp.tag_index_[tag] = lp.tags_.size();
      lp.tags_.push_back(tag);
      advance();
    }
    if (lp.tags_.empty())
      error("loop_ must be followed by at least one tag");
    // Read values
    while (cur_.type == TOKEN_VALUE) {
      lp.values_.push_back(sv());
      advance();
    }
    if (!lp.values_.empty() && lp.values_.size() % lp.tags_.size() != 0)
      error("loop value count is not a multiple of tag count");
    // Register loop tags in block index
    std::size_t li = blk.loops_.size();
    for (std::size_t i = 0; i < lp.tags_.size(); ++i)
      blk.loop_tag_index_[lp.tags_[i]] = li;
    blk.loops_.push_back(std::move(lp));
  }
};

} // namespace detail

// ── parse() ────────────────────────────────────────────────────────

Document parse(const char* data, std::size_t length, const char* source,
               bool strict) {
  Document doc;
  doc.owned_buf_.assign(data, data + length);
  detail::Parser parser(doc.owned_buf_.data(), doc.owned_buf_.size(),
                        source, strict);
  parser.run(doc);
  return doc;
}

Document parse_file(const char* path, bool strict) {
  Document doc;
  doc.mapped_file_ = MappedFile(path);
  detail::Parser parser(doc.mapped_file_.data(), doc.mapped_file_.size(),
                        path, strict);
  parser.run(doc);
  return doc;
}

} // namespace xcif
