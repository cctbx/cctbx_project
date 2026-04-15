// cctbx_project/xcif/tokenizer.cpp
#include "xcif/tokenizer.h"
#include "xcif/mapped_file.h"
#include <cstring>

namespace xcif {

namespace {

bool is_whitespace(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

bool is_ordinary(char c) {
  // Characters that may appear in an unquoted value / tag / keyword.
  // Excludes whitespace, quotes, semicolons at line-start (handled separately).
  return c != '\0' && !is_whitespace(c);
}

// Case-insensitive prefix match (CIF keywords are case-insensitive).
bool iprefix(const char* buf, std::size_t buflen,
             const char* prefix, std::size_t preflen) {
  if (buflen < preflen) return false;
  for (std::size_t i = 0; i < preflen; ++i) {
    char a = buf[i];
    char b = prefix[i];
    if (a >= 'A' && a <= 'Z') a += 32;
    if (b >= 'A' && b <= 'Z') b += 32;
    if (a != b) return false;
  }
  return true;
}

} // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

Tokenizer::Tokenizer(const char* data, std::size_t length,
                     const char* source_name)
  : start_(data), cur_(data), end_(data + length),
    source_name_(source_name), line_(1), col_(1)
{
  // Skip UTF-8 BOM (\xEF\xBB\xBF) if present
  if (length >= 3
      && (unsigned char)data[0] == 0xEF
      && (unsigned char)data[1] == 0xBB
      && (unsigned char)data[2] == 0xBF) {
    cur_ += 3;
    col_ += 3;
  }
}

// ---------------------------------------------------------------------------
// MappedFile constructor — delegates to the (data, length, name) overload
// ---------------------------------------------------------------------------
Tokenizer::Tokenizer(const MappedFile& mf)
  : Tokenizer(mf.data(),
              mf.size(),
              mf.path())
{}

Tokenizer::Tokenizer(const char* data, const char* source_name)
  : Tokenizer(data, std::strlen(data), source_name)
{}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

void Tokenizer::advance_line() {
  // Called after consuming a '\n' (or '\r' not followed by '\n').
  ++line_;
  col_ = 1;
}

Token Tokenizer::make_token(TokenType type, const char* ptr,
                            std::size_t len, int line, int col) const {
  Token t;
  t.type = type;
  t.ptr  = ptr;
  t.len  = len;
  t.line = line;
  t.col  = col;
  return t;
}

void Tokenizer::skip_whitespace_and_comments() {
  while (!at_end()) {
    char c = *cur_;
    if (c == '#') {
      // Comment: skip to end of line
      while (!at_end() && *cur_ != '\n' && *cur_ != '\r') { ++cur_; }
    } else if (c == '\r') {
      ++cur_;
      if (!at_end() && *cur_ == '\n') ++cur_; // consume \r\n as one
      advance_line();
    } else if (c == '\n') {
      ++cur_;
      advance_line();
    } else if (c == ' ' || c == '\t') {
      ++cur_;
      ++col_;
    } else {
      break;
    }
  }
}

// ---------------------------------------------------------------------------
// Quoted string:  'xxx' or "xxx"
// ---------------------------------------------------------------------------
Token Tokenizer::read_quoted(char delim, int tok_line, int tok_col) {
  ++cur_; ++col_; // consume opening delimiter
  const char* start = cur_;
  while (!at_end()) {
    char c = *cur_;
    if (c == delim) {
      // Check for CIF2 triple-quote: we should not be here for triple-quoted
      // (handled separately), so a lone delimiter ends the string.
      std::size_t len = static_cast<std::size_t>(cur_ - start);
      ++cur_; ++col_; // consume closing delimiter
      return make_token(TOKEN_VALUE, start, len, tok_line, tok_col);
    } else if (c == '\r') {
      ++cur_;
      if (!at_end() && *cur_ == '\n') ++cur_;
      advance_line();
    } else if (c == '\n') {
      ++cur_;
      advance_line();
    } else {
      ++cur_; ++col_;
    }
  }
  // Unterminated quoted string — return what we have (error path)
  return make_token(TOKEN_VALUE, start,
                    static_cast<std::size_t>(cur_ - start),
                    tok_line, tok_col);
}

// ---------------------------------------------------------------------------
// CIF2 triple-quoted string:  """xxx"""  or  '''xxx'''
// ---------------------------------------------------------------------------
Token Tokenizer::read_triple_quoted(char delim, int tok_line, int tok_col) {
  cur_ += 3; col_ += 3; // consume opening triple delimiter
  const char* start = cur_;
  while (!at_end()) {
    char c = *cur_;
    if (c == delim && (cur_ + 2 < end_) && cur_[1] == delim && cur_[2] == delim) {
      std::size_t len = static_cast<std::size_t>(cur_ - start);
      cur_ += 3; col_ += 3;
      return make_token(TOKEN_VALUE, start, len, tok_line, tok_col);
    } else if (c == '\r') {
      ++cur_;
      if (!at_end() && *cur_ == '\n') ++cur_;
      advance_line();
    } else if (c == '\n') {
      ++cur_;
      advance_line();
    } else {
      ++cur_; ++col_;
    }
  }
  return make_token(TOKEN_VALUE, start,
                    static_cast<std::size_t>(cur_ - start),
                    tok_line, tok_col);
}

// ---------------------------------------------------------------------------
// Semicolon text field:  \n;content\n;
// cur_ is sitting on the ';' that starts the field (which is at column 1).
// ---------------------------------------------------------------------------
Token Tokenizer::read_semicolon_field(int tok_line, int tok_col) {
  ++cur_; ++col_; // consume the opening ';'
  // Skip to end of the opening line (there may be text after ';' on that line,
  // per CIF spec the field text starts on the NEXT line — but many files
  // have content on the same line as ';'. We follow GEMMI/ucif: content
  // starts immediately after the opening ';'.
  const char* start = cur_;
  while (!at_end()) {
    char c = *cur_;
    if (c == '\r' || c == '\n') {
      // Consume the newline
      if (c == '\r') { ++cur_; if (!at_end() && *cur_ == '\n') ++cur_; }
      else { ++cur_; }
      advance_line();
      // Check if this line starts with ';' (the closing delimiter)
      if (!at_end() && *cur_ == ';') {
        std::size_t len = static_cast<std::size_t>(cur_ - start);
        ++cur_; col_ = 2; // consume closing ';'
        return make_token(TOKEN_VALUE, start, len, tok_line, tok_col);
      }
    } else {
      ++cur_; ++col_;
    }
  }
  return make_token(TOKEN_VALUE, start,
                    static_cast<std::size_t>(cur_ - start),
                    tok_line, tok_col);
}

// ---------------------------------------------------------------------------
// Unquoted token (keywords, tags, plain values)
// ---------------------------------------------------------------------------
Token Tokenizer::read_unquoted(int tok_line, int tok_col) {
  const char* start = cur_;
  while (!at_end() && is_ordinary(*cur_)) { ++cur_; ++col_; }
  std::size_t len = static_cast<std::size_t>(cur_ - start);

  // Classify by prefix (case-insensitive)
  if (len > 0 && start[0] == '_') {
    return make_token(TOKEN_TAG, start, len, tok_line, tok_col);
  }
  if (iprefix(start, len, "data_", 5)) {
    return make_token(TOKEN_BLOCK_HEADER, start, len, tok_line, tok_col);
  }
  if (iprefix(start, len, "loop_", 5) && len == 5) {
    return make_token(TOKEN_LOOP, start, len, tok_line, tok_col);
  }
  if (iprefix(start, len, "save_", 5)) {
    if (len == 5)
      return make_token(TOKEN_SAVE_END, start, len, tok_line, tok_col);
    else
      return make_token(TOKEN_SAVE_HEADER, start, len, tok_line, tok_col);
  }
  return make_token(TOKEN_VALUE, start, len, tok_line, tok_col);
}

// ---------------------------------------------------------------------------
// next() — main dispatch
// ---------------------------------------------------------------------------
Token Tokenizer::next() {
  skip_whitespace_and_comments();

  if (at_end()) {
    return make_token(TOKEN_EOF, cur_, 0, line_, col_);
  }

  char     c        = *cur_;
  int      tok_line = line_;
  int      tok_col  = col_;

  // Semicolon text field: only when ';' appears at column 1
  if (c == ';' && tok_col == 1) {
    return read_semicolon_field(tok_line, tok_col);
  }

  // Triple-quoted strings (CIF2)
  if ((c == '"' || c == '\'') && peek2() == c && peek3() == c) {
    return read_triple_quoted(c, tok_line, tok_col);
  }

  // Single/double quoted strings
  if (c == '"' || c == '\'') {
    return read_quoted(c, tok_line, tok_col);
  }

  // Everything else: unquoted token
  return read_unquoted(tok_line, tok_col);
}

} // namespace xcif
