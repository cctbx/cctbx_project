// cctbx_project/xcif/include/xcif/tokenizer.h
#ifndef XCIF_TOKENIZER_H
#define XCIF_TOKENIZER_H

#include <string>
#include <cstddef>

namespace xcif {

// ---------------------------------------------------------------------------
// Token types
// ---------------------------------------------------------------------------
enum TokenType {
  TOKEN_EOF          = 0,
  TOKEN_TAG          = 1,  // _category.item
  TOKEN_VALUE        = 2,  // any value: unquoted, quoted, semicolon-field
  TOKEN_BLOCK_HEADER = 3,  // data_NAME
  TOKEN_LOOP         = 4,  // loop_
  TOKEN_SAVE_HEADER  = 5,  // save_NAME
  TOKEN_SAVE_END     = 6   // save_  (alone)
};

// ---------------------------------------------------------------------------
// Token — zero-copy view into the source buffer
// ---------------------------------------------------------------------------
struct Token {
  TokenType   type;
  const char* ptr;   // points into original source buffer (not NUL-terminated)
  std::size_t len;
  int         line;  // 1-based
  int         col;   // 1-based

  // Convenience: materialise as std::string (copies; use only when needed)
  std::string as_str() const { return std::string(ptr, len); }
};

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------
class Tokenizer {
 public:
  // Construct from an in-memory buffer (zero-copy; buffer must outlive tokenizer)
  Tokenizer(const char* data, std::size_t length,
            const char* source_name = "<input>");

  // Convenience: NUL-terminated string (length = strlen)
  Tokenizer(const char* data, const char* source_name = "<input>");

  // Return the next token.  Returns TOKEN_EOF indefinitely at end of input.
  Token next();

  const char* source_name() const { return source_name_; }

 private:
  const char* start_;
  const char* cur_;
  const char* end_;
  const char* source_name_;
  int         line_;
  int         col_;

  // Internal helpers
  void        skip_whitespace_and_comments();
  Token       make_token(TokenType type, const char* ptr,
                         std::size_t len, int line, int col) const;
  Token       read_quoted(char delim, int line, int col);
  Token       read_triple_quoted(char delim, int line, int col);
  Token       read_semicolon_field(int line, int col);
  Token       read_unquoted(int line, int col);
  void        advance_line();
  bool        at_end() const { return cur_ >= end_; }
  char        peek()   const { return at_end() ? '\0' : *cur_; }
  char        peek2()  const { return (cur_ + 1 < end_) ? cur_[1] : '\0'; }
  char        peek3()  const { return (cur_ + 2 < end_) ? cur_[2] : '\0'; }
};

} // namespace xcif

#endif // XCIF_TOKENIZER_H
