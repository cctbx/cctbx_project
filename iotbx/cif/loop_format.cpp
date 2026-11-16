#include <iotbx/cif/loop_format.h>
#include <stdexcept>
#include <sstream>

namespace iotbx { namespace cif {

namespace {

  bool
  is_continuation(unsigned char c)
  {
    return (c & 0xC0) == 0x80;
  }

  std::size_t
  code_points(std::string const& s)
  {
    std::size_t n = 0;
    for (std::size_t i = 0; i < s.size(); i++) {
      if (!is_continuation(s[i])) n++;
    }
    return n;
  }

  // what str.isspace() reports: Unicode White_Space plus U+001C..U+001F
  bool
  python_is_space(unsigned long cp)
  {
    if (cp < 0x80) {
      return (cp >= 0x09 && cp <= 0x0D) || (cp >= 0x1C && cp <= 0x20);
    }
    return cp == 0x85 || cp == 0xA0 || cp == 0x1680
        || (cp >= 0x2000 && cp <= 0x200A)
        || cp == 0x2028 || cp == 0x2029 || cp == 0x202F
        || cp == 0x205F || cp == 0x3000;
  }

  // bytes in the UTF-8 sequence a lead byte starts; 0 for no lead byte
  std::size_t
  sequence_length(unsigned char b0)
  {
    if (b0 < 0x80) return 1;
    if ((b0 & 0xE0) == 0xC0) return 2;
    if ((b0 & 0xF0) == 0xE0) return 3;
    if ((b0 & 0xF8) == 0xF0) return 4;
    return 0;
  }

  // the code point of s[i..i+n); false unless that is one UTF-8 sequence
  bool
  decode(std::string const& s, std::size_t i, std::size_t n, unsigned long& cp)
  {
    static const unsigned char lead_mask[5] = {0, 0x7F, 0x1F, 0x0F, 0x07};
    unsigned char b0 = s[i];
    if (n == 0 || n != sequence_length(b0) || i + n > s.size()) return false;
    cp = b0 & lead_mask[n];
    for (std::size_t k = 1; k < n; k++) {
      unsigned char b = s[i + k];
      if (!is_continuation(b)) return false;
      cp = (cp << 6) | (b & 0x3F);
    }
    return true;
  }

  // the size of s after str.rstrip(), not cutting before begin
  std::size_t
  rstripped_size(std::string const& s, std::size_t begin)
  {
    std::size_t end = s.size();
    while (end > begin) {
      std::size_t lead = end - 1;
      while (lead > begin && is_continuation(s[lead]) && end - lead < 4) {
        lead--;
      }
      unsigned long cp;
      if (!decode(s, lead, end - lead, cp) || !python_is_space(cp)) break;
      end = lead;
    }
    return end;
  }

  // byte index past the leading whitespace of s, as str.lstrip() counts it
  std::size_t
  leading_spaces(std::string const& s)
  {
    std::size_t i = 0;
    while (i < s.size()) {
      unsigned long cp;
      std::size_t n = sequence_length(s[i]);
      if (!decode(s, i, n, cp) || !python_is_space(cp)) break;
      i += n;
    }
    return i;
  }

  bool
  is_quote(char c)
  {
    return c == '\'' || c == '"';
  }

  // the rules of format_value, in the order it applies them; any newline
  // too, as loop.show keeps text fields out of the column width
  bool
  value_needs_format_value(std::string const& v)
  {
    static const std::string whitespace(" \t\n\r\x0b\x0c");
    if (v.empty() || v.find('\n') != std::string::npos) return true;
    if (is_quote(v[0])) {
      // already quoted: a second quote after printable characters only
      for (std::size_t j = 1; j < v.size(); j++) {
        unsigned char c = v[j];
        if (is_quote(c)) return false;
        if (c != '\t' && (c < 0x20 || c > 0x7E)) break;
      }
      return v.find_first_of(whitespace) != std::string::npos;
    }
    std::size_t i = leading_spaces(v);
    if (i < v.size() && v[i] == ';'
        && v.find(';', i + 1) != std::string::npos) {
      return true;   // a semicolon text field
    }
    if (v[0] == '#' || v[0] == '$' || v[0] == '[' || v[0] == ']'
        || v[0] == '_') {
      return true;
    }
    return v.find_first_of(whitespace) != std::string::npos;
  }

  bool
  is_ascii_digit(char c)
  {
    return c >= '0' && c <= '9';
  }

  // the grammar of first_non_plain_number in loop_format.h
  bool
  is_plain_number(std::string const& s)
  {
    std::size_t n = s.size();
    if (n == 0 || n > 40) return false;
    std::size_t i = 0;
    if (s[i] == '+' || s[i] == '-') i++;
    std::size_t digits = 0;
    while (i < n && is_ascii_digit(s[i])) { i++; digits++; }
    if (i < n && s[i] == '.') {
      i++;
      while (i < n && is_ascii_digit(s[i])) { i++; digits++; }
    }
    if (digits == 0) return false;
    if (i < n && (s[i] == 'e' || s[i] == 'E')) {
      i++;
      if (i < n && (s[i] == '+' || s[i] == '-')) i++;
      std::size_t exponent_digits = 0;
      while (i < n && is_ascii_digit(s[i])) { i++; exponent_digits++; }
      if (exponent_digits == 0 || exponent_digits > 2) return false;
    }
    return i == n;
  }

} // namespace <anonymous>

  std::string
  format_loop_rows(
    af::const_ref<af::shared<std::string> > const& columns,
    std::size_t start,
    std::size_t stop,
    std::string const& prefix,
    af::const_ref<int> const& widths)
  {
    std::size_t n_columns = columns.size();
    std::size_t n_rows = (n_columns == 0 ? 0 : columns[0].size());
    for (std::size_t j = 1; j < n_columns; j++) {
      if (columns[j].size() != n_rows) {
        std::ostringstream msg;
        msg << "format_loop_rows: column " << j << " has "
            << columns[j].size() << " rows, column 0 has " << n_rows;
        throw std::invalid_argument(msg.str());
      }
    }
    if (start > stop || stop > n_rows) {
      std::ostringstream msg;
      msg << "format_loop_rows: row range " << start << ".." << stop
          << " outside 0.." << n_rows;
      throw std::invalid_argument(msg.str());
    }
    bool aligned = (widths.size() != 0);
    if (aligned && widths.size() != n_columns) {
      std::ostringstream msg;
      msg << "format_loop_rows: " << widths.size() << " widths for "
          << n_columns << " columns";
      throw std::invalid_argument(msg.str());
    }
    // one reservation for the whole range
    std::size_t separators
      = (n_columns == 0 ? 0 : n_columns - 1) * (aligned ? 2 : 1);
    std::size_t per_row = prefix.size() + separators + 1;
    std::size_t size = (stop - start) * per_row;
    for (std::size_t j = 0; j < n_columns; j++) {
      std::size_t w = (aligned ? static_cast<std::size_t>(
        widths[j] < 0 ? -widths[j] : widths[j]) : 0);
      for (std::size_t i = start; i < stop; i++) {
        std::size_t len = columns[j][i].size();
        size += (len > w ? len : w);
      }
    }
    std::string result;
    result.reserve(size);
    for (std::size_t i = start; i < stop; i++) {
      std::size_t line_begin = result.size();
      result += prefix;
      for (std::size_t j = 0; j < n_columns; j++) {
        std::string const& value = columns[j][i];
        if (!aligned) {
          if (j != 0) result += ' ';
          result += value;
          continue;
        }
        if (j != 0) result += "  ";
        int width = widths[j];
        if (width == 0) {
          result += value;
          continue;
        }
        std::size_t field = static_cast<std::size_t>(
          width < 0 ? -width : width);
        std::size_t n = code_points(value);
        std::size_t pad = (field > n ? field - n : 0);
        if (width > 0) result.append(pad, ' ');
        result += value;
        if (width < 0) result.append(pad, ' ');
      }
      if (aligned) result.resize(rstripped_size(result, line_begin));
      result += '\n';
    }
    return result;
  }

  bool
  column_needs_format_value(af::const_ref<std::string> const& column)
  {
    for (std::size_t i = 0; i < column.size(); i++) {
      if (value_needs_format_value(column[i])) return true;
    }
    return false;
  }

  std::size_t
  first_non_plain_number(
    af::const_ref<std::string> const& column,
    std::size_t start)
  {
    for (std::size_t i = start; i < column.size(); i++) {
      std::string const& value = column[i];
      if (value.size() == 1 && (value[0] == '.' || value[0] == '?')) {
        continue;
      }
      if (!is_plain_number(value)) return i;
    }
    return column.size();
  }

}} // namespace iotbx::cif
