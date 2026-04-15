// cctbx_project/xcif/numeric.cpp
#include <xcif/numeric.h>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace xcif {

static const double XCIF_NAN = std::numeric_limits<double>::quiet_NaN();

// ── as_double ──────────────────────────────────────────────────────

double as_double(const string_view& sv) {
  if (sv.empty() || is_null(sv)) return XCIF_NAN;

  const char* p = sv.data();
  std::size_t len = sv.size();

  // Strip SU if present: "1.542(3)" → parse only "1.542"
  const char* paren = static_cast<const char*>(
    std::memchr(p, '(', len));
  std::size_t parse_len = paren
    ? static_cast<std::size_t>(paren - p)
    : len;

  // Stack buffer + NUL-terminate for strtod
  char buf[32];
  if (parse_len >= sizeof(buf)) return XCIF_NAN;
  std::memcpy(buf, p, parse_len);
  buf[parse_len] = '\0';

  char* end;
  double val = std::strtod(buf, &end);
  if (end == buf) return XCIF_NAN;
  return val;
}

// ── as_int ─────────────────────────────────────────────────────────

int as_int(const string_view& sv) {
  if (sv.empty())
    throw std::invalid_argument("as_int: empty string");

  const char* p = sv.data();
  std::size_t len = sv.size();
  std::size_t i = 0;
  bool neg = false;

  if (p[0] == '-')      { neg = true; ++i; }
  else if (p[0] == '+') { ++i; }

  if (i == len)
    throw std::invalid_argument("as_int: no digits");

  int result = 0;
  for (; i < len; ++i) {
    char c = p[i];
    if (c < '0' || c > '9')
      throw std::invalid_argument("as_int: non-digit character");
    result = result * 10 + (c - '0');
  }
  return neg ? -result : result;
}

// ── as_double_with_su ──────────────────────────────────────────────

std::pair<double, double> as_double_with_su(const string_view& sv) {
  if (sv.empty() || is_null(sv))
    return std::make_pair(XCIF_NAN, 0.0);

  const char* p = sv.data();
  std::size_t len = sv.size();

  // Fast path: no SU
  const char* paren = static_cast<const char*>(
    std::memchr(p, '(', len));
  if (!paren)
    return std::make_pair(as_double(sv), 0.0);

  // Parse value portion (before '(')
  std::size_t val_len = static_cast<std::size_t>(paren - p);
  char buf[32];
  if (val_len >= sizeof(buf))
    return std::make_pair(XCIF_NAN, 0.0);
  std::memcpy(buf, p, val_len);
  buf[val_len] = '\0';

  char* end;
  double value = std::strtod(buf, &end);
  if (end == buf)
    return std::make_pair(XCIF_NAN, 0.0);

  // Count decimal places (stop at 'e'/'E')
  int decimals = 0;
  const char* dot = static_cast<const char*>(
    std::memchr(p, '.', val_len));
  if (dot) {
    for (const char* q = dot + 1; q < p + val_len; ++q) {
      if (*q == 'e' || *q == 'E') break;
      ++decimals;
    }
  }

  // Extract SU digits between '(' and ')'
  const char* su_start = paren + 1;
  std::size_t remaining = len - static_cast<std::size_t>(su_start - p);
  const char* close = static_cast<const char*>(
    std::memchr(su_start, ')', remaining));
  if (!close)
    return std::make_pair(value, 0.0);

  int su_int = 0;
  for (const char* q = su_start; q < close; ++q) {
    if (*q >= '0' && *q <= '9')
      su_int = su_int * 10 + (*q - '0');
  }

  // Scale SU: "1.542(3)" → 3 decimals → su = 3 / 10^3 = 0.003
  double su = static_cast<double>(su_int);
  for (int i = 0; i < decimals; ++i)
    su /= 10.0;

  return std::make_pair(value, su);
}

} // namespace xcif
