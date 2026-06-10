// cctbx_project/xcif/include/xcif/numeric.h
//
// Free functions for CIF numeric value interpretation.
// Predicates are inlined; heavy parsers are declared here,
// implemented in numeric.cpp (Step 8).
#ifndef XCIF_NUMERIC_H
#define XCIF_NUMERIC_H

#include <xcif/string_view.h>
#include <utility>

namespace xcif {

// ── Predicates (trivial, always inlined) ───────────────────────────

inline bool is_unknown(const string_view& sv) {
  return sv.size() == 1 && sv[0] == '.';
}

inline bool is_inapplicable(const string_view& sv) {
  return sv.size() == 1 && sv[0] == '?';
}

inline bool is_null(const string_view& sv) {
  return is_unknown(sv) || is_inapplicable(sv);
}

// ── Numeric conversions (implemented in numeric.cpp) ───────────────

// Parse a CIF value as a double.
// - Strips standard uncertainty: "1.542(3)" → 1.542
// - Returns NaN for "." and "?"
// - Returns NaN for non-numeric strings
// - Handles optional sign, decimal point, and scientific notation
double as_double(const string_view& sv);

// Parse a CIF value as an integer.
// - Handles optional sign: [+-]?[0-9]+
// - Throws std::invalid_argument on non-integer input
int as_int(const string_view& sv);

// Parse a CIF value as a double with standard uncertainty.
// - Returns (value, su). If no SU present, su = 0.0
// - SU is scaled by decimal position: "1.542(3)" → (1.542, 0.003)
// - Returns (NaN, 0.0) for "." and "?"
std::pair<double, double> as_double_with_su(const string_view& sv);

} // namespace xcif
#endif // XCIF_NUMERIC_H
