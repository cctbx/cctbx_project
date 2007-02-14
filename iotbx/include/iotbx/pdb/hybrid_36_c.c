/*! C port of the hy36encode() and hy36decode() functions in the
    hybrid_36.py Python prototype/reference implementation.
    See the Python script for more information.

    This file has no dependencies other than ANSI C headers.
    Optionally, use hybrid_36_c.h, or simply copy the declarations
    into your code.

    This file is unrestricted Open Source (cctbx.sf.net).
    Please send corrections and enhancements to cctbx@cci.lbl.gov .

    Ralf W. Grosse-Kunstleve, Feb 2007.
 */

/* The following #include may be commented out.
   It is here only to enforce consistency of the declarations
   and the definitions.
 */
#include <iotbx/pdb/hybrid_36_c.h>

#include <stdio.h>
#include <ctype.h>

static const char* digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char* digits_lower = "0123456789abcdefghijklmnopqrstuvwxyz";

static const char* value_out_of_range = "value out of range.";
static const char* invalid_number_literal = "invalid number literal.";
static const char* unsupported_width = "unsupported width.";

static
const char*
encode_pure(
  const char* digits,
  unsigned digits_size,
  unsigned value,
  char* result)
{
  if (value < 0) return value_out_of_range;
  if (value == 0) {
    result[0] = digits[0];
    result[1] = '\0';
    return NULL;
  }
  char buf[16];
  unsigned i = 0;
  while (value != 0) {
    unsigned rest = value / digits_size;
    buf[i++] = digits[value - rest * digits_size];
    value = rest;
  }
  while (i != 0) {
    *result++ = buf[--i];
  }
  *result = '\0';
  return NULL;
}

static
const char*
decode_pure(
  const int* digits_values,
  unsigned digits_size,
  const char* s,
  unsigned s_size,
  int* result)
{
  unsigned value = 0;
  unsigned i = 0;
  for(;i<s_size;i++) {
    value *= digits_size;
    int si = s[i];
    if (si < 0 || si > 127) return invalid_number_literal;
    int dv = digits_values[si];
    if (dv < 0) return invalid_number_literal;
    value += (unsigned) dv;
  }
  *result = value;
  return NULL;
}

const char*
hy36encode(unsigned width, int value, char* result)
{
  int i = value;
  if (width == 4U) {
    if (i >= -999) {
      if (i < 10000) {
        sprintf(result, "%4d", i);
        return NULL;
      }
      i -= 10000;
      if (i < 1213056 /* 26*36**3 */) {
        i += 466560 /* 10*36**3 */;
        return encode_pure(digits_upper, 36U, (unsigned) i, result);
      }
      i -= 1213056;
      if (i < 1213056) {
        i += 466560;
        return encode_pure(digits_lower, 36U, (unsigned) i, result);
      }
    }
  }
  else if (width == 5U) {
    if (i >= -9999) {
      if (i < 100000) {
        sprintf(result, "%5d", i);
        return NULL;
      }
      i -= 100000;
      if (i < 43670016 /* 26*36**4 */) {
        i += 16796160 /* 10*36**4 */;
        return encode_pure(digits_upper, 36U, (unsigned) i, result);
      }
      i -= 43670016;
      if (i < 43670016) {
        i += 16796160;
        return encode_pure(digits_lower, 36U, (unsigned) i, result);
      }
    }
  }
  else {
    return unsupported_width;
  }
  return value_out_of_range;
}

const char*
hy36decode(unsigned width, const char* s, unsigned s_size, int* result)
{
  static int first_call = 1;
  static int digits_values_upper[128U];
  static int digits_values_lower[128U];
  static const char*
    ie_range = "internal error: int value of base-36 digits out of range.";
  unsigned i;
  int di, f, neg;
  const char* diag;
  if (first_call) {
    first_call = 0;
    for(i=0;i<128U;i++) digits_values_upper[i] = -1;
    for(i=0;i<128U;i++) digits_values_lower[i] = -1;
    for(i=0;i<36U;i++) {
      di = digits_upper[i];
      if (di < 0 || di > 127) return ie_range;
      digits_values_upper[di] = i;
    }
    for(i=0;i<36U;i++) {
      di = digits_lower[i];
      if (di < 0 || di > 127) return ie_range;
      digits_values_lower[di] = i;
    }
  }
  if (s_size == width) {
    f = s[0];
    if (f == '-' || f == ' ' || isdigit(f)) {
      i = 0;
      while (i < s_size && s[i] == ' ') i++;
      if (i == s_size) {
        *result = 0;
        return NULL;
      }
      neg = 0;
      if (s[i] == '-') {
        neg = 1;
        i++;
        if (i == s_size) return invalid_number_literal;
      }
      diag = decode_pure(digits_values_upper, 10U, s+i, s_size-i, result);
      if (diag) return diag;
      if (neg) *result = -(*result);
      return NULL;
    }
    else if (isupper(f)) {
      diag = decode_pure(digits_values_upper, 36U, s, s_size, result);
      if (diag == 0) {
        /* result - 10*36**(width-1) + 10**width */
        if      (width == 4U) (*result) -= 456560;
        else if (width == 5U) (*result) -= 16696160;
        else return unsupported_width;
        return NULL;
      }
    }
    else if (islower(f)) {
      diag = decode_pure(digits_values_lower, 36U, s, s_size, result);
      if (diag == 0) {
        /* result + 16*36**(width-1) + 10**width */
        if      (width == 4U) (*result) += 756496;
        else if (width == 5U) (*result) += 26973856;
        else return unsupported_width;
        return NULL;
      }
    }
  }
  return invalid_number_literal;
}
