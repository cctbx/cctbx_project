// (jEdit options) :folding=explicit:
//{{{ Package, imports
//package ....;
//}}}
/** Java port of the hy36encode() and hy36decode() functions in the
    hybrid_36.py Python prototype/reference implementation.
    See the Python script for more information.

    This file is unrestricted Open Source (https://github.com/cctbx/cctbx_project).
    Please send corrections and enhancements to cctbx@cci.lbl.gov .

    See also: http://cci.lbl.gov/hybrid_36/
    git version at:
    https://github.com/cctbx/cctbx_project/blob/master/iotbx/pdb/Hy36.java

    Ralf W. Grosse-Kunstleve, Vincent B. Chen, Jeff J. Headd, Sep 2007.
 */
public
class Hy36 {

  //{{{ Constants
  private static String digits_base10 = "0123456789";
  private static String digits_upper  = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  private static String digits_lower  = "0123456789abcdefghijklmnopqrstuvwxyz";

  private static boolean first_call = true;
  private static int[] digits_values_upper = new int[128];
  private static int[] digits_values_lower = new int[128];

  private static String value_out_of_range = "value out of range.";
  private static String invalid_number_literal = "invalid number literal.";
  private static String unsupported_width = "unsupported width.";
  //}}}

  //{{{ private encode/decode_pure functions
  private static
  String
  encode_pure(
    String digits,
    int width,
    int value)
  {
    boolean neg = false;
    if (value < 0) {
      neg = true;
      value = -value;
    }
    String buf = "";
    while (true) {
      int rest = value / digits.length();
      buf += digits.charAt(value - rest * digits.length());
      if (rest == 0) break;
      value = rest;
    }
    if (neg) buf += '-';
    String result = "";
    for(int i=buf.length();i<width;i++) result += " ";
    for(int i=buf.length()-1;i>=0;i--) result += buf.charAt(i);
    return result;
  }

  private static
  int
  decode_pure(
    int[] digits_values,
    int digits_size,
    String s)
  {
    boolean have_minus = false;
    boolean have_non_blank = false;
    int value = 0;
    for(int i=0;i<s.length();i++) {
      char si = s.charAt(i);
      if (si < 0 || si > 127) {
        throw new Error(invalid_number_literal);
      }
      if (si == ' ') {
        if (!have_non_blank) continue;
        value *= digits_size;
      }
      else if (si == '-') {
        if (have_non_blank) {
          throw new Error(invalid_number_literal);
        }
        have_non_blank = true;
        have_minus = true;
        continue;
      }
      else {
        have_non_blank = true;
        int dv = digits_values[si];
        if (dv < 0 || dv >= digits_size) {
          throw new Error(invalid_number_literal);
        }
        value *= digits_size;
        value += dv;
      }
    }
    if (have_minus) value = -value;
    return value;
  }
  //}}}

  //{{{ encode
  /** hybrid-36 encoder: converts integer value to string result

        width: must be 4 (e.g. for residue sequence numbers)
                    or 5 (e.g. for atom serial numbers)

        value: integer value to be converted

        return value: String of size width
   */
  public static
  String
  encode(int width, int value)
  {
    int i = value;
    if (width == 4) {
      if (i >= -999) {
        if (i < 10000) {
          return encode_pure(digits_base10, 4, i);
        }
        i -= 10000;
        if (i < 1213056 /* 26*36**3 */) {
          i += 466560 /* 10*36**3 */;
          return encode_pure(digits_upper, 0, i);
        }
        i -= 1213056;
        if (i < 1213056) {
          i += 466560;
          return encode_pure(digits_lower, 0, i);
        }
      }
    }
    else if (width == 5) {
      if (i >= -9999) {
        if (i < 100000) {
          return encode_pure(digits_base10, 5, i);
        }
        i -= 100000;
        if (i < 43670016 /* 26*36**4 */) {
          i += 16796160 /* 10*36**4 */;
          return encode_pure(digits_upper, 0, i);
        }
        i -= 43670016;
        if (i < 43670016) {
          i += 16796160;
          return encode_pure(digits_lower, 0, i);
        }
      }
    }
    else {
      throw new Error(unsupported_width);
    }
    throw new Error(value_out_of_range);
  }
  //}}}

  //{{{ decode
  /** hybrid-36 decoder: converts string s to integer result

        width: must be 4 (e.g. for residue sequence numbers)
                    or 5 (e.g. for atom serial numbers)

        s: string to be converted

        return value: conversion result
   */
  public static
  int
  decode(int width, String s)
  {
    String ie_range="internal error hy36.decode: integer value out of range.";
    if (first_call) {
      first_call = false;
      for(int i=0;i<128;i++) digits_values_upper[i] = -1;
      for(int i=0;i<128;i++) digits_values_lower[i] = -1;
      for(int i=0;i<36;i++) {
        int di = (int) digits_upper.charAt(i);
        if (di < 0 || di > 127) {
          throw new Error(ie_range);
        }
        digits_values_upper[di] = i;
      }
      for(int i=0;i<36;i++) {
        int di = (int) digits_lower.charAt(i);
        if (di < 0 || di > 127) {
          throw new Error(ie_range);
        }
        digits_values_lower[di] = i;
      }
    }
    if (s.length() == width) {
      int di = (int) s.charAt(0);
      if (di >= 0 && di <= 127) {
        if (digits_values_upper[di] >= 10) {
          int result = decode_pure(digits_values_upper, 36, s);
          /* result - 10*36**(width-1) + 10**width */
          if      (width == 4) result -= 456560;
          else if (width == 5) result -= 16696160;
          else {
            throw new Error(unsupported_width);
          }
          return result;
        }
        else if (digits_values_lower[di] >= 10) {
          int result = decode_pure(digits_values_lower, 36, s);
          /* result + 16*36**(width-1) + 10**width */
          if      (width == 4) result += 756496;
          else if (width == 5) result += 26973856;
          else {
            throw new Error(unsupported_width);
          }
          return result;
        }
        else {
          int result = decode_pure(digits_values_upper, 10, s);
          if (!(width == 4 || width == 5)) {
            throw new Error(unsupported_width);
          }
          return result;
        }
      }
    }
    throw new Error(invalid_number_literal);
  }
  //}}}

  //{{{ check functions
  private static
  void
  check_str(String result, String expected)
  {
    if (!result.equals(expected)) {
      System.out.println("ERROR: \"" + result + "\" != \"" + expected +"\"");
    }
  }

  private static
  void
  check_int(int result, int expected)
  {
    if (result != expected) {
      System.out.println("ERROR: " + result + " != " + expected);
    }
  }
  //}}}

  //{{{ recycle functions
  private static
  void
  recycle4(int value, String encoded)
  {
    String s = encode(4, value);
    check_str(s, encoded);
    int d = decode(4, s);
    check_int(d, value);
  }

  private static
  void
  recycle5(int value, String encoded)
  {
    String s = encode(5, value);
    check_str(s, encoded);
    int d = decode(5, s);
    check_int(d, value);
  }
  //}}}

  //{{{ check exception functions
  private static
  void
  check_encode_exception(int width, int value, String expected_msg)
  {
    String msg = "";
    try { encode(width, value); }
    catch(Error e) { msg = e.toString(); }
    check_str(msg, "java.lang.Error: " + expected_msg);
  }

  private static
  void
  check_decode_exception(int width, String s, String expected_msg)
  {
    String msg = "";
    try { decode(width, s); }
    catch(Error e) { msg = e.toString(); }
    check_str(msg, "java.lang.Error: " + expected_msg);
  }
  //}}}

  //{{{ kernighan_and_ritchie_rand
  private static int random_seed = 13;
  private static
  int
  kernighan_and_ritchie_rand()
  {
    random_seed = random_seed * 1103515245 + 12345;
    int result = (random_seed / 65536) % 32768;
    if (result < 0) result += 32768;
    return result;
  }
  //}}}

  //{{{ main
  /** for running test cases **/
  public static
  void
  main(String[] args)
  {
    check_int(decode(4, "    "), 0);
    check_int(decode(4, "  -0"), 0);
    recycle4(-999, "-999");
    recycle4(-78, " -78");
    recycle4(-6, "  -6");
    recycle4(0, "   0");
    recycle4(9999, "9999");
    recycle4(10000, "A000");
    recycle4(10001, "A001");
    recycle4(10002, "A002");
    recycle4(10003, "A003");
    recycle4(10004, "A004");
    recycle4(10005, "A005");
    recycle4(10006, "A006");
    recycle4(10007, "A007");
    recycle4(10008, "A008");
    recycle4(10009, "A009");
    recycle4(10010, "A00A");
    recycle4(10011, "A00B");
    recycle4(10012, "A00C");
    recycle4(10013, "A00D");
    recycle4(10014, "A00E");
    recycle4(10015, "A00F");
    recycle4(10016, "A00G");
    recycle4(10017, "A00H");
    recycle4(10018, "A00I");
    recycle4(10019, "A00J");
    recycle4(10020, "A00K");
    recycle4(10021, "A00L");
    recycle4(10022, "A00M");
    recycle4(10023, "A00N");
    recycle4(10024, "A00O");
    recycle4(10025, "A00P");
    recycle4(10026, "A00Q");
    recycle4(10027, "A00R");
    recycle4(10028, "A00S");
    recycle4(10029, "A00T");
    recycle4(10030, "A00U");
    recycle4(10031, "A00V");
    recycle4(10032, "A00W");
    recycle4(10033, "A00X");
    recycle4(10034, "A00Y");
    recycle4(10035, "A00Z");
    recycle4(10036, "A010");
    recycle4(10046, "A01A");
    recycle4(10071, "A01Z");
    recycle4(10072, "A020");
    recycle4(10000+36*36-1, "A0ZZ");
    recycle4(10000+36*36, "A100");
    recycle4(10000+36*36*36-1, "AZZZ");
    recycle4(10000+36*36*36, "B000");
    recycle4(10000+26*36*36*36-1, "ZZZZ");
    recycle4(10000+26*36*36*36, "a000");
    recycle4(10000+26*36*36*36+35, "a00z");
    recycle4(10000+26*36*36*36+36, "a010");
    recycle4(10000+26*36*36*36+36*36-1, "a0zz");
    recycle4(10000+26*36*36*36+36*36, "a100");
    recycle4(10000+26*36*36*36+36*36*36-1, "azzz");
    recycle4(10000+26*36*36*36+36*36*36, "b000");
    recycle4(10000+2*26*36*36*36-1, "zzzz");
    //
    check_int(decode(5, "     "), 0);
    check_int(decode(5, "   -0"), 0);
    recycle5(-9999, "-9999");
    recycle5(-123, " -123");
    recycle5(-45, "  -45");
    recycle5(-6, "   -6");
    recycle5(0, "    0");
    recycle5(12, "   12");
    recycle5(345, "  345");
    recycle5(6789, " 6789");
    recycle5(99999, "99999");
    recycle5(100000, "A0000");
    recycle5(100010, "A000A");
    recycle5(100035, "A000Z");
    recycle5(100036, "A0010");
    recycle5(100046, "A001A");
    recycle5(100071, "A001Z");
    recycle5(100072, "A0020");
    recycle5(100000+36*36-1, "A00ZZ");
    recycle5(100000+36*36, "A0100");
    recycle5(100000+36*36*36-1, "A0ZZZ");
    recycle5(100000+36*36*36, "A1000");
    recycle5(100000+36*36*36*36-1, "AZZZZ");
    recycle5(100000+36*36*36*36, "B0000");
    recycle5(100000+2*36*36*36*36, "C0000");
    recycle5(100000+26*36*36*36*36-1, "ZZZZZ");
    recycle5(100000+26*36*36*36*36, "a0000");
    recycle5(100000+26*36*36*36*36+36-1, "a000z");
    recycle5(100000+26*36*36*36*36+36, "a0010");
    recycle5(100000+26*36*36*36*36+36*36-1, "a00zz");
    recycle5(100000+26*36*36*36*36+36*36, "a0100");
    recycle5(100000+26*36*36*36*36+36*36*36-1, "a0zzz");
    recycle5(100000+26*36*36*36*36+36*36*36, "a1000");
    recycle5(100000+26*36*36*36*36+36*36*36*36-1, "azzzz");
    recycle5(100000+26*36*36*36*36+36*36*36*36, "b0000");
    recycle5(100000+2*26*36*36*36*36-1, "zzzzz");
    //
    check_encode_exception(4, -1000, "value out of range.");
    check_encode_exception(4, 2436112, "value out of range.");
    check_encode_exception(5, -10000, "value out of range.");
    check_encode_exception(5, 87440032, "value out of range.");
    //
    check_decode_exception(4, "", "invalid number literal.");
    check_decode_exception(4, "    0", "invalid number literal.");
    check_decode_exception(4, " abc", "invalid number literal.");
    check_decode_exception(4, "abc-", "invalid number literal.");
    check_decode_exception(4, "A=BC", "invalid number literal.");
    check_decode_exception(4, "40a0", "invalid number literal.");
    check_decode_exception(4, "40A0", "invalid number literal.");
    check_decode_exception(5, "", "invalid number literal.");
    check_decode_exception(5, "     0", "invalid number literal.");
    check_decode_exception(5, " abcd", "invalid number literal.");
    check_decode_exception(5, "ABCD-", "invalid number literal.");
    check_decode_exception(5, "a=bcd", "invalid number literal.");
    check_decode_exception(5, "410b0", "invalid number literal.");
    check_decode_exception(5, "410B0", "invalid number literal.");
    //
    check_encode_exception(3, 0, "unsupported width.");
    check_encode_exception(6, 0, "unsupported width.");
    check_decode_exception(3, "AAA", "unsupported width.");
    check_decode_exception(6, "zzzzzz", "unsupported width.");
    //
    int value = -9999;
    while (value < 100000+2*26*36*36*36*36) {
      check_int(decode(5, encode(5, value)), value);
      value += kernighan_and_ritchie_rand() % 10000;
    }
    //
    System.out.println("OK");
  }
  //}}}
}
