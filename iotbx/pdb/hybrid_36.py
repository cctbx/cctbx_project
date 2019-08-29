"""
               Prototype/reference implementation for
                      encoding and decoding
                     atom serial numbers and
                     residue sequence numbers
                          in PDB files.

PDB ATOM and HETATM records reserve columns 7-11 for the atom serial
number. This 5-column number is used as a reference in the CONECT
records, which also reserve exactly five columns for each serial
number.

With the decimal counting system only up to 99999 atoms can be stored
and uniquely referenced in a PDB file. A simple extension to enable
processing of more atoms is to adopt a counting system with more than
ten digits. To maximize backward compatibility, the counting system is
only applied for numbers greater than 99999. The "hybrid-36" counting
system implemented in this file is:

  ATOM      1
  ...
  ATOM  99999
  ATOM  A0000
  ATOM  A0001
  ...
  ATOM  A0009
  ATOM  A000A
  ...
  ATOM  A000Z
  ATOM  ZZZZZ
  ATOM  a0000
  ...
  ATOM  zzzzz

I.e. the first 99999 serial numbers are represented as usual. The
following atoms use a base-36 system (10 digits + 26 letters) with
upper-case letters. 43670016 (26*36**4) additional atoms can be
numbered this way. If there are more than 43770015 (99999+43670016)
atoms, a base-36 system with lower-case letters is used, allowing for
43670016 additional atoms. I.e. in total 87440031 (99999+2*43670016)
atoms can be stored and uniquely referenced via CONECT records.

The counting system is designed to avoid lower-case letters until the
range of numbers addressable by upper-case letters is exhausted.
Importantly, with this counting system the distinction between
"traditional" and "extended" PDB files becomes evident only if there
are more than 99999 atoms to be stored. Programs that are
updated to support the hybrid-36 counting system will continue to
interoperate with programs that do not as long as there are less than
100000 atoms.

PDB ATOM and HETATM records also reserve columns 23-26 for the residue
sequence number. This 4-column number is used as a reference in other
record types (SSBOND, LINK, HYDBND, SLTBRG, CISPEP), which also reserve
exactly four columns for each sequence number.

With the decimal counting system only up to 9999 residues per chain can
be stored and uniquely referenced in a PDB file. If the hybrid-36
system is adopted, 1213056 (26*36**3) additional residues can be
numbered using upper-case letters, and the same number again using
lower-case letters. I.e. in total each chain may contain up to 2436111
(9999+2*1213056) residues that can be uniquely referenced from the
other record types given above.

The implementation in this file should run with Python 2.6 or higher.
There are no other requirements. Run this script without arguments to
obtain usage examples.

Note that there are only about 60 lines of "real" code. The rest is
documentation and unit tests.

To update an existing program to support the hybrid-36 counting system,
simply replace the existing read/write source code for integer values
with equivalents of the hy36decode() and hy36encode() functions below.

This file is unrestricted Open Source (cctbx.sf.net).
Please send corrections and enhancements to cctbx@cci.lbl.gov .

See also:
  http://cci.lbl.gov/hybrid_36/
  http://www.pdb.org/ "Dictionary & File Formats"

Ralf W. Grosse-Kunstleve, Feb 2007.
"""
from __future__ import absolute_import, division, print_function
try:
  from six.moves import range
  from six.moves import zip
except ImportError:
  pass

digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
digits_lower = digits_upper.lower()
digits_upper_values = dict([pair for pair in zip(digits_upper, range(36))])
digits_lower_values = dict([pair for pair in zip(digits_lower, range(36))])

def encode_pure(digits, value):
  "encodes value using the given digits"
  assert value >= 0
  if (value == 0): return digits[0]
  n = len(digits)
  result = []
  while (value != 0):
    rest = value // n
    result.append(digits[value - rest * n])
    value = rest
  result.reverse()
  return "".join(result)

def decode_pure(digits_values, s):
  "decodes the string s using the digit, value associations for each character"
  result = 0
  n = len(digits_values)
  for c in s:
    result *= n
    result += digits_values[c]
  return result

def hy36encode(width, value):
  "encodes value as base-10/upper-case base-36/lower-case base-36 hybrid"
  i = value
  if (i >= 1-10**(width-1)):
    if (i < 10**width):
      return ("%%%dd" % width) % i
    i -= 10**width
    if (i < 26*36**(width-1)):
      i += 10*36**(width-1)
      return encode_pure(digits_upper, i)
    i -= 26*36**(width-1)
    if (i < 26*36**(width-1)):
      i += 10*36**(width-1)
      return encode_pure(digits_lower, i)
  raise ValueError("value out of range.")

def hy36decode(width, s):
  "decodes base-10/upper-case base-36/lower-case base-36 hybrid"
  if (len(s) == width):
    f = s[0]
    if (f == "-" or f == " " or f.isdigit()):
      try: return int(s)
      except ValueError: pass
      if (s == " "*width): return 0
    elif (f in digits_upper_values):
      try: return decode_pure(
        digits_values=digits_upper_values, s=s) - 10*36**(width-1) + 10**width
      except KeyError: pass
    elif (f in digits_lower_values):
      try: return decode_pure(
        digits_values=digits_lower_values, s=s) + 16*36**(width-1) + 10**width
      except KeyError: pass
  raise ValueError("invalid number literal.")

def exercise(hy36enc=hy36encode, hy36dec=hy36decode):
  for digits,digits_values in [(digits_upper, digits_upper_values),
                               (digits_lower, digits_lower_values)]:
    for value in range(1000):
      s = encode_pure(digits=digits_upper, value=value)
      d = decode_pure(digits_values=digits_upper_values, s=s)
      assert d == value
  #
  def recycle4(value, encoded):
    s = hy36enc(width=4, value=value)
    assert s == encoded
    d = hy36dec(width=4, s=s)
    assert d == value
  #
  assert hy36dec(width=4, s="    ") == 0
  assert hy36dec(width=4, s="  -0") == 0
  recycle4(-999, "-999")
  recycle4(-78, " -78")
  recycle4(-6, "  -6")
  recycle4(0, "   0")
  recycle4(9999, "9999")
  recycle4(10000, "A000")
  recycle4(10001, "A001")
  recycle4(10002, "A002")
  recycle4(10003, "A003")
  recycle4(10004, "A004")
  recycle4(10005, "A005")
  recycle4(10006, "A006")
  recycle4(10007, "A007")
  recycle4(10008, "A008")
  recycle4(10009, "A009")
  recycle4(10010, "A00A")
  recycle4(10011, "A00B")
  recycle4(10012, "A00C")
  recycle4(10013, "A00D")
  recycle4(10014, "A00E")
  recycle4(10015, "A00F")
  recycle4(10016, "A00G")
  recycle4(10017, "A00H")
  recycle4(10018, "A00I")
  recycle4(10019, "A00J")
  recycle4(10020, "A00K")
  recycle4(10021, "A00L")
  recycle4(10022, "A00M")
  recycle4(10023, "A00N")
  recycle4(10024, "A00O")
  recycle4(10025, "A00P")
  recycle4(10026, "A00Q")
  recycle4(10027, "A00R")
  recycle4(10028, "A00S")
  recycle4(10029, "A00T")
  recycle4(10030, "A00U")
  recycle4(10031, "A00V")
  recycle4(10032, "A00W")
  recycle4(10033, "A00X")
  recycle4(10034, "A00Y")
  recycle4(10035, "A00Z")
  recycle4(10036, "A010")
  recycle4(10046, "A01A")
  recycle4(10071, "A01Z")
  recycle4(10072, "A020")
  recycle4(10000+36**2-1, "A0ZZ")
  recycle4(10000+36**2, "A100")
  recycle4(10000+36**3-1, "AZZZ")
  recycle4(10000+36**3, "B000")
  recycle4(10000+26*36**3-1, "ZZZZ")
  recycle4(10000+26*36**3, "a000")
  recycle4(10000+26*36**3+35, "a00z")
  recycle4(10000+26*36**3+36, "a010")
  recycle4(10000+26*36**3+36**2-1, "a0zz")
  recycle4(10000+26*36**3+36**2, "a100")
  recycle4(10000+26*36**3+36**3-1, "azzz")
  recycle4(10000+26*36**3+36**3, "b000")
  recycle4(10000+2*26*36**3-1, "zzzz")
  #
  def recycle5(value, encoded):
    s = hy36enc(width=5, value=value)
    assert s == encoded
    d = hy36dec(width=5, s=s)
    assert d == value
  #
  assert hy36dec(width=5, s="     ") == 0
  assert hy36dec(width=5, s="   -0") == 0
  recycle5(-9999, "-9999")
  recycle5(-123, " -123")
  recycle5(-45, "  -45")
  recycle5(-6, "   -6")
  recycle5(0, "    0")
  recycle5(12, "   12")
  recycle5(345, "  345")
  recycle5(6789, " 6789")
  recycle5(99999, "99999")
  recycle5(100000, "A0000")
  recycle5(100010, "A000A")
  recycle5(100035, "A000Z")
  recycle5(100036, "A0010")
  recycle5(100046, "A001A")
  recycle5(100071, "A001Z")
  recycle5(100072, "A0020")
  recycle5(100000+36**2-1, "A00ZZ")
  recycle5(100000+36**2, "A0100")
  recycle5(100000+36**3-1, "A0ZZZ")
  recycle5(100000+36**3, "A1000")
  recycle5(100000+36**4-1, "AZZZZ")
  recycle5(100000+36**4, "B0000")
  recycle5(100000+2*36**4, "C0000")
  recycle5(100000+26*36**4-1, "ZZZZZ")
  recycle5(100000+26*36**4, "a0000")
  recycle5(100000+26*36**4+36-1, "a000z")
  recycle5(100000+26*36**4+36, "a0010")
  recycle5(100000+26*36**4+36**2-1, "a00zz")
  recycle5(100000+26*36**4+36**2, "a0100")
  recycle5(100000+26*36**4+36**3-1, "a0zzz")
  recycle5(100000+26*36**4+36**3, "a1000")
  recycle5(100000+26*36**4+36**4-1, "azzzz")
  recycle5(100000+26*36**4+36**4, "b0000")
  recycle5(100000+2*26*36**4-1, "zzzzz")
  #
  for width in [4,5]:
    for value in [-(10**(width-1)), 10**width+2*26*36**(width-1)]:
      try: hy36enc(width=width, value=value)
      except (ValueError, RuntimeError) as e:
        assert str(e) == "value out of range."
      else: raise RuntimeError("Exception expected.")
  #
  for width,ss in [
      (4, ["", "    0", " abc", "abc-", "A=BC", "40a0", "40A0"]),
      (5, ["", "     0", " abcd", "ABCD-", "a=bcd", "410b0", "410B0"])]:
    for s in ss:
      try: hy36dec(width, s=s)
      except (ValueError, RuntimeError) as e:
        assert str(e) == "invalid number literal."
      else: raise RuntimeError("Exception expected.")
  #
  import random
  value = -9999
  while value < 100000+2*26*36**4:
    try:
      s = hy36enc(width=5, value=value)
      d = hy36dec(width=5, s=s)
    except Exception:
      print("value:", value)
      raise
    assert d == value
    value += random.randint(0, 10000)

def run():
  import sys
  def usage():
    c = sys.argv[0]
    print("usage examples:")
    print("  python %s info" % c)
    print("  python %s exercise" % c)
    print("  python %s hy36encode 4 10000" % c)
    print("  python %s hy36decode 4 zzzz" % c)
    print("  python %s hy36encode 5 100000" % c)
    print("  python %s hy36decode 5 zzzzz" % c)
    sys.exit(1)
  if (len(sys.argv) < 2):
    usage()
  task = sys.argv[1]
  if (task == "info"):
    sys.stdout.write(__doc__)
    return
  if (task == "exercise"):
    exercise()
    print("OK")
  else:
    if (   len(sys.argv) != 4
        or task not in ["hy36encode", "hy36decode"]):
      usage()
    f = globals()[task]
    width = int(sys.argv[2])
    assert width > 0
    s = sys.argv[3]
    if (task == "hy36encode"):
      print(f(value=int(s), width=width))
    else:
      s = " "*max(0, width-len(s)) + s
      print(f(s=s, width=width))

if (__name__ == "__main__"):
  run()
