import fable
from libtbx.test_utils import Exception_expected, show_diff

def try_code_none(f):
  try: f(code=None)
  except TypeError, e: pass
  else: raise Exception_expected

def exercise_unsigned_integer_scan(f):
  try_code_none(f=f)
  assert f(code="") == -1
  assert f(code="0") == 1
  assert f(code="12") == 2
  assert f(code="1a") == 1
  assert f(code="a1") == -1
  assert f(code="a1", start=1) == 2
  assert f(code="a12", start=1) == 3
  assert f(code="a1b", start=1) == 2
  assert f(code="0123456789*2", start=1) == 10
  assert f(code="3456", start=1, stop=-1) == 4
  assert f(code="3456", start=1, stop=3) == 3

def exercise_floating_point_scan_after_exponent_char(f):
  try_code_none(f=f)
  assert f(code="") == -1
  assert f(code="e") == -1
  assert f(code="+") == -1
  assert f(code="-") == -1
  assert f(code="0") == 1
  assert f(code="123456789") == 9
  assert f(code="+0") == 2
  assert f(code="-12") == 3
  assert f(code="-+0") == -1
  assert f(code="-+0", start=1) == 3
  assert f(code="-+34", start=1, stop=-1) == 4
  assert f(code="-+34", start=1, stop=3) == 3

def exercise_floating_point_scan_after_dot(f):
  try_code_none(f=f)
  assert f(code="") == 0
  assert f(code="0") == 1
  assert f(code="123456789") == 9
  assert f(code="e1") == 2
  assert f(code="e+1") == 3
  assert f(code="d-12") == 4
  assert f(code="0e+1") == 4
  assert f(code="0e") == -1
  assert f(code="456d003") == 7
  assert f(code="456d-0123456789") == 15
  assert f(code="xe+1", start=1) == 4
  assert f(code="xe+1") == 0
  assert f(code="xe+56", start=1, stop=-1) == 5
  assert f(code="xe+56", start=1, stop=4) == 4

def exercise_identifier_scan(f):
  try_code_none(f=f)
  assert f(code="") == -1
  assert f(code=")") == -1
  assert f(code="0") == -1
  assert f(code="a") == 1
  assert f(code="a(") == 1
  assert f(code="a_") == 2
  assert f(code="a_(") == 2
  assert f(code="a0_") == 3
  assert f(code="a0_=") == 3
  assert f(code="a0", start=1) == -1
  assert f(code="x*y*z_(i)", start=4) == 6
  assert f(code="_abcdefghijklmnopqrstuvwxyz0123456789=0") == 37
  assert f(code="*abc", start=1, stop=-1) == 4
  assert f(code="*abc", start=1, stop=3) == 3

def exercise_find_closing_parenthesis(f):
  try_code_none(f=f)
  assert f(code="") == -1
  assert f(code=")") == 0
  assert f(code="x)") == 1
  assert f(code="x)", start=1) == 1
  assert f(code="x)", start=2) == -1
  assert f(code="x)", start=2000) == -1
  assert f(code="()") == -1
  assert f(code="())") == 2
  assert f(code="(())") == -1
  assert f(code="(()))") == 4
  assert f(code="a(b(c(d(e),f(g),h(0:3)),g(4),h,j))") == -1
  assert f(code="a(b(c(d(e),f(g),h(0:3)),g(4),h,j)))") == 34
  assert f(code="x)", start=1, stop=-1) == 1
  assert f(code="x)", start=1, stop=1) == -1

def exercise_fem_utils_int_types():
  int_sizes = fable.exercise_fem_utils_int_types()
  expected = [1, 2, 4, 8]
  if (int_sizes != expected):
    hpp = "fem/utils/int_sizes.hpp"
    print "FATAL: %s: sizes are" % hpp, int_sizes, "but should be", expected
    raise RuntimeError(
      "%s needs to be adjusted for this platform." % hpp)

def exercise_fem_real_types():
  real_sizes = fable.exercise_fem_real_types()
  expected = [4, 8, 16]
  if (real_sizes != expected):
    hpp = "fem/data_types_star.hpp"
    print "FATAL: %s: sizes are" % hpp, real_sizes, "but should be", expected
    raise RuntimeError(
      "%s needs to be adjusted for this platform." % hpp)

def exercise_fem_format_tokenizer():
  f = fable.exercise_fem_format_tokenizer
  try: f("")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: empty string"
  else: raise Exception_expected
  try: f("(")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: ("
  else: raise Exception_expected
  try: f("x")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: x"
  else: raise Exception_expected
  try: f("  x  ")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: x"
  else: raise Exception_expected
  try: f("'")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: '"
  else: raise Exception_expected
  try: f("''")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: ''"
  else: raise Exception_expected
  try: f("'''")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: '''"
  else: raise Exception_expected
  try: f("''''")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: ''''"
  else: raise Exception_expected
  try: f("' '")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: ' '"
  else: raise Exception_expected
  try: f('" "')
  except RuntimeError, e:
    assert str(e) == 'Invalid FORMAT specification: " "'
  else: raise Exception_expected
  assert f("()") == []
  assert f(" (\t ) \t") == []
  assert f("(x)") == [("format", "x")]
  assert f(" ( X ) ") == [("format", "x")]
  assert f(" ( X , X ) ") == [("format", "x")]*2
  assert f("((:,/,$))") == [
    ("op", "("),
    ("op", ":"),
    ("op", "/"),
    ("op", "$"),
    ("op", ")")]
  assert f(" ( '' ) ") == [("string", "")]
  assert f(" ( '''' ) ") == [("string", "'")]
  assert f(" ( ''' ' ) ") == [("string", "' ")]
  assert f(" ( ' ''' ) ") == [("string", " '")]
  assert f(" ( ' '' ' ) ") == [("string", " ' ")]
  assert f(' ( """""" ) ') == [("string", '""')]
  assert f(' ( " X "" \'\t" ) ') == [("string", ' X " \'\t')]
  assert f(" ( + 1 P ) ") == [("format", "+1p")]
  assert f(" ( - 2 3 p ) ") == [("format", "-23p")]
  try: f(" ( + 1 X ) ")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: (+1x)"
  else: raise Exception_expected
  try: f(" ( 8 H ) ")
  except RuntimeError, e:
    assert str(e) == \
      "FATAL: Not supported: FORMAT Hollerith edit descriptor: (8h)"
  else: raise Exception_expected
  assert f(" ( 3 4 x ) ") == [("format", "34x")]
  assert f(" ( 5 6 7 p ) ") == [("format", "567p")]
  assert f(" ( 2 3 4 5 ) ") == [("integer", "2345")]
  assert f("(D)") == [("format", "d")]
  assert f("(E)") == [("format", "e")]
  assert f("(F)") == [("format", "f")]
  assert f("(G)") == [("format", "g")]
  assert f("(I)") == [("format", "i")]
  assert f("(Z)") == [("format", "z")]
  assert f("(D20)") == [("format", "d20")]
  assert f("(D30.10)") == [("format", "d30.10")]
  try: f("(D30.)")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: (d30.)"
  else: raise Exception_expected
  assert f("(A)") == [("format", "a")]
  assert f("(A37)") == [("format", "a37")]
  assert f("(L)") == [("format", "l")]
  assert f("(L46)") == [("format", "l46")]
  assert f("( B n )") == [("format", "bn")]
  assert f("( b z )") == [("format", "bz")]
  assert f("(S)") == [("format", "s")]
  assert f("( S P )") == [("format", "sp")]
  assert f("( s S )") == [("format", "ss")]
  try: f("(T)")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: (t)"
  else: raise Exception_expected
  try: f("(TL)")
  except RuntimeError, e:
    assert str(e) == "Invalid FORMAT specification: (tl)"
  else: raise Exception_expected
  assert f("(T12)") == [("format", "t12")]
  assert f("(TL23)") == [("format", "tl23")]
  assert f("(TR34)") == [("format", "tr34")]
  assert f("(21A83)") == [("integer", "21"), ("format", "a83")]
  assert f("(32('Uq'))") == [
    ("integer", "32"), ("op", "("), ("string", "Uq"), ("op", ")")]
  assert f("""\
(/1X,'Xk aqkIe',/1X,'UxKw',3(/,1X,3F12.6),/,' RxpO',/,\
        3F8.2,3F7.2,\
     /1X,'McR',3(/,1X,3F12.6),/1X,'HPi',\
     'EBx,Qyn ',2F6.1,' Qs>:-I ',3F6.2,/1X,'We,JkC,HIoR,Nb,JRG',\
     'RgMhBn JHyUK<>:lK',6F9.6,/1X,'K^dT ',3(/,1X,3F12.6))""") == [
    ("op", "/"),
    ("format", "1x"),
    ("string", "Xk aqkIe"),
    ("op", "/"),
    ("format", "1x"),
    ("string", "UxKw"),
    ("integer", "3"),
    ("op", "("),
    ("op", "/"),
    ("format", "1x"),
    ("integer", "3"),
    ("format", "f12.6"),
    ("op", ")"),
    ("op", "/"),
    ("string", " RxpO"),
    ("op", "/"),
    ("integer", "3"),
    ("format", "f8.2"),
    ("integer", "3"),
    ("format", "f7.2"),
    ("op", "/"),
    ("format", "1x"),
    ("string", "McR"),
    ("integer", "3"),
    ("op", "("),
    ("op", "/"),
    ("format", "1x"),
    ("integer", "3"),
    ("format", "f12.6"),
    ("op", ")"),
    ("op", "/"),
    ("format", "1x"),
    ("string", "HPi"),
    ("string", "EBx,Qyn "),
    ("integer", "2"),
    ("format", "f6.1"),
    ("string", " Qs>:-I "),
    ("integer", "3"),
    ("format", "f6.2"),
    ("op", "/"),
    ("format", "1x"),
    ("string", "We,JkC,HIoR,Nb,JRG"),
    ("string", "RgMhBn JHyUK<>:lK"),
    ("integer", "6"),
    ("format", "f9.6"),
    ("op", "/"),
    ("format", "1x"),
    ("string", "K^dT "),
    ("integer", "3"),
    ("op", "("),
    ("op", "/"),
    ("format", "1x"),
    ("integer", "3"),
    ("format", "f12.6"),
    ("op", ")")]

def compare_floats(v, val):
  from math import frexp
  vm, ve = frexp(val)
  scale = 2.0**(-ve)
  delta = abs(v*scale - val*scale)
  return delta < 1.e-14

def exercise_fem_utils_string_to_double():
  stream_end = 256
  err_eoi = "End of input while reading floating-point value."
  err_inv = "Invalid character while reading floating-point value: "
  err_oor = "Out-of-range error while reading floating-point value."
  f = fable.exercise_fem_utils_string_to_double
  def check(str, err=None):
    v,e,n = f(str=str)
    if (err is not None):
      assert e == err
      assert v == 0
    else:
      assert e is None
      val = float(str)
      assert compare_floats(v, val)
      assert n == stream_end
      vx,e,n = f(str=str+"x")
      assert e is None
      assert n == ord("x")
    return v,e,n
  for str in ["", " ", "\t", " \t"]:
    check(str, err=err_eoi)
  check("0")
  check("1")
  check("12")
  check(".1")
  check("1.")
  check("-1")
  check("-.1")
  check("-1.")
  check("+1.")
  check("0.1")
  check(".01")
  check("0.01")
  check("0.012")
  check("0.0012")
  check("00.0012")
  check("1234567890")
  check("12345678901234567890")
  check("123456789012345678901234567890")
  check(".12")
  check("1.2")
  check("12.")
  check("012.")
  check("12.34")
  check("123.4")
  check("1234.")
  check("2.3e0")
  check("2.345e1")
  check("2.345e2")
  check("2.345e3")
  check("2.345e4")
  check("2.345e5")
  check("2.345e+305")
  check("2.345e-1")
  check("2.345e-2")
  check("2.345e-3")
  check("2.345e-4")
  check("2.345e-5")
  check("2.345e-305")
  v,e,n = check("000.0000000000000000000002345e-5")
  assert "%.3f" % (v * 1e27) == "2.345"
  v,e,n = f("2.345e45")
  assert v == f("2.345E45")[0]
  assert v == f("2.345d45")[0]
  assert v == f("2.345D45")[0]
  v,e,n = f("1e305.")
  def chk():
    assert v != 0
    assert e is None
    assert n == ord(".")
  chk()
  v,e,n = f(".0001e309.")
  chk()
  v,e,n = f("1e-305.")
  chk()
  v,e,n = f("10000e-309.")
  chk()
  v,e,n = f("1e3091")
  def chk():
    assert v == 0
    assert e == err_oor
    assert n == ord("1")
  chk()
  v,e,n = f("1e-3091")
  chk()
  v,e,n = f("xy")
  assert v == 0
  assert e == err_inv + '"x" (ordinal=120)'
  assert n == ord("y")
  v,e,n = f("3e")
  assert v == 0
  assert e == err_eoi
  assert n == stream_end
  v,e,n = f("3eyz")
  assert v == 0
  assert e == err_inv + '"y" (ordinal=121)'
  assert n == ord("z")
  v,e,n = f("eq")
  assert v == 0
  assert e == err_inv + '"e" (ordinal=101)'
  assert n == ord("q")
  v,e,n = f('"\'')
  assert v == 0
  assert e == err_inv + "'\"' (double quote, ordinal=34)"
  assert n == ord("'")
  v,e,n = f("'\"")
  assert v == 0
  assert e == err_inv + '"\'" (single quote, ordinal=39)'
  assert n == ord('"')
  v,e,n = f(chr(134)+'x')
  assert v == 0
  assert e == err_inv + 'ordinal=134'
  assert n == ord("x")

def exercise_fem_utils_string_to_double_fmt():
  f = fable.exercise_fem_utils_string_to_double_fmt
  def check(str, w, d, blanks_zero, exp_scale, val):
    v,e,n = f(str=str, w=w, d=d, blanks_zero=blanks_zero, exp_scale=exp_scale)
    assert e is None
    assert n == ord("4")
    assert compare_floats(v, val)
  check("   4", 3, 0, False, 0, 0.)
  check("   4", 3, 0, True, 0, 0.)
  #
  # compare with test/valid/read_fmt_double.f
  check("1234", 3, 0, False, 0, 123.)
  check("1234", 3, 1, False, 0, 12.3)
  check("1234", 3, 2, False, 0, 1.23)
  check(".234", 3, 2, False, 0, 0.23)
  check("1.34", 3, 2, False, 0, 1.3)
  check("12.4", 3, 2, False, 0, 12.)
  check("1 34", 3, 0, False, 0, 13.)
  check("1 34", 3, 0, True, 0, 103.)
  check("1234", 3, 0, False, 1, 12.3)
  check("1234", 3, 0, False, 2, 1.23)
  check("1234", 3, 0, False, -1, 1230.)
  check("1234", 3, 0, False, -2, 12300.)
  check("1e34", 3, 0, False, 1, 1000.)
  check("3 5 e 1 24", 9, 0, False, 0, 35000000000000.)
  check("3 5 e 1 24", 9, 0, True, 0, 3.05e105)
  check("3 5 e 1 24", 9, 1, False, 0, 3500000000000.)
  check("3 5 e 1 24", 9, 1, True, 0, 3.05e104)

def run(args):
  assert len(args) == 0
  for f in [fable.py_unsigned_integer_scan,
               fable.unsigned_integer_scan]:
    exercise_unsigned_integer_scan(f)
  for f in [fable.py_floating_point_scan_after_exponent_char,
               fable.floating_point_scan_after_exponent_char]:
    exercise_floating_point_scan_after_exponent_char(f)
  for f in [fable.py_floating_point_scan_after_dot,
               fable.floating_point_scan_after_dot]:
    exercise_floating_point_scan_after_dot(f)
  for f in [fable.py_identifier_scan,
               fable.identifier_scan]:
    exercise_identifier_scan(f)
  for f in [fable.py_find_closing_parenthesis,
               fable.find_closing_parenthesis]:
    exercise_find_closing_parenthesis(f)
  if (fable.ext is not None):
    exercise_fem_utils_int_types()
    exercise_fem_real_types()
    exercise_fem_format_tokenizer()
    exercise_fem_utils_string_to_double()
    exercise_fem_utils_string_to_double_fmt()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
