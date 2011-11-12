from cctbx.sgtbx import subgroups
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx import uctbx
from libtbx import complex_math
from libtbx.utils import format_cpu_times
from libtbx.test_utils import Exception_expected, approx_equal, \
  not_approx_equal, show_diff
import libtbx.load_env
import math
import weakref
import pickle
from cStringIO import StringIO
import os

ad_hoc_1992_pairs = """\
Abm2 Aem2
Bma2 Bme2
B2cm B2em
C2mb C2me
Cm2a Cm2e
Ac2m Ae2m
Aba2 Aea2
Bba2 Bbe2
B2cb B2eb
C2cb C2ce
Cc2a Cc2e
Ac2a Ae2a
Cmca Cmce
Ccmb Ccme
Abma Aema
Acam Aeam
Bbcm Bbem
Bmab Bmeb
Cmma Cmme
Abmm Aemm
Bmcm Bmem
Ccca Ccce
Abaa Aeaa
Bbcb Bbeb
""".splitlines()

def exercise_symbols():
  s = sgtbx.space_group_symbols("p 2")
  assert s.number() == 3
  assert s.schoenflies() == "C2^1"
  assert s.qualifier() == "b"
  assert s.hermann_mauguin() == "P 1 2 1"
  assert s.extension() == "\0"
  assert s.change_of_basis_symbol() == ""
  assert s.universal_hermann_mauguin() == "P 1 2 1"
  assert s.hall() == " P 2y"
  s = sgtbx.space_group_symbols("p 2", "A")
  assert s.hall() == " P 2y"
  s = sgtbx.space_group_symbols("p 2", "I")
  assert s.hall() == " P 2"
  s = sgtbx.space_group_symbols(146)
  assert s.hall() == " R 3"
  s = sgtbx.space_group_symbols(146, "r")
  assert s.hall() == " P 3*"
  s = sgtbx.space_group_symbols(146, "", "A1983")
  assert s.hall() == " R 3"
  s = sgtbx.space_group_symbols(146, "", "I1952")
  assert s.hall() == " P 3*"
  assert s.point_group_type() == "3"
  assert s.laue_group_type() == "-3"
  assert s.crystal_system() == "Trigonal"
  o = StringIO()
  for hm in ["P2(1)", "P3112", "P3212", "P6122", "P6522", "P6222", "P6422",
             "Pnnn:1"]:
    s = sgtbx.space_group_symbols("%s (z, x, y)" % hm)
    assert s.change_of_basis_symbol() == "z,x,y"
    assert s.universal_hermann_mauguin().endswith(" (z,x,y)")
    print >> o, s.hall()
  assert not show_diff(o.getvalue(), """\
 P 2yb (z,x,y)
 P 31 2 (z+1/3,x,y)
 P 32 2 (z+1/6,x,y)
 P 61 2 (z+5/12,x,y)
 P 65 2 (z+1/12,x,y)
 P 62 2 (z+1/3,x,y)
 P 64 2 (z+1/6,x,y)
 P 2 2 -1n (z,x,y)
""")
  s = sgtbx.space_group_symbols("P222(1)")
  assert s.hall() == " P 2c 2"
  s = sgtbx.space_group_symbols("P 31 1 2 (x,y,z-1/3)")
  assert s.hall() == " P 31 2"
  symbols = {
    "  hall p 3  ": "p 3  ",
    " hall:   p 4": "p 4",
    " hall :   p 6": "p 6",
    " p    1 21/n     1": "-P 2yn",
    "C  6  v  #  1": " P 6 -2",
    "1  2": "-C 2y",
  }
  for symbol,hall in symbols.items():
    assert sgtbx.space_group_symbols(symbol).hall() == hall, symbol
  i = sgtbx.space_group_symbol_iterator()
  assert id(i) == id(iter(i))
  assert i.next().universal_hermann_mauguin() == "P 1"
  assert i.next().universal_hermann_mauguin() == "P -1"
  assert len(tuple(i)) == 528
  i = sgtbx.space_group_symbol_iterator()
  n = 0
  for s in i:
    n += 1
    hm = s.universal_hermann_mauguin()
    assert sgtbx.space_group_symbols(hm).universal_hermann_mauguin() == hm
  assert n == 530
  #
  short_symbols_i = [
    "P2", "P21", "B2", "C2", "Pm", "Pb", "Pc", "Bm", "Cm", "Bb", "Cc", "P2/m",
    "P21/m", "B2/m", "C2/m", "P2/b", "P2/c", "P21/b", "P21/c", "B2/b", "C2/c"]
  expected_numbers_i=[3,4,5,5,6,7,7,8,8,9,9,10,11,12,12,13,13,14,14,15,15]
  assert [sgtbx.space_group_symbols(short, "I").number()
    for short in short_symbols_i] == expected_numbers_i
  short_symbols_a = [
   "P2", "P21", "C2", "A2", "I2", "Pm", "Pc", "Pn", "Pa", "Cm", "Am", "Im",
   "Cc", "An", "Ia", "Aa", "Cn", "Ic", "P2/m", "P21/m", "C2/m", "A2/m",
   "I2/m", "P2/c", "P2/n", "P2/a", "P21/c", "P21/n", "P21/a", "C2/c", "A2/n",
   "I2/a", "A2/a", "C2/n", "I2/c"]
  expected_numbers_a = [
    3,4,5,5,5,6,7,7,7,8,8,8,9,9,9,9,9,9,10,11,
    12,12,12,13,13,13,14,14,14,15,15,15,15,15,15]
  assert [sgtbx.space_group_symbols(short, "A").number()
    for short in short_symbols_a] == expected_numbers_a
  symbols_cpp = libtbx.env.under_dist("cctbx", "sgtbx/symbols.cpp")
  if (not os.path.isfile(symbols_cpp)):
    print "Skipping checks based on %s: file not available" % symbols_cpp
  else:
    f = iter(open(symbols_cpp))
    for volume,table_name,short_symbols in [
          ("I", "vol_i_short_mono_hm_dict", short_symbols_i),
          ("A", "vol_a_short_mono_hm_dict", short_symbols_a)]:
      for line in f:
        if (line.find(table_name) > 0): break
      else: raise AssertionError
      symbol_pairs = []
      for line in f:
        if (line.find("{ 0, 0 },") > 0): break
        for c in '{},"': line = line.replace(c,"")
        symbol_pairs.append(line.split())
      else: raise AssertionError
      assert len(symbol_pairs) == len(short_symbols)
      for expected_short,(short,long) in zip(short_symbols,symbol_pairs):
        assert expected_short == short
        assert sgtbx.space_group_symbols(short, volume).hall() \
            == sgtbx.space_group_symbols(long, volume).hall()
  #
  s = sgtbx.space_group_symbols
  try: s("(x,y,z)")
  except RuntimeError, e:
    assert str(e) == "cctbx Error: Space group symbol not recognized: (x,y,z)"
  else: raise Exception_expected
  try: s("P3:2")
  except RuntimeError, e:
    assert str(e) == "cctbx Error: Space group symbol not recognized: P3:2"
  else: raise Exception_expected
  try: s(300)
  except RuntimeError, e:
    assert str(e) == "cctbx Error: Space group number out of range: 300"
  else: raise Exception_expected
  try: s(space_group_number=1, table_id="x")
  except RuntimeError, e:
    assert str(e) == "cctbx Error: table_id not recognized: x"
  else: raise Exception_expected
  for extension in ["1", ":1"]:
    try: s(space_group_number=75, extension=extension)
    except RuntimeError, e:
      assert str(e) == "cctbx Error: Space group symbol not recognized: 75:1"
    else: raise Exception_expected
  try: s(space_group_number=75, extension="x")
  except RuntimeError, e:
    assert str(e) == "cctbx Error: Space group symbol not recognized: 75:x"
  else: raise Exception_expected
  #
  for cabc,cxyz in [("2/3a+1/3b+1/3c, -1/3a+1/3b+1/3c, -1/3a-2/3b+1/3c",
                     "x+z,-x+y+z,-y+z"),
                    ("-1/3a-2/3b+1/3c, 2/3a+1/3b+1/3c, -1/3a+1/3b+1/3c",
                     "x+z,-x+y+z,-y+z")]:
    symbols = s("R 3 m :h (%s)" % cabc)
    assert str(sgtbx.space_group_info(group=sgtbx.space_group(symbols))) \
        == "R 3 m :R"
  #
  a83_symbols = []
  for symbol in "Aem2 Aea2 Cmce Cmme Ccce CCCE:1".split():
    a83_symbols.append(str(sgtbx.space_group_info(symbol=symbol)))
  assert a83_symbols == [
    "A b m 2", "A b a 2", "C m c a", "C m m a", "C c c a :2", "C c c a :1"]
  for pair in ad_hoc_1992_pairs:
    o,n = [sgtbx.space_group_info(symbol=symbol) for symbol in pair.split()]
    assert o.group() == n.group()
  #
  def check(symbols, uhm):
    for symbol in symbols:
      s = sgtbx.space_group_symbols(symbol=symbol, table_id="")
      assert s.universal_hermann_mauguin() == uhm
  check(["R3", "H3", " h 3 "], "R 3 :H")
  check(["R32", "H32", "_h_3_2_"], "R 3 2 :H")

def exercise_tr_vec():
  tr_vec = sgtbx.tr_vec
  t = tr_vec()
  assert t.num() == (0,0,0)
  assert t.den() == sgtbx.sg_t_den
  t = tr_vec(13)
  assert t.num() == (0,0,0)
  assert t.den() == 13
  t = tr_vec((1,2,3))
  assert t.num() == (1,2,3)
  assert t.den() == sgtbx.sg_t_den
  t = tr_vec((1,2,3), 13)
  assert t.num() == (1,2,3)
  assert t.den() == 13
  assert t == t
  assert not t != t
  s = tr_vec((1,2,3), 12)
  assert not t == s
  assert t != s
  s = tr_vec((1,2,4), 13)
  assert not t == s
  assert t != s
  assert t.is_valid()
  assert not t.is_zero()
  t = tr_vec(0)
  assert not t.is_valid()
  assert t.is_zero()
  t = tr_vec((2,4,-2)).new_denominator(6)
  assert t.den() == 6
  assert t.num() == (1,2,-1)
  t = t.scale(3)
  assert t.den() == 18
  assert t.num() == (3,6,-3)
  t = t.mod_positive()
  assert t.den() == 18
  assert t.num() == (3,6,15)
  t = t.mod_short()
  assert t.den() == 18
  assert t.num() == (3,6,-3)
  t = t.cancel()
  assert t.den() == 6
  assert t.num() == (1,2,-1)
  assert t.as_double() == (1./6,2./6,-1./6)
  a = tr_vec((1,2,3), 12)
  b = tr_vec((9,0,1), 12)
  c = a.plus(b)
  assert c.den() == 6
  assert c.num() == (5,1,2)
  c = a.minus(b)
  assert c.den() == 6
  assert c.num() == (-4,1,1)
  a = tr_vec((-1,0,2), 4)
  assert str(a) == "-1/4,0,1/2"
  assert a.as_string() == "-1/4,0,1/2"
  assert a.as_string(True) == "-.25,0,.5"
  assert a.as_string(False, ";") == "-1/4;0;1/2"

def exercise_rot_mx():
  tr_vec = sgtbx.tr_vec
  rot_mx = sgtbx.rot_mx
  rot_mx_info = sgtbx.rot_mx_info
  r = rot_mx()
  assert r.den() == 1
  assert r.num() == (1,0,0,0,1,0,0,0,1)
  r = rot_mx(2)
  assert r.den() == 2
  assert r.num() == (2,0,0,0,2,0,0,0,2)
  r = rot_mx(2, 3)
  assert r.den() == 2
  assert r.num() == (6,0,0,0,6,0,0,0,6)
  r = rot_mx((1,0,-1,1,1,0,0,-1,-1))
  assert r.den() == 1
  assert r.num() == (1,0,-1,1,1,0,0,-1,-1)
  r = rot_mx((1,0,-1,1,1,0,0,-1,-1), 2)
  assert r.den() == 2
  assert r.num() == (1,0,-1,1,1,0,0,-1,-1)
  assert r == r
  assert not r != r
  s = rot_mx((1,0,-1,1,1,0,0,-1,-1), 3)
  assert not r == s
  assert r != s
  s = rot_mx((1,0,-1,2,1,0,0,-1,-1), 2)
  assert not r == s
  assert r != s
  assert r.is_valid()
  assert not r.is_unit_mx()
  r = rot_mx(0)
  assert not r.is_valid()
  r = rot_mx()
  assert r.is_unit_mx()
  r = rot_mx(2, 1)
  assert r.is_unit_mx()
  r = rot_mx().minus_unit_mx()
  assert r.den() == 1
  assert r.num() == (0,0,0,0,0,0,0,0,0)
  r = rot_mx((1,0,-1,1,1,0,0,-1,-1)).minus_unit_mx()
  assert r.den() == 1
  assert r.num() == (0,0,-1,1,0,0,0,-1,-2)
  r = rot_mx(12, 3).new_denominator(4)
  assert r.den() == 4
  assert r.num() == (12,0,0,0,12,0,0,0,12)
  r = r.scale(3)
  assert r.den() == 12
  assert r.num() == (36,0,0,0,36,0,0,0,36)
  assert r.determinant() == 27
  r = rot_mx(range(9), 2).transpose()
  assert r.num() == (0, 3, 6, 1, 4, 7, 2, 5, 8)
  assert r.den() == 2
  assert rot_mx().inverse().num() == rot_mx().num()
  r3 = (0,-1,0,1,-1,0,0,0,1)
  r = rot_mx(r3).inverse()
  assert r.den() == 1
  assert r.num() == (-1,1,0,-1,0,0,0,0,1)
  assert rot_mx(r3).inverse().inverse().num() == r3
  r = rot_mx((3,0,0,0,3,0,0,0,3), 6).cancel()
  assert r.den() == 2
  assert r.num() == (1,0,0,0,1,0,0,0,1)
  r = rot_mx((0,-1,1,1,0,1,-1,1,1)).scale(12)
  s = r.inverse()
  assert s.den() == 12
  assert s.num() == (-4,8,-4,-8,4,4,4,4,4)
  s = r.inverse_cancel()
  assert s.den() == 3
  assert s.num() == (-1,2,-1,-2,1,1,1,1,1)
  assert rot_mx().multiply(rot_mx()).is_unit_mx()
  assert rot_mx(r3).multiply(rot_mx(r3)) == rot_mx(r3).inverse()
  assert rot_mx(r3).multiply(rot_mx(r3).multiply(rot_mx(r3))).is_unit_mx()
  assert rot_mx().multiply(tr_vec((1,2,3), 12)) == tr_vec((1,2,3), 12)
  assert rot_mx().multiply(tr_vec((2,4,6), 12)) == tr_vec((1,2,3), 6)
  assert rot_mx(r3).multiply(tr_vec((1,2,3), 8)) == tr_vec((-2,-1,3), 8)
  assert rot_mx(r3).multiply(tr_vec((4,-4,2), 8)) == tr_vec((2,4,1), 4)
  assert rot_mx((6,0,0,0,6,0,0,0,6), 3).divide(2) == rot_mx()
  assert rot_mx(r3).divide(2).as_double() == (0,-1./2,0,1./2,-1./2,0,0,0,1./2)
  assert rot_mx(r3).type() == 3
  assert rot_mx(r3).order() == 3
  assert rot_mx(r3).order(3) == 3
  assert rot_mx(1, -1).type() == -1
  assert rot_mx(1, -1).order() == 2
  assert rot_mx(1, -1).order(-1) == 2
  assert rot_mx().accumulate() == rot_mx((1,0,0,0,1,0,0,0,1))
  assert rot_mx((-1,0,0,0,-1,0,0,0,-1)).accumulate() \
      == rot_mx((0,0,0,0,0,0,0,0,0))
  assert rot_mx(r3).accumulate() == rot_mx((0,0,0,0,0,0,0,0,3))
  assert rot_mx((0,1,0,-1,0,0,0,0,-1)).accumulate() \
      == rot_mx((0,0,0,0,0,0,0,0,0))
  i = rot_mx_info(rot_mx())
  assert i.type() == 1
  assert i.ev() == (0,0,0)
  assert i.basis_of_invariant() == ((1,0,0), (0,1,0), (0,0,1))
  assert i.sense() == 0
  i = rot_mx().info()
  assert i.type() == 1
  assert i.ev() == (0,0,0)
  assert i.basis_of_invariant() == ((1,0,0), (0,1,0), (0,0,1))
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,0,0,0,-1,0,0,0,1)))
  assert i.type() == 2
  assert i.ev() == (0,0,1), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,0,0,0,1,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (0,1,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((1,0,0,0,-1,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (1,0,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  i = rot_mx_info(rot_mx((1,0,0,0,1,0,0,0,-1)))
  assert i.type() == -2
  assert i.ev() == (0,0,1), i.ev()
  assert i.basis_of_invariant() == ((1,0,0), (0,1,0))
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((1,0,0,0,-1,0,0,0,1)))
  assert i.type() == -2
  assert i.ev() == (0,1,0), i.ev()
  assert i.basis_of_invariant() == ((1,0,0), (0,0,1))
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,0,0,0,1,0,0,0,1)))
  assert i.type() == -2
  assert i.ev() == (1,0,0), i.ev()
  assert i.basis_of_invariant() == ((0,1,0), (0,0,1))
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((0,-1,0,-1,0,0,0,0,1)))
  assert i.type() == -2
  assert i.ev() == (1,1,0), i.ev()
  assert i.basis_of_invariant() == ((0,0,1), (-1,1,0))
  i = rot_mx_info(rot_mx((1,0,0,0,0,1,0,1,0)))
  assert i.type() == -2
  assert i.ev() == (0,-1,1), i.ev()
  assert i.basis_of_invariant() == ((1,0,0), (0,1,1))
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((0,-1,0,1,-1,0,0,0,1)))
  assert i.type() == 3
  assert i.ev() == (0,0,1), i.ev()
  assert i.basis_of_invariant() == (i.ev(), )
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,-1,0,1,-1,0,0,0,1)).inverse())
  assert i.type() == 3
  assert i.ev() == (0,0,1), i.ev()
  assert i.basis_of_invariant() == (i.ev(), )
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,1,0,-1,1,0,0,0,-1)))
  assert i.type() == -3
  assert i.ev() == (0,0,1), i.ev()
  assert i.basis_of_invariant() == ()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,1,0,-1,1,0,0,0,-1)).inverse())
  assert i.type() == -3
  assert i.ev() == (0,0,1), i.ev()
  assert i.basis_of_invariant() == ()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,0,1,0,1,0,-1,0,0)))
  assert i.type() == 4
  assert i.ev() == (0,1,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,0,1,0,1,0,-1,0,0)).inverse())
  assert i.type() == 4
  assert i.ev() == (0,1,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,0,-1,0,-1,0,1,0,0)))
  assert i.type() == -4
  assert i.ev() == (0,1,0), i.ev()
  assert i.basis_of_invariant() == ()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,0,-1,0,-1,0,1,0,0)).inverse())
  assert i.type() == -4
  assert i.ev() == (0,1,0), i.ev()
  assert i.basis_of_invariant() == ()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((1,0,0,0,1,-1,0,1,0)))
  assert i.type() == 6
  assert i.ev() == (1,0,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((1,0,0,0,1,-1,0,1,0)).inverse())
  assert i.type() == 6
  assert i.ev() == (1,0,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((-1,0,0,0,-1,1,0,-1,0)))
  assert i.type() == -6
  assert i.ev() == (1,0,0), i.ev()
  assert i.basis_of_invariant() == ()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((-1,0,0,0,-1,1,0,-1,0)).inverse())
  assert i.type() == -6
  assert i.ev() == (1,0,0), i.ev()
  assert i.basis_of_invariant() == ()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,0,1,1,0,0,0,1,0)))
  assert i.type() == 3
  assert i.ev() == (1,1,1), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,0,1,1,0,0,0,1,0)).inverse())
  assert i.type() == 3
  assert i.ev() == (1,1,1), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,-1,0,-1,0,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (-1,1,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,1,0,0,1,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (1,2,0), i.ev()
  assert i.basis_of_invariant() == (i.ev(),)
  r = rot_mx((-1,1,0,0,1,0,0,0,-1))
  assert r.as_xyz() == "-x+y,y,-z"
  assert r.as_hkl() == "-h,h+k,-l"
  r = rot_mx([-2,3,5,1,4,-2,-3,2,-1],6)
  assert approx_equal(r * [0.2,0.1,-0.5], [-0.4333333, 0.2666667, 0.01666667])
  assert approx_equal([0.2,0.1,-0.5] * r, [0.2, 0, 0.21666667])
  assert str(r * sgtbx.vec3_rat_from_str("1/5,1/10,-1/2")) \
      == "(-13/30, 4/15, 1/60)"

def exercise_rt_mx():
  tr_vec = sgtbx.tr_vec
  rot_mx = sgtbx.rot_mx
  rt_mx = sgtbx.rt_mx
  s = rt_mx()
  assert s.r().den() == 1
  assert s.t().den() == sgtbx.sg_t_den
  s = rt_mx(2)
  assert s.r().den() == 2
  assert s.t().den() == sgtbx.sg_t_den
  s = rt_mx(2, 3)
  assert s.r().den() == 2
  assert s.t().den() == 3
  r = rot_mx((-1,1,0,0,1,0,0,0,-1))
  t = tr_vec((1,2,3))
  s = rt_mx(r, t)
  assert s.r() == r
  assert s.t() == t
  s = rt_mx(r)
  assert s.r() == r
  assert s.t() == tr_vec()
  s = rt_mx(r, 3)
  assert s.r() == r
  assert s.t() == tr_vec(3)
  s = rt_mx(t)
  assert s.r() == rot_mx()
  assert s.t() == t
  s = rt_mx(t, 3)
  assert s.r() == rot_mx(3)
  assert s.t() == t
  p = sgtbx.parse_string("x,y,z")
  assert p.string() == "x,y,z"
  s = rt_mx(p)
  assert s.is_unit_mx()
  p = sgtbx.parse_string("-x,-y+1/2,z")
  s = rt_mx(p)
  assert s.as_xyz() == "-x,-y+1/2,z"
  p = sgtbx.parse_string("-x,-y+1/2,z;")
  s = rt_mx(p, ";")
  assert s.as_xyz() == "-x,-y+1/2,z"
  p = sgtbx.parse_string("-  x,- y + 1/2, z ;")
  s = rt_mx(p, ";", 2)
  assert s.r().den() == 2
  assert s.as_xyz() == "-x,-y+1/2,z"
  p = sgtbx.parse_string("-  x,-1* y + .5, z ;")
  s = rt_mx(p, ";", 2, 4)
  assert s.r().den() == 2
  assert s.t().den() == 4
  assert s.as_xyz() == "-x,-y+1/2,z"
  s = rt_mx("1/2+x,-y,z")
  assert s.as_xyz() == "x+1/2,-y,z"
  s = rt_mx("x-y,x,z+5/6#", "#")
  assert s.as_xyz() == "x-y,x,z+5/6"
  s = rt_mx("-y,x-y,z+.66666666#", "#", 3)
  assert s.r().den() == 3
  assert s.as_xyz() == "-y,x-y,z+2/3"
  s = rt_mx("y,y-x,z+1/6#", "#", 3, 6)
  assert s.r().den() == 3
  assert s.t().den() == 6
  assert s.as_xyz() == "y,-x+y,z+1/6"
  assert rt_mx("y,y-x,z+1/6") == rt_mx("y,-x+y,z+1/6")
  assert not rt_mx("y,y-x,z+1/6") != rt_mx("y,-x+y,z+1/6")
  assert not rt_mx("y,y-x,z-1/6") == rt_mx("y,-x+y,z+1/6")
  assert rt_mx("y,y-x,z-1/6") != rt_mx("y,-x+y,z+1/6")
  assert rt_mx().is_valid()
  assert not rt_mx(0).is_valid()
  assert rt_mx("-x,-y,-z", "", 2, 3).unit_mx() == rt_mx(2, 3)
  assert rt_mx().is_unit_mx()
  assert not rt_mx("-x,-y,-z").is_unit_mx()
  assert rt_mx("0.667+x, 0.333+y, 0.333+z").as_xyz() == 'x+2/3,y+1/3,z+1/3'
  s = rt_mx("y,y-x,z+1/4", "", 3, 8)
  assert str(s) == "y,-x+y,z+1/4"
  assert s.as_xyz() == str(s)
  assert s.as_xyz(False) == str(s)
  assert s.as_xyz(True) == "y,-x+y,z+.25"
  assert s.as_xyz(False, False) == str(s)
  assert s.as_xyz(False, True) == "y,-x+y,1/4+z"
  assert s.as_xyz(False, False, "xyz") == str(s)
  assert s.as_xyz(False, False, "abc") == "b,-a+b,c+1/4"
  assert s.as_xyz(False, False, "xyz", ",") == str(s)
  assert s.as_xyz(False, False, "xyz", " ; ") == "y ; -x+y ; z+1/4"
  assert s.as_int_array() == (0,3,0,-3,3,0,0,0,3,0,0,2)
  assert s.as_double_array() == (0,1,0,-1,1,0,0,0,1,0,0,0.25)
  assert s.new_denominators(6).r().den() == 6
  assert s.new_denominators(9,16).r().den() == 9
  assert s.new_denominators(9,16).t().den() == 16
  assert s.new_denominators(rt_mx("x,y,z", "", 9, 16)).r().den() == 9
  assert s.new_denominators(rt_mx("x,y,z", "", 9, 16)).t().den() == 16
  assert rt_mx("x-1/2,y+4/3,z-1/3").mod_positive() \
      == rt_mx("x+1/2,y+1/3,z+2/3")
  assert rt_mx("x-1/2,y+5/3,z-2/3").mod_short() \
      == rt_mx("x+1/2,y-1/3,z+1/3")
  assert rt_mx("x-y,x,z+5/6").inverse() == rt_mx("y,-x+y,z-5/6")
  s = rt_mx("y,y-x,z+1/6")
  assert s.refine_gridding((1,1,1)) == (1,1,6)
  assert rt_mx("y,y-x,z+1/6", "", 2, 24).cancel() \
      == rt_mx("y,y-x,z+1/6", "", 1, 6)
  assert rt_mx("x-y,x,z+5/6").inverse_cancel() \
      != rt_mx("y,-x+y,z-5/6")
  assert rt_mx("x-y,x,z+5/6").inverse_cancel() \
      == rt_mx("y,-x+y,z-5/6").cancel()
  assert rt_mx("x-y,x,z+5/6").multiply(rt_mx("y,-x+y,z-5/6")).is_unit_mx()
  i = sgtbx.translation_part_info(rt_mx())
  assert i.intrinsic_part().is_zero()
  assert i.location_part().is_zero()
  assert i.origin_shift().is_zero()
  i = sgtbx.translation_part_info(rt_mx("-y,x-y,z+1/3"))
  assert i.intrinsic_part().as_double() == (0,0,1./3)
  assert i.location_part().is_zero()
  assert i.origin_shift().is_zero()
  i = sgtbx.translation_part_info(rt_mx("-x+y,-x+1/4,z+2/3"))
  assert i.intrinsic_part().as_double() == (0,0,2./3)
  assert i.location_part().as_double() == (0,-1./4,0)
  assert i.origin_shift().as_double() == (1./12,1./6,0)
  assert approx_equal(rt_mx("z,x,y") * (2,3,4), (4,2,3))
  assert str(rt_mx("z,x,y") * sgtbx.vec3_rat_from_str("2,-1/7,10/-4")) \
      == "(-5/2, 2, -1/7)"
  r = (-1,1,0,0,1,0,0,0,-1)
  t = (1/12.,2/12.,3/12.)
  assert str(rt_mx(r, t)) == "-x+y+1/12,y+1/6,-z+1/4"
  assert str(rt_mx(r, t, 2)) == "-x+y+1/12,y+1/6,-z+1/4"
  assert str(rt_mx(r, t, 2, 24)) == "-x+y+1/12,y+1/6,-z+1/4"
  assert str(rt_mx("+0*x-1*y+0*z,+0*x+0*y+1*z,-1*x+0*y-1*z")) == "-y,z,-x-z"
  assert str(rt_mx("-x+y,-x+1/4,z+2/3") + (2,-1,3)) == "-x+y+2,-x-3/4,z+11/3"
  for m in [rt_mx("x,y,z"), rt_mx(0,0), rt_mx("-x+y+2,-x-3/4,z+11/3")]:
    p = pickle.dumps(m)
    l = pickle.loads(p)
    assert l == m
  v = sgtbx.stl_vector_rt_mx()
  v.append(rt_mx("-x+y+2,-x-3/4,z+11/3"))
  v.append(rt_mx("x,y,z"))
  assert [str(elem) for elem in v] == ["-x+y+2,-x-3/4,z+11/3", "x,y,z"]
  p = pickle.dumps(v)
  l = pickle.loads(p)
  assert [str(elem) for elem in l] == ["-x+y+2,-x-3/4,z+11/3", "x,y,z"]
  l.extend(v)
  assert l.size() == 4
  l.clear()
  assert l.size() == 0
  assert str(rt_mx("z,x,y") + (1,2,-3)) == "z+1,x+2,y-3"
  assert str(rt_mx("z+1/6,x,y") + tr_vec((1,2,-3),12)) == "z+1/4,x+1/6,y-1/4"
  assert str(rt_mx("-x+y,y,-z")) == "-x+y,y,-z"
  assert str(rt_mx("-h,h+k,-l")) == "-x+y,y,-z"
  try: rt_mx("h,x,z")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: mix of x,y,z and h,k,l notation:
  h,x,z
  __^""")
  else: raise Exception_expected
  try: rt_mx("")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: unexpected end of input:
  \n\
  ^""")
  else: raise Exception_expected
  try: rt_mx("x, ")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: unexpected end of input:
  x, \n\
  ___^""")
  else: raise Exception_expected
  try: rt_mx("x")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: not enough row expressions:
  x
  _^""")
  else: raise Exception_expected
  try: rt_mx("x,y,x,z")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: too many row expressions:
  x,y,x,z
  _____^""")
  try: rt_mx("x++")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: unexpected character:
  x++
  __^""")
  try: rt_mx("a, b, c")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: a,b,c notation not supported in this context:
  a, b, c
  ^""")
  else: raise Exception_expected
  #
  s = rt_mx("y-13/2,y-x+8/3,z+1/4")
  site_frac_1 = (4.2, -5.2, 8.9)
  site_frac_2 = (-2.1, -7.2, 12.9)
  assert s.unit_shifts_minimum_distance(
    site_frac_1=site_frac_1,
    site_frac_2=site_frac_2) == (18, -3, -4)
  su = s.add_unit_shifts_minimum_distance(
    site_frac_1=site_frac_1,
    site_frac_2=site_frac_2)
  assert str(su) == "y+23/2,-x+y-1/3,z-15/4"
  assert approx_equal(su*site_frac_2, (4.3, -5.4333333, 9.15))

def exercise_change_of_basis_op():
  rt_mx = sgtbx.rt_mx
  change_of_basis_op = sgtbx.change_of_basis_op
  c = change_of_basis_op(rt_mx(), rt_mx())
  assert c.is_identity_op()
  c = change_of_basis_op(rt_mx())
  assert c.is_identity_op()
  c = change_of_basis_op("z,x,y")
  c = change_of_basis_op("z,x,y;junk", ";")
  c = change_of_basis_op("z,x,y;junk", ";", 3)
  c = change_of_basis_op("z,x,y;junk", ";", 3, 4)
  assert not c.is_identity_op()
  assert c.c().r().den() == 3
  assert c.c().t().den() == 4
  assert c.c_inv().r().den() == 3
  assert c.c_inv().t().den() == 4
  assert str(c.c()) == "z,x,y"
  assert str(c.c_inv()) == "y,z,x"
  d = c.identity_op()
  assert d.is_identity_op()
  assert d.c().r().den() == 3
  assert d.c().t().den() == 4
  assert d.c_inv().r().den() == 3
  assert d.c_inv().t().den() == 4
  assert str(d.c()) == "x,y,z"
  assert str(d.c_inv()) == "x,y,z"
  assert d.is_valid()
  d = change_of_basis_op(0, 0)
  assert not d.is_valid()
  d = c.new_denominators(4, 5)
  assert d.c().r().den() == 4
  assert d.c().t().den() == 5
  d = d.new_denominators(c)
  assert d.c().r().den() == 3
  assert d.c().t().den() == 4
  c = change_of_basis_op("z,x,y")
  assert c.select(False) == c.c()
  assert c.select(True) == c.c_inv()
  d = c.inverse()
  assert c.c() == d.c_inv()
  assert d.c() == c.c_inv()
  c = change_of_basis_op("x+3/2,y-5/4,z+1/3")
  c.mod_positive_in_place()
  assert str(c.c()) == "x+1/2,y+3/4,z+1/3"
  assert str(c.mod_short().c()) == "x+1/2,y-1/4,z+1/3"
  assert str(c.mod_short().c_inv()) == "x+1/2,y+1/4,z-1/3"
  c.mod_short_in_place()
  assert str(c.c()) == "x+1/2,y-1/4,z+1/3"
  assert str(c.c_inv()) == "x+1/2,y+1/4,z-1/3"
  c = change_of_basis_op("z,x,y")
  assert str(c.apply(rt_mx("-x,-y,z"))) == "x,-y,-z"
  xf = c((0.1,0.2,0.3))
  for i in xrange(3): assert approx_equal(xf[i], (0.3,0.1,0.2)[i])
  c.update(change_of_basis_op("z,x,y"))
  assert str(c.c()) == "y,z,x"
  assert (c * c.inverse()).is_identity_op()
  c = change_of_basis_op("z,x,y")
  u = uctbx.unit_cell((2,3,5))
  assert approx_equal(c.apply(u).parameters(), (5,2,3,90,90,90))
  assert c.apply((1,2,3)) == (3,1,2)
  i = flex.miller_index(((1,2,3),(2,3,4)))
  assert tuple(c.apply(i)) == ((3,1,2),(4,2,3))
  assert c.apply_results_in_non_integral_indices(miller_indices=i).size() == 0
  d = change_of_basis_op("-x,-y,2z")
  assert list(d.apply_results_in_non_integral_indices(miller_indices=i)) == [0]
  i = i.select(flex.size_t([1,0]))
  assert list(d.apply_results_in_non_integral_indices(miller_indices=i)) == [1]
  s = pickle.dumps(c)
  l = pickle.loads(s)
  assert str(c.c()) == str(l.c())
  for s in ["-x+y,y,-z", "-h,h+k,-l", "-a,a+b,-c"]:
    c = change_of_basis_op(s)
    assert c.as_xyz() == "-x+y,y,-z"
    assert c.as_hkl() == "-h,h+k,-l"
    assert c.as_abc() == "-a,a+b,-c"
  c = change_of_basis_op("-a+1/8,a+b-1/3,-c+2/3")
  assert c.as_xyz() == "-x+y+11/24,y+1/3,-z+2/3"
  assert c.as_abc() == "-a+1/8,a+b-1/3,-c+2/3"
  assert c.inverse().as_xyz() == "-x+y+1/8,y-1/3,-z+2/3"
  a = c.c().as_4x4_rational()
  b = c.c_inv().as_4x4_rational()
  assert str(a.elems).replace(" ","") \
      == "(-1,1,0,11/24,0,1,0,1/3,0,0,-1,2/3,0,0,0,1)"
  assert (a*b).elems == (1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1)
  assert str(b.transpose().elems).replace(" ","") \
      == "(-1,0,0,0,1,1,0,0,0,0,-1,0,1/8,-1/3,2/3,1)"
  assert c.symbol() == "-a+1/8,a+b-1/3,-c+2/3"
  assert c.inverse().symbol() == "-x+y+1/8,y-1/3,-z+2/3"
  assert str(c) == c.symbol()
  assert str(c.inverse()) == c.inverse().symbol()

def exercise_space_group():
  tr_vec = sgtbx.tr_vec
  rt_mx = sgtbx.rt_mx
  space_group = sgtbx.space_group
  g = space_group()
  p = sgtbx.parse_string("P 1")
  g = space_group(p)
  p = sgtbx.parse_string("P 1")
  g = space_group(p, False)
  p = sgtbx.parse_string("P 1")
  g = space_group(p, False, False)
  p = sgtbx.parse_string("P 1")
  g = space_group(p, False, False, False)
  assert g.t_den() == sgtbx.sg_t_den
  p = sgtbx.parse_string("P 1")
  g = space_group(p, False, False, False, 2*sgtbx.sg_t_den)
  assert g.t_den() == 2*sgtbx.sg_t_den
  g = space_group("P 1")
  g = space_group("P 1", True)
  g = space_group("1", False, True)
  g = space_group("1", False, True, True)
  g = space_group("1", False, True, True, 3*sgtbx.sg_t_den)
  assert g.t_den() == 3*sgtbx.sg_t_den
  g = space_group(sgtbx.space_group_symbols(1))
  g = space_group(g)
  g.reset()
  assert g.r_den() == 1
  assert g.t_den() == sgtbx.sg_t_den
  g.reset(6)
  assert g.r_den() == 1
  assert g.t_den() == 6
  g.reset()
  assert g.expand_ltr(tr_vec()) is g
  assert g.n_ltr() == 1
  g.reset()
  g.expand_ltr(tr_vec((6,0,0)))
  assert g.n_ltr() == 2
  g.reset()
  g.expand_ltr(tr_vec((4,0,0)))
  assert g.n_ltr() == 3
  g.reset()
  assert g.f_inv() == 1
  assert not g.is_centric()
  assert not g.is_origin_centric()
  assert g.expand_inv(tr_vec((0,0,0))) is g
  assert g.f_inv() == 2
  assert g.is_centric()
  assert g.is_origin_centric()
  g.reset()
  g.expand_inv(tr_vec((1,2,3)))
  assert g.is_centric()
  assert not g.is_origin_centric()
  g.reset()
  assert g.n_smx() == 1
  assert g.expand_smx(rt_mx("-x,-y,z")) is g
  assert g.n_smx() == 2
  g = sgtbx.space_group()
  assert g.expand_smx("-x,y,-z") is g
  assert not g.is_tidy()
  assert g.make_tidy() is g
  assert g.is_tidy()
  assert g.order_z() == 2
  g.expand_smx("-x,-y,-z")
  assert g.order_z() == 4
  assert not g.is_tidy()
  for s in g:
    assert g.contains(smx=s)
  assert not g.contains(smx=sgtbx.rt_mx("y,z,x"))
  for z,n in (("P",1), ("a",2), ("B",2),("c",2), ("I",2), ("r",3), ("F",4)):
    g.reset()
    g.expand_conventional_centring_type(z)
    assert g.n_ltr() == n
  p = sgtbx.parse_string("P 4")
  g.reset()
  g.parse_hall_symbol(p)
  assert g.order_z() == 4
  assert len(g) == 4
  assert g.n_equivalent_positions() == 4
  p = sgtbx.parse_string("P 2x")
  g.parse_hall_symbol(p, True)
  assert g.order_z() == 8
  p = sgtbx.parse_string("-1")
  g.parse_hall_symbol(p, False, True)
  assert g.order_z() == 16
  g = space_group("P 4")
  c = sgtbx.change_of_basis_op("x,y,z")
  h = g.change_basis(c)
  assert g == h
  c = sgtbx.change_of_basis_op("z,x,y")
  h = g.change_basis(c)
  assert g != h
  g = space_group("P 4x")
  assert g == h
  g = space_group("P 6")
  assert g.order_p() == 6
  assert g.order_z() == 6
  g = space_group("-F 2 2")
  assert g.order_p() == 8
  assert g.order_z() == 32
  j = 0
  for i_ltr in xrange(g.n_ltr()):
    for i_inv in xrange(g.f_inv()):
      for i_smx in xrange(g.n_smx()):
        assert g(i_ltr, i_inv, i_smx) == g(j)
        assert g[j] == g(j)
        j += 1
  assert len(tuple(g)) == g.order_z()
  g = space_group("P 3")
  g.expand_smx(rt_mx("-x,-y,z"))
  h = space_group("P 2")
  h.expand_smx(rt_mx("-x+y,-x,z"))
  gx = [str(s) for s in g]
  hx = [str(s) for s in h]
  assert gx != hx
  g.make_tidy()
  h.make_tidy()
  gx = [str(s) for s in g]
  hx = [str(s) for s in h]
  assert gx == hx
  assert g == h
  assert not g != h
  for z in "PABCIRF":
    assert space_group(z + " 1").conventional_centring_type_symbol() == z
  assert space_group("P 1").z2p_op().is_identity_op()
  assert not space_group("R 1").z2p_op().is_identity_op()
  assert not space_group("R 1").construct_z2p_op().is_identity_op()
  assert space_group("P 3").is_chiral()
  assert not space_group("-P 3").is_chiral()
  g = space_group("P 41 (1,2,-1)")
  m = flex.miller_index(((0,0,1), (0,0,4), (1,0,0)))
  assert g.is_sys_absent((0,0,1))
  assert not g.is_sys_absent((0,0,4))
  assert tuple(g.is_sys_absent(m)) == (True,False,False)
  assert g.is_centric((1,0,0))
  assert not g.is_centric((0,0,4))
  assert tuple(g.is_centric(m)) == (False,False,True)
  p = g.phase_restriction((1,0,0))
  assert not p.sys_abs_was_tested()
  assert p.ht() == 2
  assert g.is_valid_phase((1,0,0), p.ht_angle())
  assert g.is_valid_phase((1,0,0), p.ht_angle(), False)
  assert g.is_valid_phase((1,0,0), p.ht_angle(True), True)
  assert not g.is_valid_phase((1,0,0), p.ht_angle()+math.pi/180, False)
  assert g.is_valid_phase((1,0,0), p.ht_angle()+math.pi/180, False, 1.e6)
  assert approx_equal(g.nearest_valid_phases(
    miller_indices=flex.miller_index([(1,2,3),(1,0,0)]),
    phases=flex.double([0.123,0.234])), [0.123, 0.52359877559829882])
  assert approx_equal(g.nearest_valid_phases(
    miller_indices=flex.miller_index([(1,2,3),(1,0,0)]),
    phases=flex.double([123,234]),
    deg=True), [123, 210])
  assert g.multiplicity((1,2,3), False) == 8
  assert g.multiplicity((1,2,3), True) == 4
  assert tuple(g.multiplicity(m, False)) == (2,2,4)
  assert tuple(g.multiplicity(m, True)) == (1,1,4)
  assert g.epsilon((1,2,3)) == 1
  assert tuple(g.epsilon(m)) == (4,4,1)
  gm = space_group("F 2 2 -1d")
  def check(s, m):
    assert gm.multiplicity(sgtbx.vec3_rat_from_str(s)) == m
  check("1/4,1/4,1/4", 8)
  check("1/4,1/4,-1/4", 8)
  check("1/8,1/8,1/8", 16)
  check("5/8,5/8,5/8", 16)
  check("1/9,1/4,1/4", 16)
  check("1/4,1/9,1/4", 16)
  check("1/4,1/4,1/9", 16)
  check("1/7,1/5,1/9", 32)
  u = uctbx.unit_cell((3, 3, 4, 90, 90, 120))
  assert u.is_similar_to(space_group("P 6").average_unit_cell(u))
  assert space_group("P 6").is_compatible_unit_cell(u)
  assert space_group("P 6").is_compatible_unit_cell(u, 0.01)
  assert space_group("P 6").is_compatible_unit_cell(u, 0.01, 1)
  assert not space_group("P 3*").is_compatible_unit_cell(u)
  assert space_group("P 3*").average_unit_cell(u).is_similar_to(
    uctbx.unit_cell((3.3665, 3.3665, 3.3665, 97.6056, 97.6056, 97.6056)),
    1.e-5, 1.e-3)
  assert space_group("P 3*").is_compatible_unit_cell(u, 1, 30)
  u = uctbx.unit_cell((95.2939, 95.2939, 98.4232, 94.3158, 115.226, 118.822))
  g = space_group("C 2y (x+y,-x+y+z,z)")
  assert g.is_compatible_unit_cell(u)
  assert approx_equal(g.average_u_star(u_star=range(6,0,-1)),
    (6.5, 5.5, 4.0, 3.5, 2.5, 1.5))
  g = space_group("C 2 -2c")
  h = g.build_derived_reflection_intensity_group(anomalous_flag=True)
  assert h == space_group("C 2 -2")
  h = h.build_derived_acentric_group()
  assert h == space_group("C 2 -2")
  h = g.build_derived_reflection_intensity_group(anomalous_flag=False)
  assert h == space_group("-C 2 2")
  h = h.build_derived_acentric_group()
  assert h == space_group("C 2 2")
  h = g.build_derived_patterson_group()
  assert h == space_group("-C 2 2")
  h = g.build_derived_point_group()
  assert h == space_group("P 2 -2")
  h = g.build_derived_laue_group()
  assert h == space_group("-P 2 2")
  assert g.point_group_type() == "mm2"
  assert g.laue_group_type() == "mmm"
  assert g.crystal_system() == "Orthorhombic"
  for s in sgtbx.space_group_symbol_iterator():
    g = space_group(s)
    m = g.match_tabulated_settings()
    assert s.number() == m.number()
    assert not g.build_derived_acentric_group().is_centric()
  g = space_group("P 31")
  assert g.gridding() == (1,1,3)
  assert g.refine_gridding((2,1,4)) == (2,2,12)
  mod_p = ["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3"]
  mod_s = ["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z-1/3"]
  assert [str(s) for s in g.all_ops()] == mod_p
  assert [str(s) for s in g.all_ops(0)] == mod_p
  assert [str(s) for s in g.all_ops(1)] == mod_p
  assert [str(s) for s in g.all_ops(-1)] == mod_s
  assert [str(s) for s in g.all_ops(-1, False)] == mod_s
  assert [str(s) for s in g.all_ops(-1, True)] == mod_s
  assert g.all_ops(0, False)[1].t().den() == sgtbx.sg_t_den
  assert g.all_ops(0, True)[1].t().den() == 3
  g = space_group("P 31")
  assert g.type().number() == 144
  p = pickle.dumps(g)
  l = pickle.loads(p)
  assert g == l
  #
  g = space_group("-P 4 2")
  for c in ["c-1/3,a+1/4,b-1/8", "z-1/3,x+1/2,y-1/4"]:
    c = sgtbx.change_of_basis_op(c)
    for cc in [c.as_xyz(), c.as_abc()]:
      h = space_group("-P 4 2 (%s)" % cc)
      cc = sgtbx.change_of_basis_op(cc)
      assert g.change_basis(cc) == h
      assert h.change_basis(cc.inverse()) == g
      for ccc in [cc.inverse().as_xyz(), cc.inverse().as_abc()]:
        ccc = sgtbx.change_of_basis_op(ccc)
        assert h.change_basis(ccc) == g
  cxyz = "2/3*x-1/3*y-1/3*z,1/3*x+1/3*y-2/3*z,1/3*x+1/3*y+1/3*z"
  chkl = "h-k,k-l,h+k+l"
  cabc = "a-b,b-c,a+b+c"
  for c in [cxyz, chkl, cabc]:
    assert sgtbx.change_of_basis_op(c).as_xyz() == cxyz
    assert sgtbx.change_of_basis_op(c).as_hkl() == chkl
    assert sgtbx.change_of_basis_op(c).as_abc() == cabc
  r = space_group("P 3* -2") # R 3 m :r
  h = space_group('R 3 -2"') # R 3 m :h
  for c in [cxyz, chkl, cabc]:
    t = space_group("P 3* -2 (%s)" % c)
    c = sgtbx.change_of_basis_op(c)
    assert r.change_basis(c) == t
    assert t == h
    assert c.inverse().as_xyz() == "x+z,-x+y+z,-y+z"
    assert c.inverse().as_hkl() \
        == "2/3*h+1/3*k+1/3*l,-1/3*h+1/3*k+1/3*l,-1/3*h-2/3*k+1/3*l"
    assert c.inverse().as_abc() \
        == "2/3*a+1/3*b+1/3*c,-1/3*a+1/3*b+1/3*c,-1/3*a-2/3*b+1/3*c"
    for ci in [c.inverse().as_xyz(),
               c.inverse().as_hkl(),
               c.inverse().as_abc()]:
      t = space_group('R 3 -2" (%s)' % ci)
      ci = sgtbx.change_of_basis_op(ci)
      assert h.change_basis(ci) == t
      assert t == r
  #
  try: space_group("-P 4 2 (p)")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Parse error: unexpected character:
  -P 4 2 (p)
  ________^""")
  else: raise Exception_expected
  try: space_group("-P 4 2 (x,x,x)")
  except RuntimeError, e:
    assert not show_diff(str(e), """\
cctbx Error: Rotation matrix is not invertible.""")
  else: raise Exception_expected
  sg = space_group("-P 2ac 2ab")
  cb = sgtbx.change_of_basis_op("x+1/12,y+1/12,z-1/12")
  sg1 = sg.change_basis(cb)
  assert sg1.is_centric() and not sg1.is_origin_centric()
  icb = sg1.change_of_origin_realising_origin_centricity()
  sg2 = sg1.change_basis(icb)
  assert sg2.is_origin_centric()
  assert str(cb * icb) == "a,b,c+1/2"
  #
  sg = space_group()
  s = sgtbx.rt_mx(2, 12)
  assert s.r().den() == 2
  assert s.t().den() == 12
  try: sg.expand_smx(s)
  except RuntimeError, e:
    assert str(e) == "cctbx Error: sgtbx::space_group::expand_smx():" \
      " rotation-part denominator must be 1 (implementation limitation)."
  else: raise Exception_expected
  s = sgtbx.rt_mx(1, 24)
  assert s.r().den() == 1
  assert s.t().den() == 24
  try: sg.expand_smx(s)
  except RuntimeError, e:
    assert str(e) == "cctbx Error: sgtbx::space_group::expand_smx():" \
      " incompatible translation-part denominator."
  else: raise Exception_expected

def exercise_space_group_type():
  space_group = sgtbx.space_group
  space_group_type = sgtbx.space_group_type
  t = space_group_type()
  assert t.group() == space_group()
  assert t.number() == 1
  assert t.cb_op().is_identity_op()
  t = space_group_type("P 2")
  assert t.hall_symbol() == " P 2y"
  assert t.universal_hermann_mauguin_symbol() == "P 1 2 1"
  t = space_group_type("P 2", "a1983")
  assert t.hall_symbol() == " P 2y"
  assert t.universal_hermann_mauguin_symbol() == "P 1 2 1"
  t = space_group_type("P 2", "i1952")
  assert t.hall_symbol() == " P 2y (z,x,y)"
  assert t.universal_hermann_mauguin_symbol() == "P 1 2 1 (c,a,b)"
  t = space_group_type("P 3")
  assert len(t.addl_generators_of_euclidean_normalizer(False, False)) == 0
  assert len(t.addl_generators_of_euclidean_normalizer(True, False)) == 1
  assert len(t.addl_generators_of_euclidean_normalizer(False, True)) == 2
  assert len(t.addl_generators_of_euclidean_normalizer(True, True)) == 3
  g = t.expand_addl_generators_of_euclidean_normalizer(False, False)
  assert g == t.group()
  j = t.expand_addl_generators_of_euclidean_normalizer(True, False).type()
  assert j.lookup_symbol() == "P -3"
  j = t.expand_addl_generators_of_euclidean_normalizer(False, True).type()
  assert j.lookup_symbol() == "P 6 m m"
  j = t.expand_addl_generators_of_euclidean_normalizer(True, True).type()
  assert j.lookup_symbol() == "P 6/m m m"
  assert not space_group_type("P 3").is_enantiomorphic()
  assert space_group_type("P 31").is_enantiomorphic()
  assert space_group_type("P -1").change_of_hand_op().is_identity_op()
  assert str(space_group_type("P 1").change_of_hand_op().c()) == "-x,-y,-z"
  assert str(space_group_type("P 31").change_of_hand_op().c()) == "-x,-y,-z"
  assert str(space_group_type("I41").change_of_hand_op().c()) == "-x+1/2,-y,-z"
  g = space_group(sgtbx.space_group_symbols("C c c a :1"))
  t = sgtbx.space_group_type(group=g, tidy_cb_op=False)
  assert not t.cb_op_is_tidy()
  assert t.hall_symbol(tidy_cb_op=False) == "-C 2a 2ac (x+1/2,-y+1/4,-z-1/4)"
  assert t.hall_symbol(tidy_cb_op=True) == "-C 2a 2ac (x+1/2,y-1/4,z+1/4)"
  assert t.hall_symbol() == "-C 2a 2ac (x+1/2,y-1/4,z+1/4)"
  assert t.universal_hermann_mauguin_symbol(tidy_cb_op=False) \
      == "C c c a :2 (a+1/2,-b+1/4,-c-1/4)"
  assert t.universal_hermann_mauguin_symbol(tidy_cb_op=True) \
      == "C c c a :2 (a+1/2,b+1/4,c-1/4)"
  assert t.universal_hermann_mauguin_symbol() \
      == "C c c a :2 (a+1/2,b+1/4,c-1/4)"
  p = pickle.dumps(t)
  l = pickle.loads(p)
  assert l.group() == t.group()
  assert not l.cb_op_is_tidy()
  t = g.change_basis(g.z2p_op()).type()
  assert t.cb_op_is_tidy()
  h = "-C 2a 2ac (x-y-1/4,x+y+1/4,z+1/4)"
  assert t.hall_symbol() == h
  assert t.hall_symbol(tidy_cb_op=True) == h
  assert t.hall_symbol(tidy_cb_op=False) == h
  u = "C c c a :2 (x-y-1/4,x+y+1/4,z+1/4)"
  assert t.universal_hermann_mauguin_symbol() == u
  assert t.universal_hermann_mauguin_symbol(tidy_cb_op=True) == u
  assert t.universal_hermann_mauguin_symbol(tidy_cb_op=False) == u
  assert t.lookup_symbol() == u
  assert g.type().lookup_symbol() == "C c c a :1"
  p = pickle.dumps(t)
  l = pickle.loads(p)
  assert l.group() == t.group()
  assert l.cb_op_is_tidy()
  #
  for s in sgtbx.space_group_symbol_iterator():
    uhm_in = s.universal_hermann_mauguin() + " (z+1/4, x-1/4, y+3/4)"
    g_out = sgtbx.space_group(sgtbx.space_group_symbols(uhm_in))
    uhm_out = g_out.type().universal_hermann_mauguin_symbol()
    if (uhm_out != uhm_in):
      g_out2 = sgtbx.space_group(sgtbx.space_group_symbols(uhm_out))
      assert g_out2 == g_out
      uhm_out2 = g_out2.type().universal_hermann_mauguin_symbol()
      assert uhm_out2 == uhm_out
    if (uhm_out.find("(") >= 0):
      ehm, cb = uhm_out.split("(")
      assert cb[-1] == ")"
      g_out2 = sgtbx.space_group(sgtbx.space_group_symbols(ehm)) \
        .change_basis(sgtbx.change_of_basis_op(cb[:-1]))
      assert g_out2 == g_out
    g_out = g_out.change_basis(g_out.z2p_op()) # primitive setting
    uhm_out = g_out.type().universal_hermann_mauguin_symbol()
    g_out2 = sgtbx.space_group(sgtbx.space_group_symbols(uhm_out))
    assert g_out2 == g_out
    uhm_out2 = g_out2.type().universal_hermann_mauguin_symbol()
    assert uhm_out2 == uhm_out
    if (uhm_out.find("(") >= 0):
      ehm, cb = uhm_out.split("(")
      assert cb[-1] == ")"
      g_out2 = sgtbx.space_group(sgtbx.space_group_symbols(ehm)) \
        .change_basis(sgtbx.change_of_basis_op(cb[:-1]))
      assert g_out2 == g_out
  #
  # ensure universal_hermann_mauguin_symbol ":H" is consistently upper-case
  t = sgtbx.space_group_type("r -3 M :h (3*z,-x+2*z,-y+z)")
  assert t.universal_hermann_mauguin_symbol() == "R -3 m :H (3*z,-x+2*z,-y+z)"
  assert t.lookup_symbol() == "R -3 m :H (3*z,-x+2*z,-y+z)"
  t = sgtbx.space_group_type("r -3 M :h")
  assert t.universal_hermann_mauguin_symbol() == "R -3 m :H"
  assert t.lookup_symbol() == "R -3 m :H"
  #
  for pair in [line.split() for line in ad_hoc_1992_pairs]:
    o,n = [sgtbx.space_group_info(symbol=symbol) for symbol in pair]
    for ad_hoc_1992 in [False, True]:
      l1 = o.type().lookup_symbol(ad_hoc_1992=ad_hoc_1992)
      l2 = n.type().lookup_symbol(ad_hoc_1992=ad_hoc_1992)
      assert not show_diff(l1, l2)
      assert l1.replace(" ", "").startswith(pair[int(ad_hoc_1992)])
  sgi = sgtbx.space_group_info(symbol="Ae2a").primitive_setting()
  assert not show_diff(
    sgi.type().lookup_symbol(), "A b a 2 (x,y-z,y+z)")
  assert not show_diff(
    sgi.type().lookup_symbol(ad_hoc_1992=True), "A e a 2 (x,y-z,y+z)")
  #
  def all_t_zero(g): # this simple test is correct only for reference settings
    assert g.make_tidy()
    for s in g.smx():
      if (not s.t().is_zero()): return False
    return True
  n_symmorphic = 0
  for sgi in sgtbx.reference_space_group_infos():
    is_symmorphic = sgi.type().is_symmorphic()
    assert is_symmorphic == all_t_zero(sgi.group())
    if (is_symmorphic):
      n_symmorphic += 1
  assert n_symmorphic == 73

def exercise_phase_info():
  phase_info = sgtbx.phase_info
  g = sgtbx.space_group("P 61 (1 2 0)")
  p = phase_info(g, (1,2,3))
  assert p.sys_abs_was_tested()
  p = phase_info(g, (1,2,3), False)
  assert p.sys_abs_was_tested()
  p = phase_info(g, (1,2,3), True)
  assert not p.sys_abs_was_tested()
  p = phase_info(g, (0,0,6))
  assert not p.is_sys_absent()
  p = phase_info(g, (0,0,1))
  assert p.is_sys_absent()
  p = phase_info(g, (1,0,0))
  assert p.is_centric()
  p = phase_info(g, (0,0,2))
  assert not p.is_centric()
  p = phase_info(g, (1,0,0))
  assert p.ht() == 2
  assert p.t_den() == g.t_den()
  assert approx_equal(p.ht_angle(), float(p.ht()) / p.t_den() * math.pi)
  assert p.ht_angle(False) == p.ht_angle()
  assert approx_equal(p.ht_angle(True), p.ht_angle()*180/math.pi)
  p = phase_info(sgtbx.space_group("P 2ac 2ab"), (0,10,8))
  assert approx_equal(p.nearest_valid_phase(1.e-15), 0)
  assert approx_equal(p.nearest_valid_phase(-1.e-15), 0)
  assert approx_equal(p.nearest_valid_phase(math.pi/2-1.e-6), 0)
  assert approx_equal(p.nearest_valid_phase(math.pi/2+1.e-6), math.pi)
  assert approx_equal(p.nearest_valid_phase(math.pi-1.e-15), math.pi)
  assert approx_equal(p.nearest_valid_phase(math.pi+1.e-15), math.pi)
  p = phase_info(g, (1,0,0))
  phi = p.ht_angle()
  assert approx_equal(p.nearest_valid_phase(phi+1.e-15), phi)
  assert approx_equal(p.nearest_valid_phase(phi-1.e-15), phi)
  assert approx_equal(p.nearest_valid_phase(phi+math.pi/2-1.e-6), phi)
  assert approx_equal(p.nearest_valid_phase(phi+math.pi/2+1.e-6), phi+math.pi)
  assert approx_equal(p.nearest_valid_phase(phi+math.pi+1.e-15), phi+math.pi)
  assert approx_equal(p.nearest_valid_phase(phi+math.pi-1.e-15), phi+math.pi)
  for i in xrange(-3, 4):
    phi = p.ht_angle() + i * math.pi
    assert p.is_valid_phase(phi)
    assert p.is_valid_phase(phi*180/math.pi, True)
    assert p.is_valid_phase(phi, False, 1.e-6)
    phi = p.ht_angle() + i * math.pi + math.pi / 180
    assert not p.is_valid_phase(phi)
    assert not p.is_valid_phase(phi*180/math.pi, True)
    assert not p.is_valid_phase(phi, False, 1.e-6)
    assert p.is_valid_phase(phi, False, 1.e6)
    assert p.is_valid_phase(p.nearest_valid_phase(phi))
    assert p.is_valid_phase(p.nearest_valid_phase(phi, False), False)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+1.e-15, False), False)
    assert p.is_valid_phase(p.nearest_valid_phase(phi-1.e-15, False), False)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+math.pi+1.e-15, False), False)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+math.pi-1.e-15, False), False)
    assert p.is_valid_phase(p.nearest_valid_phase(phi, True), True)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+1.e-15, True), True)
    assert p.is_valid_phase(p.nearest_valid_phase(phi-1.e-15, True), True)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+180+1.e-15, True), True)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+180-1.e-15, True), True)
  for i in xrange(-3, 4):
    phi = p.ht_angle() + i * math.pi
    f = complex_math.polar((1, phi))
    assert approx_equal(f, p.valid_structure_factor(f))
    f = complex_math.polar((1, phi+math.pi/2))
    assert approx_equal(abs(p.valid_structure_factor(f)), 0)
    f = complex_math.polar((1, phi+math.pi/3))
    assert approx_equal(abs(p.valid_structure_factor(f)), 0.5)
    f = complex_math.polar((1, phi+math.pi/4))
    assert approx_equal(abs(p.valid_structure_factor(f)), math.sqrt(2)/2)
    f = complex_math.polar((1, phi+math.pi/6))
    assert approx_equal(abs(p.valid_structure_factor(f)), math.sqrt(3)/2)
    f = complex_math.polar((1, phi-math.pi/6))
    assert approx_equal(abs(p.valid_structure_factor(f)), math.sqrt(3)/2)
    f = complex_math.polar((2, phi+math.pi/6))
    assert approx_equal(abs(p.valid_structure_factor(f)), math.sqrt(3))
    f = complex_math.polar((3, phi+math.pi/6))
    assert approx_equal(abs(p.valid_structure_factor(f)), 3*math.sqrt(3)/2)

def exercise_reciprocal_space_asu():
  reciprocal_space_asu = sgtbx.reciprocal_space_asu
  t = sgtbx.space_group_type("P 1 2 1")
  a = reciprocal_space_asu(t)
  assert a.cb_op().c() == t.cb_op().c()
  assert a.is_reference()
  assert a.reference_as_string() == "k>=0 and (l>0 or (l==0 and h>=0))"
  assert a.is_inside((1,2,3))
  assert not a.is_inside((-1,-2,-3))
  assert a.which((1,2,3)) == 1
  assert a.which((-1,-2,-3)) == -1
  t = sgtbx.space_group_type("P 1 1 2")
  a = reciprocal_space_asu(t)
  assert not a.is_reference()
  t = sgtbx.space_group_type("P 3 1 2")
  a = reciprocal_space_asu(t)
  assert a.which((-1,0,1)) == 0

def exercise_brick():
  tr_vec = sgtbx.tr_vec
  brick = sgtbx.brick
  t = sgtbx.space_group_type("P 1 2 1")
  b = brick(t)
  assert b.as_string() == "0<=x<=1/2; 0<=y<1; 0<=z<1"
  assert b.as_string() == str(b)
  assert b.is_inside(tr_vec((6,2,3)))
  assert not b.is_inside(tr_vec((7,2,3)))

def exercise_site_symmetry():
  site_symmetry = sgtbx.site_symmetry
  u = uctbx.unit_cell((3,4,5,80,100,110))
  g = sgtbx.space_group("P 2")
  s = site_symmetry(u, g, (0,0,0))
  s = site_symmetry(u, g, (0,0,0), 0.5)
  s = site_symmetry(
    unit_cell=u,
    space_group=g,
    original_site=(0.05,0,0),
    min_distance_sym_equiv=0.5,
    assert_min_distance_sym_equiv=True)
  assert s.unit_cell().parameters() == u.parameters()
  assert s.space_group() == g
  assert s.original_site() == (0.05,0,0)
  assert approx_equal(s.min_distance_sym_equiv(), 0.5)
  assert s.exact_site() == (0,0,0)
  assert approx_equal(s.distance_moved(), u.distance((0.05,0,0), (0,0,0)))
  assert s.shortest_distance() > 0.5
  assert s.check_min_distance_sym_equiv()
  assert s.multiplicity() == 1
  assert str(s.special_op()) == "0,0,z"
  assert not s.is_point_group_1()
  assert s.point_group_type() == "2"
  assert [str(o) for o in s.unique_ops()] == ["0,0,z"]
  assert sgtbx.rt_mx("-x,-y,z") in s
  cb_op = sgtbx.change_of_basis_op("z,x,y")
  sc = s.change_basis(cb_op=cb_op)
  assert str(sc.special_op()) == "x,0,0"
  assert str(sc.matrices()[0]) == "x,y,z"
  assert str(sc.matrices()[1]) == "x,-y,-z"
  assert str(s.matrices()[1]) == "-x,-y,z"
  assert sgtbx.rt_mx("x,-y,-z") in sc
  assert sc.multiplicity() == 1
  cb_op = sgtbx.change_of_basis_op("-x,-y,-z") # negative determinant
  sc = s.change_basis(cb_op=cb_op)
  assert str(sc.special_op()) == "0,0,z"
  assert sc.multiplicity() == 1
  u = uctbx.unit_cell((10,10,15,90,90,120))
  g = sgtbx.space_group("-P 3 2")
  s = site_symmetry(unit_cell=u, space_group=g, original_site=(0,0,0))
  assert s.multiplicity() == 1
  a = (5,2,3,-1,2,-2)
  assert not s.is_compatible_u_star(u_star=a)
  assert s.is_compatible_u_star(u_star=a, tolerance=1.e6)
  a = s.average_u_star(u_star=a)
  assert s.is_compatible_u_star(u_star=a)
  assert len(s.matrices()) == 12
  assert str(s.matrices()[1]) == "-y,x-y,z"
  t = sgtbx.site_symmetry_table()
  t.reserve(5)
  assert t.indices().size() == 0
  t.process(site_symmetry_ops=s)
  assert t.indices().size() == 1
  assert t.n_special_positions() == 1
  assert list(t.special_position_indices()) == [0]
  assert t.n_unique() == 2
  t.process(site_symmetry_ops=s)
  assert t.indices().size() == 2
  assert t.n_unique() == 2
  s2 = site_symmetry(unit_cell=u, space_group=g, original_site=(0.5,0.5,0.5))
  assert s2.multiplicity() == 3
  t.process(site_symmetry_ops=s2)
  assert t.indices().size() == 3
  assert t.n_unique() == 3
  t.process(site_symmetry_ops=s)
  assert t.indices().size() == 4
  assert t.n_unique() == 3
  s3 = site_symmetry(unit_cell=u, space_group=g, original_site=(0.25,0.66,0.0))
  assert s3.multiplicity() == 12
  assert s3.is_point_group_1()
  assert s3 == s3
  assert not (s3 != s3)
  assert not (s3 == s2)
  assert s3 != s2
  t.process(site_symmetry_ops=s3)
  assert t.indices().size() == 5
  assert t.n_unique() == 3
  assert len(t.table()) == 3
  assert list(t.indices()) == [1,1,2,1,0]
  assert t.n_special_positions() == 4
  assert list(t.special_position_indices()) == [0,1,2,3]
  assert t.is_special_position(i_seq=0)
  assert not t.is_special_position(i_seq=4)
  assert t.get(0).multiplicity() == s.multiplicity()
  assert t.get(1).multiplicity() == s.multiplicity()
  assert t.get(2).multiplicity() == s2.multiplicity()
  assert t.get(3).multiplicity() == s.multiplicity()
  assert t.get(4).multiplicity() == s3.multiplicity()
  assert t.get(0).special_op() == s.special_op()
  assert t.get(1).special_op() == s.special_op()
  assert t.get(2).special_op() == s2.special_op()
  assert t.get(3).special_op() == s.special_op()
  assert t.get(4).special_op() == s3.special_op()
  assert t.get(0).n_matrices() == 12
  assert t.get(1).n_matrices() == 12
  assert t.get(2).n_matrices() == 4
  assert t.get(3).n_matrices() == 12
  assert t.get(4).n_matrices() == 1
  for i in xrange(5):
    assert t.get(i).n_matrices() * t.get(i).multiplicity() == g.order_z()
  assert not t.get(0).is_point_group_1()
  assert t.get(4).is_point_group_1()
  assert str(t.get(0).matrices()[0]) == "x,y,z"
  assert str(t.get(4).matrices()[0]) == "x,y,z"
  tc = t.deep_copy()
  assert list(tc.indices()) == [1,1,2,1,0]
  tc.process(site_symmetry_ops=s2)
  assert list(tc.indices()) == [1,1,2,1,0,2]
  assert list(t.indices()) == [1,1,2,1,0]
  cb_op = sgtbx.change_of_basis_op("y,z,x")
  tcb = t.change_basis(cb_op=cb_op)
  assert list(tcb.indices()) == list(t.indices())
  assert list(tcb.special_position_indices()) \
      == list(t.special_position_indices())
  for tcbt,tt in zip(tcb.table(), t.table()):
    for a,b in zip(tcbt.matrices(), tt.matrices()):
      assert str(a) == str(cb_op.apply(b))
  t.process(site_symmetry_ops=s3)
  assert list(t.indices()) == [1,1,2,1,0,0]
  assert list(tc.indices()) == [1,1,2,1,0,2]
  assert list(t.special_position_indices()) == [0,1,2,3]
  assert list(tc.special_position_indices()) == [0,1,2,3,5]
  ts = tc.select(selection=flex.size_t([5,0,4,1]))
  assert list(ts.indices()) == [1,2,0,2]
  assert ts.get(0).special_op() == tc.get(5).special_op()
  assert ts.get(1).special_op() == tc.get(0).special_op()
  assert ts.get(2).special_op() == tc.get(4).special_op()
  assert ts.get(3).special_op() == tc.get(1).special_op()
  ts = tc.select(selection=flex.bool([True,True,False,False,True,True]))
  assert ts.get(0).special_op() == tc.get(0).special_op()
  assert ts.get(1).special_op() == tc.get(1).special_op()
  assert ts.get(2).special_op() == tc.get(4).special_op()
  assert ts.get(3).special_op() == tc.get(5).special_op()
  ts = tc.select(selection=flex.size_t())
  assert ts.indices().size() == 0
  for i in xrange(tc.indices().size()):
    s = tc.get(i)
    p = pickle.dumps(s)
    l = pickle.loads(p)
    assert l.multiplicity() == s.multiplicity()
    assert l.special_op() == s.special_op()
    if (i == 0):
      assert len(s.matrices()) == 12
      assert len(l.matrices()) == 12
    for unpickled,original in zip(l.matrices(), s.matrices()):
      assert unpickled == original
  p = pickle.dumps(ts)
  l = pickle.loads(p)
  assert l.indices().size() == 0
  assert l.table() == ()
  assert l.special_position_indices().size() == 0
  p = pickle.dumps(tc)
  l = pickle.loads(p)
  assert list(l.indices()) == [1,1,2,1,0,2]
  assert len(l.table()) == 3
  assert str(l.table()[0].special_op()) == "x,y,z"
  for unpickled,original in zip(l.table(), tc.table()):
    assert unpickled.special_op() == original.special_op()
    for unp,orig in zip(unpickled.matrices(), original.matrices()):
      assert unp == orig
  assert list(l.special_position_indices()) == [0,1,2,3,5]
  original_sites_frac=flex.vec3_double([
    (0,0,0),
    (0.25,0.66,0.0),
    (0,0,0),
    (0.5,0.5,0.5)])
  t = sgtbx.site_symmetry_table()
  t.process(
    unit_cell=u,
    space_group=g,
    original_sites_frac=original_sites_frac)
  assert list(t.indices()) == [1,0,1,2]
  #
  d = t.discard_symmetry()
  assert list(d.indices()) == [0,0,0,0]
  assert len(d.table())==1
  assert d.table()[0].is_point_group_1()
  assert d.special_position_indices().size()==0
  #
  t = sgtbx.site_symmetry_table()
  t.process(
    unit_cell=u,
    space_group=g,
    original_sites_frac=original_sites_frac,
    min_distance_sym_equiv=100,
    assert_min_distance_sym_equiv=False)
  assert list(t.indices()) == [1,2,1,3]
  for ugpf,ti in [
        ([False,False,False,False], [1,0,1,2]),
        ([False,False,True,False], [1,0,0,2]),
        ([True,False,False,False], [0,0,1,2]),
        ([True,False,True,False], [0,0,0,1]),
        ([True,False,True,True], [0,0,0,0]),
        ([True,True,True,True], [0,0,0,0])]:
    t = sgtbx.site_symmetry_table()
    t.process(
      unit_cell=u,
      space_group=g,
      original_sites_frac=original_sites_frac,
      unconditional_general_position_flags=flex.bool(ugpf))
    assert list(t.indices()) == ti
  #
  u = uctbx.unit_cell((3,4,5,80,100,110))
  g = sgtbx.space_group("P 2")
  #
  ss = site_symmetry(
    unit_cell=u,
    space_group=g,
    original_site=(0,0,0))
  assert str(ss.special_op()) == "0,0,z"
  ss = site_symmetry(
    unit_cell=u,
    space_group=g,
    original_site=(0,0,0),
    min_distance_sym_equiv=0)
  assert str(ss.special_op()) == "x,y,z"
  #
  ss = site_symmetry(
    unit_cell=u,
    space_group=g,
    original_site=(1.05,-2,0.123))
  assert str(ss.special_op()) == "1,-2,z"
  sc = ss.site_constraints()
  assert sc.row_echelon_lcm == 12
  assert list(sc.row_echelon_form()) == [24, 0, 0, 0, 24, 0]
  assert approx_equal(sc.row_echelon_constants, [24, -48])
  assert list(sc.independent_indices) == [2]
  assert sc.n_independent_params() == 1
  assert sc.n_dependent_params() == 2
  ip = sc.independent_params(all_params=(7,8,0.123))
  assert approx_equal(ip, [0.123])
  ap = sc.all_params(independent_params=ip)
  assert sc.all_shifts([0.128]) == (0, 0, 0.128)
  assert approx_equal(ap, [1.0,-2.0,0.123])
  assert sc.gradient_sum_matrix().focus() == (1,3)
  assert approx_equal(sc.gradient_sum_matrix(), [0,0,1])
  ig = sc.independent_gradients(
    all_gradients=flex.double([0.1,0.2,0.5]))
  assert approx_equal(ig, [0.5])
  ig = sc.independent_gradients((0.1,0.2,0.5))
  assert approx_equal(ig, [0.5])
  ic = sc.independent_curvatures(
    all_curvatures=flex.double([0.1,0.2,0.5,0.3,0.2,-0.1]))
  assert approx_equal(ic, [-0.1])
  #
  u = uctbx.unit_cell((3,4,5,80,100,110))
  g = sgtbx.space_group("P 2")
  s1 = site_symmetry(u, g, (0.1,0.2,0.3))
  s2 = site_symmetry(u, g, (0,0,0))
  s3 = site_symmetry(u, g, (0.5,0.5,0))
  assert str(s1.special_op()) == "x,y,z"
  assert str(s2.special_op()) == "0,0,z"
  assert str(s3.special_op()) == "1/2,1/2,z"
  t = sgtbx.site_symmetry_table()
  t.process(insert_at_index=0, site_symmetry_ops=s1)
  assert list(t.indices()) == [0]
  assert list(t.special_position_indices()) == []
  t.process(insert_at_index=0, site_symmetry_ops=s2)
  assert list(t.indices()) == [1,0]
  assert list(t.special_position_indices()) == [0]
  t = sgtbx.site_symmetry_table()
  t.process(insert_at_index=0, site_symmetry_ops=s2)
  assert list(t.indices()) == [1]
  assert list(t.special_position_indices()) == [0]
  t.process(insert_at_index=0, site_symmetry_ops=s1)
  assert list(t.indices()) == [0,1]
  assert list(t.special_position_indices()) == [1]
  t.process(insert_at_index=1, site_symmetry_ops=s3)
  assert list(t.indices()) == [0,2,1]
  assert list(t.special_position_indices()) == [1,2]
  t.process(insert_at_index=0, site_symmetry_ops=s2)
  assert list(t.indices()) == [1,0,2,1]
  assert list(t.special_position_indices()) == [0,2,3]
  t.process(insert_at_index=4, site_symmetry_ops=s3)
  assert list(t.indices()) == [1,0,2,1,2]
  assert list(t.special_position_indices()) == [0,2,3,4]
  t.process(insert_at_index=0, site_symmetry_ops=s1)
  assert list(t.indices()) == [0,1,0,2,1,2]
  assert list(t.special_position_indices()) == [1,3,4,5]
  t.process(insert_at_index=6, site_symmetry_ops=s1)
  assert list(t.indices()) == [0,1,0,2,1,2,0]
  assert list(t.special_position_indices()) == [1,3,4,5]
  #
  u = uctbx.unit_cell((3,4,5,80,100,110))
  g = sgtbx.space_group("-P")
  site_c = site_symmetry(u, g, (0, 1, 0.5)).site_constraints()
  assert site_c.all_shifts(()) == (0, 0, 0)
  #
  u = uctbx.unit_cell((2, 2, 2, 80, 80, 80))
  g = sgtbx.space_group("R 3 (-y+z, x+z, -x+y+z)")
  site_c = site_symmetry(u, g, (0.1, 0.1, 0.1)).site_constraints()
  assert site_c.all_shifts((0.128,)) == (0.128, 0.128, 0.128)

def exercise_wyckoff():
  space_group_type = sgtbx.space_group_type
  wyckoff_table = sgtbx.wyckoff_table
  assert wyckoff_table(space_group_type("P 1")).size() == 1
  sg_type = space_group_type("P m m m")
  w = wyckoff_table(sg_type)
  assert w.space_group_type().group() == sg_type.group()
  assert w.size() == 27
  assert w.position(0).special_op().is_unit_mx()
  assert w.position("@").special_op().is_unit_mx()
  assert w.position("n").special_op() == w.position(13).special_op()
  assert w.lookup_index("@") == 0
  assert w.lookup_index("k") == 16
  for i in xrange(w.size()):
    assert w.lookup_index(w.position(i).letter()) == i
  p = w.position("t")
  assert p.multiplicity() == 2
  assert p.point_group_type() == "mm2"
  assert [str(s) for s in p.unique_ops(sg_type.group())] \
      == ["1/2,1/2,z", "1/2,1/2,-z"]
  u = uctbx.unit_cell((2,3,4))
  g = sg_type.group()
  x = (0.5,0.,0.123)
  ss = sgtbx.site_symmetry(u, g, x)
  m = w.mapping(ss)
  ww = weakref.ref(w)
  del w # to verify life-time support
  assert u.is_similar_to(m.unit_cell())
  assert approx_equal(m.original_site(), x)
  assert m.position().letter() == "s"
  assert m.position().point_group_type() == "mm2"
  assert str(m.sym_op()) == "x,y,z"
  assert approx_equal(m.representative_site(), x)
  assert approx_equal(m.exact_site(), x)
  assert m.distance_moved() < 1.e-6
  assert str(m.special_op()) == "1/2,0,z"
  assert ww().size() == 27
  del m
  assert ww().size() == 27
  del p
  assert ww() is None
  w = wyckoff_table(sg_type)
  x = (-2.52, 1.123, -4.97)
  ss = sgtbx.site_symmetry(u, g, x)
  m = w.mapping(ss)
  assert u.is_similar_to(m.unit_cell())
  assert approx_equal(m.original_site(), x)
  assert m.position().letter() == "o"
  assert m.position().point_group_type() == "mm2"
  assert str(m.sym_op()) == "x+3,y,z+5"
  assert approx_equal(m.representative_site(), (0.5,1.123,0))
  assert approx_equal(m.exact_site(), (-2.5,1.123,-5))
  assert approx_equal(m.distance_moved(), u.length((-0.02,0,-0.03)))
  assert str(m.special_op()) == "-5/2,y,-5"
  m = w.mapping(u, x)
  assert u.is_similar_to(m.unit_cell())
  assert approx_equal(m.original_site(), x)
  assert m.position().letter() == "b"
  assert m.position().point_group_type() == "mmm"
  xyz = str(m.sym_op())
  assert xyz in ("x+3,y-1,z+5", "x+3,-y+1,z+5"), xyz
  assert approx_equal(m.representative_site(), (0.5,0,0))
  assert approx_equal(m.exact_site(), (-2.5, 1.0, -5.0))
  assert approx_equal(m.distance_moved(), u.length((-0.02,0.123,-0.03)))
  assert str(m.special_op()) == "-5/2,1,-5"
  m = w.mapping(u, x, 0.2)
  assert m.position().letter() == "o"
  assert str(m.sym_op()) in ("x+3,y-2,z+5", "x+3,-y+2,z+5")
  if (str(m.sym_op()) == "x+3,y-2,z+5"):
    assert approx_equal(m.representative_site(), (0.5,1.123-2,0))
  else:
    assert approx_equal(m.representative_site(), (0.5,-1.123+2,0))
  assert approx_equal(m.exact_site(), (-2.5,1.123,-5))
  assert approx_equal(m.distance_moved(), u.length((-0.02,0,-0.03)))
  assert str(m.special_op()) == "-5/2,y,-5"
  #
  sg_type = space_group_type("P 42/n c m")
  w = wyckoff_table(sg_type)
  assert str(w.position(3).special_op()) == "1/2*x-1/2*y,-1/2*x+1/2*y,1/2"
  assert str(w.position(3).special_op_simplified()) == "x,-x,1/2"

def exercise_sym_equiv_sites():
  sym_equiv_sites = sgtbx.sym_equiv_sites
  u = uctbx.unit_cell((8,8,11,90,90,120))
  sg_type = sgtbx.space_group_type("P 3 1 2")
  g = sg_type.group()
  wtab = sgtbx.wyckoff_table(sg_type)
  for x,mult,sym_i in (((0,0,0), 1, (0,)),
                      ((1./5,-1./5,1.5), 3, (0,1,2)),
                      ((-1./3,1./3,0.123), 2, (0,3)),
                      ((1./4,1./5,1./3), 6, (0,1,2,3,4,5))):
    ss = sgtbx.site_symmetry(u, g, x)
    assert ss.multiplicity() == mult
    wm = wtab.mapping(u, x)
    assert wm.position().multiplicity() == mult
    for m in (ss, wm):
      e = sym_equiv_sites(m)
      assert e.unit_cell().is_similar_to(u)
      assert e.space_group() == g
      assert e.original_site() == x
      assert e.special_op() == ss.special_op()
      assert e.max_accepted_tolerance() < 0
      assert e.coordinates().size() == mult
      if (mult == 6):
        c = tuple(e.coordinates())
        vfy = ((1./4, 0.2, 1./3),
               (-1./5, 1./20, 1./3),
               (-1./20, -1./4, 1./3),
               (-1./5, -1./4, -1./3),
               (1./4, 1./20, -1./3),
               (-1./20, 1./5, -1./3))
        for i in xrange(6):
          assert approx_equal(c[i], vfy[i])
      assert tuple(e.sym_op_indices()) == sym_i
      assert e.is_special_position() \
          == (len(e.coordinates()) < e.space_group().order_z())
      for i in xrange(e.coordinates().size()):
        assert e.sym_op(i) == g(sym_i[i])
    for i in xrange(2):
      if (i == 0):
        e = sym_equiv_sites(g, x)
      if (i == 1):
        e = sym_equiv_sites(g, x, u)
        assert e.unit_cell().is_similar_to(u)
      assert e.space_group() == g
      assert e.original_site() == x
      assert not e.special_op().is_valid()
      assert e.max_accepted_tolerance() < 0
      assert e.coordinates().size() == g.order_z()
    e = sym_equiv_sites(u, g, x, ss.special_op())
    assert e.unit_cell().is_similar_to(u)
    assert e.space_group() == g
    assert e.original_site() == x
    assert e.special_op() == ss.special_op()
    assert e.max_accepted_tolerance() < 0
    assert e.coordinates().size() == mult
    for i in xrange(4):
      if (i == 0):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x,
          site_symmetry_ops=ss)
      elif (i == 1):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x)
      elif (i == 2):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x,
          minimum_distance=0.1)
      elif (i == 3):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x,
          minimum_distance=0.1,
          tolerance=0.001)
      assert e.unit_cell().is_similar_to(u)
      assert e.space_group() == g
      assert e.original_site() == x
      if (i == 0):
        assert e.special_op().is_valid()
        assert e.max_accepted_tolerance() < 0
      else:
        assert not e.special_op().is_valid()
        assert e.max_accepted_tolerance() >= 0
      assert e.coordinates().size() == mult
    e = sym_equiv_sites(u, g, x, 100, 50)
    assert e.coordinates().size() == 1
    for i in xrange(6):
      if (i == 0):
        e = sym_equiv_sites(
          site_symmetry=ss)
      elif (i == 1):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x,
          site_symmetry_ops=ss)
      elif (i == 2):
        e = sym_equiv_sites(
          wyckoff_mapping=wm)
      elif (i == 3):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x,
          special_op=ss.special_op())
      elif (i == 4):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x,
          special_op=ss.special_op(),
          multiplicity=ss.multiplicity())
      elif (i == 5):
        e = sym_equiv_sites(
          unit_cell=u,
          space_group=g,
          original_site=x)
      d = sgtbx.min_sym_equiv_distance_info(reference_sites=e, other=x)
      assert d.i_other() == 0
      assert str(d.sym_op()) == "x,y,z"
      assert d.continuous_shifts() == (0,0,0)
      assert approx_equal(d.diff(), (0,0,0))
      assert approx_equal(d.dist(), 0)
      assert approx_equal(d.sym_op() * x, x)
      a = d.apply(sites_frac=e.coordinates())
      assert a.size() == e.coordinates().size()
      for i,y in enumerate(e.coordinates()):
        assert approx_equal(a[i], y)
      for i,y in enumerate(e.coordinates()):
        d = sgtbx.min_sym_equiv_distance_info(e, y)
        assert approx_equal(d.dist(), 0)
        assert approx_equal(d.sym_op() * y, x)
      d = sgtbx.min_sym_equiv_distance_info(
        reference_sites=e,
        others=flex.vec3_double([(x[0]+0.1,x[1]+0.2,x[2]+0.3), x]))
      assert d.i_other() == 1
      assert approx_equal(d.continuous_shifts(), (0,0,0))
      assert approx_equal(d.dist(), 0)
  x = (1.732, -1.414, 2.236)
  for hall_symbol in ("P 1",
                      "P -1 (2,-1,1)",
                      "R 1 (3,2,-1)",
                      "P -4 2 3 (2,1,-1)",
                      "-P 2ab 2bc 3 (1,2,-2)",
                      "F -4 2 3 (-2,2,-1)",
                      "-F 2uv 2vw 3 (1,-2,-1)"):
    g = sgtbx.space_group(hall_symbol)
    e = sym_equiv_sites(g, x)
    c = e.coordinates()
    assert c.size() == g.order_z()
    for i,v in enumerate(c):
      assert approx_equal(g(i) * x, v)
  g = sgtbx.space_group()
  e = sym_equiv_sites(g, (0.016, 0.895, 0.111), uctbx.unit_cell(()))
  d = sgtbx.min_sym_equiv_distance_info(e, (0.939, 0.128, 0.178))
  assert approx_equal(d.dist()**2, 0.064707)
  for hall_symbol,shift_flags in (("P 1", (True,True,True)),
                                  ("P 2", (False,False,True)),
                                  ("P -2", (True,True,False)),
                                  ("P 3x", (True,False,False)),
                                  ("P 3 2", (False,False,False))):
    for x in ((1.732, -1.414, 2.236), (0.939, 0.128, 0.178)):
      g = sgtbx.space_group(hall_symbol)
      e = sym_equiv_sites(g, x, uctbx.unit_cell(()))
      d = sgtbx.min_sym_equiv_distance_info(
        reference_sites=e,
        other=x,
        principal_continuous_allowed_origin_shift_flags=shift_flags)
      assert approx_equal(d.continuous_shifts(), (0,0,0))
      assert approx_equal(d.dist(), 0)
      shift = [s * f for s,f in zip(shift_flags, (0.123,0.234,0.345))]
      for y in e.coordinates():
        d = sgtbx.min_sym_equiv_distance_info(e, y, shift_flags)
        assert approx_equal(d.dist(), 0)
        assert approx_equal(d.sym_op() * y, x)
        z = [y[i]+shift[i] for i in (0,1,2)]
        d = sgtbx.min_sym_equiv_distance_info(e, z, shift_flags)
        assert approx_equal(d.dist(), 0)
        if (shift_flags != (0,0,0)):
          assert not_approx_equal(d.sym_op() * z, x)
        else:
          assert approx_equal(d.sym_op() * z, x)
        fz = flex.vec3_double(1, z)
        assert approx_equal(d.apply(fz)[0], x)
  g = sgtbx.space_group()
  e = sym_equiv_sites(g, (0.1,0.2,0.3), u)
  d = sgtbx.min_sym_equiv_distance_info(e, (0.2,0.4,0.1))
  assert approx_equal(d.diff(), (-0.1,-0.2,0.2))
  assert approx_equal(d.dist(), u.length((-0.1,-0.2,0.2)))
  #
  u = uctbx.unit_cell((3,4,5,80,100,110))
  g = sgtbx.space_group("P 2")
  e = sym_equiv_sites(unit_cell=u, space_group=g, original_site=(0,0,0))
  assert e.coordinates().size() == 1
  e = sym_equiv_sites(
    unit_cell=u,
    space_group=g,
    original_site=(0,0,0),
    minimum_distance=0,
    tolerance=0)
  assert e.coordinates().size() == 2

def exercise_seminvariant():
  space_group = sgtbx.space_group
  structure_seminvariants = sgtbx.structure_seminvariants
  tests = (
    ("P 1",  [((1, 0, 0), 0), ((0, 1, 0), 0), ((0, 0, 1), 0)]),
    ("-P 1", [((1, 0, 0), 2), ((0, 1, 0), 2), ((0, 0, 1), 2)]),

    ("P 2",   [((1, 0, 0), 2), ((0, 1, 0), 2), ((0, 0, 1), 0)]),
    ("P 2y",  [((1, 0, 0), 2), ((0, 1, 0), 0), ((0, 0, 1), 2)]),
    ("P 2x",  [((1, 0, 0), 0), ((0, 1, 0), 2), ((0, 0, 1), 2)]),
    ("P -2",  [((1, 0, 0), 0), ((0, 1, 0), 0), ((0, 0, 1), 2)]),
    ("P -2y", [((1, 0, 0), 0), ((0, 1, 0), 2), ((0, 0, 1), 0)]),
    ("P -2x", [((1, 0, 0), 2), ((0, 1, 0), 0), ((0, 0, 1), 0)]),

    ("P 2 2",   [((1, 0, 0), 2), ((0, 1, 0), 2), ((0, 0, 1), 2)]),
    ("P -2 2",  [((1, 0, 0), 0), ((0, 1, 0), 2), ((0, 0, 1), 2)]),
    ("P -2 -2", [((1, 0, 0), 2), ((0, 1, 0), 0), ((0, 0, 1), 2)]),
    ("P 2 -2",  [((1, 0, 0), 2), ((0, 1, 0), 2), ((0, 0, 1), 0)]),
    ("-P 2 2",  [((1, 0, 0), 2), ((0, 1, 0), 2), ((0, 0, 1), 2)]),

    ("P 3",   [((0, 0, 1), 0), ((1, 2, 0), 3)]),
    ("P 3y",  [((0, 1, 0), 0), ((1, 0, 2), 3)]),
    ("P 3x",  [((1, 0, 0), 0), ((0, 1, 2), 3)]),
    ("P 3*",  [((1, 1, 1), 0)]),
    ("-P 3",  [((0, 0, 1), 2)]),
    ("-P 3y", [((0, 1, 0), 2)]),
    ("-P 3x", [((1, 0, 0), 2)]),
    ("-P 3*", [((1, 1, 1), 2)]),

    ("P 3 2",  [((2, 4, 3), 6)]),
    ("P 3y 2", [((2, 3, 4), 6)]),
    ("P 3x 2", [((3, 2, 4), 6)]),

    ('P 3 2"',  [((0, 0, 1), 2)]),
    ('P 3y 2"', [((0, 1, 0), 2)]),
    ('P 3x 2"', [((1, 0, 0), 2)]),

    ("P 4",   [((0, 0, 1), 0), ((1, 1, 0), 2)]),
    ("P 4y",  [((0, 1, 0), 0), ((1, 0, 1), 2)]),
    ("P 4x",  [((1, 0, 0), 0), ((0, 1, 1), 2)]),
    ("-P 4",  [((0, 0, 1), 2), ((1, 1, 0), 2)]),
    ("-P 4y", [((0, 1, 0), 2), ((1, 0, 1), 2)]),
    ("-P 4x", [((1, 0, 0), 2), ((0, 1, 1), 2)]),

    ("P 6",  [((0, 0, 1), 0)]),
    ("P 6y", [((0, 1, 0), 0)]),
    ("P 6x", [((1, 0, 0), 0)]),

    ("-C 2 2",   [((1, 0, 0), 2), ((0, 0, 1), 2)]),
    ("-I 2 2",   [((1, 0, 0), 2), ((0, 1, 0), 2)]),
    ("-F 2 2",   [((1, 0, 0), 2)]),
    ("-I 4",     [((0, 0, 1), 2)]),
    ("-I 2 2 3", []),

    ("C 2y",    [((0, 1, 0), 0), ((0, 0, 1), 2)]),
    ("C -2y",   [((1, 0, 0), 0), ((0, 0, 1), 0)]),
    ("C 2 -2",  [((1, 0, 0), 2), ((0, 0, 1), 0)]),
    ("C 2 2",   [((1, 0, 0), 2), ((0, 0, 1), 2)]),
    ("A 2 -2",  [((1, 0, 0), 2), ((0, 0, 1), 0)]),
    ("I 2 -2",  [((1, 0, 0), 2), ((0, 0, 1), 0)]),
    ("I 2 2",   [((1, 0, 0), 2), ((0, 1, 0), 2)]),
    ("F 4 2 3", [((1, 0, 0), 2)]),
    ("F 2 2",   [((1, 1, 1), 4)]),
    ("I 4",     [((0, 0, 1), 0)]),
    ("I 4 2",   [((0, 0, 1), 2)]),
    ("I -4",    [((2, 0, 1), 4)]),
    ("F 2 -2",  [((0, 0, 1), 0)]),
    ("I 2 3",   []),
  )
  for hs,vfy in tests:
    g = space_group(hs)
    ss = structure_seminvariants(space_group=g)
    assert ss.size() == len(ss.vectors_and_moduli())
    assert [(vm.v, vm.m) for vm in ss.vectors_and_moduli()] == vfy
  ss = structure_seminvariants(space_group("P 2 -2"))
  assert ss.is_ss(miller_index=(2,-2,0))
  assert ss.apply_mod(miller_index=(2,-2,0)) == (0,0,0)
  assert not ss.is_ss((2,1,0))
  assert ss.apply_mod((2,1,0)) == (0,1,0)
  ss = structure_seminvariants(space_group("C 2y"))
  assert ss.is_ss((2,0,4))
  assert ss.apply_mod((2,0,4)) == (0,0)
  assert not ss.is_ss((2,1,4))
  assert ss.apply_mod((2,1,4)) == (1,0)
  ss = structure_seminvariants(space_group("C 2 2"))
  assert ss.is_ss((4,2,6))
  assert ss.apply_mod((4,2,6)) == (0,0)
  assert not ss.is_ss((3,2,4))
  assert ss.apply_mod((3,2,4)) == (1,0)
  assert ss.apply_mod((2,3,5)) == (0,1)
  assert ss.apply_mod((3,3,5)) == (1,1)
  ss = structure_seminvariants(space_group("I 4"))
  assert ss.is_ss((1,3,0))
  assert ss.apply_mod((1,3,0)) == (0,)
  assert not ss.is_ss((1,3,4))
  assert ss.apply_mod((1,3,4)) == (4,)
  ss = structure_seminvariants(space_group("I 2 3"))
  assert ss.is_ss((1,3,4))
  assert ss.apply_mod((1,3,4)) == ()
  ss = structure_seminvariants(space_group("I -4"))
  assert ss.gridding() == (2,1,4)
  assert ss.refine_gridding(grid=(3,2,10)) == (6,2,20)
  ss = structure_seminvariants(space_group("P -2x"))
  a = ss.grid_adapted_moduli(dim=(3,5,7))
  assert [vm.m for vm in a] == [2,5,7]
  ss = structure_seminvariants(space_group("P 3*"))
  a = ss.grid_adapted_moduli((3,5,7))
  assert [vm.m for vm in a] == [3*5*7]
  a = ss.grid_adapted_moduli((3,5,12))
  assert [vm.m for vm in a] == [60]
  ss = structure_seminvariants(space_group("P 2 -2"))
  ss_continuous = ss.select(discrete=False)
  assert ss_continuous.size() == 1
  assert ss_continuous.vectors_and_moduli()[0].m == 0
  assert ss_continuous.vectors_and_moduli()[0].v == (0,0,1)
  ss_discrete = ss.select(discrete=True)
  assert ss_discrete.size() == 2
  assert ss_discrete.vectors_and_moduli()[0].m == 2
  assert ss_discrete.vectors_and_moduli()[0].v == (1,0,0)
  assert ss_discrete.vectors_and_moduli()[1].m == 2
  assert ss_discrete.vectors_and_moduli()[1].v == (0,1,0)
  assert ss.continuous_shifts_are_principal()
  assert ss.principal_continuous_shift_flags() == (False,False,True)
  assert approx_equal(
    ss.subtract_principal_continuous_shifts(translation=[1,2,3]), [1,2,0])
  ss = structure_seminvariants(space_group("P 3*"))
  assert not ss.continuous_shifts_are_principal()
  assert ss.principal_continuous_shift_flags(assert_principal=False) \
      == (True,True,True)
  assert approx_equal(
    ss.subtract_principal_continuous_shifts(
      translation=[1,2,3], assert_principal=False), [0,0,0])

def exercise_lattice_symmetry():
  reduced_cell = uctbx.unit_cell((12,13,14,80,83,86))
  lattice_group = sgtbx.lattice_symmetry_group(
    reduced_cell=reduced_cell,
    max_delta=3,
    enforce_max_delta_for_generated_two_folds=False)
  assert lattice_group.order_z() == 1
  assert sgtbx.lattice_symmetry_find_max_delta(
    reduced_cell=reduced_cell, space_group=lattice_group) < 1.e-10

def exercise_lattice_symmetry_options():
  niggli_cell = uctbx.unit_cell((405.08382,411.53719,417.48903,
                                  89.97113,89.77758,89.90342))
  testing_combinations=[(0.9,True,10),(0.9,False,10),
                        (1.4,True,10),(1.4,False,30),
                        (1.8,True,30),(1.8,False,30)]
  for max_delta,enforce,n_subgrps in testing_combinations:
    group = sgtbx.lattice_symmetry_group(
      reduced_cell=niggli_cell,
      max_delta=max_delta,
      enforce_max_delta_for_generated_two_folds=enforce)
    group_info = sgtbx.space_group_info(group=group)
    subgrs = subgroups.subgroups(group_info).groups_parent_setting()
    assert len(subgrs) == n_subgrps

def exercise_find_affine():
  for sym,n in [("P 1",67704),("P 2",104),("C 2y",36),("B 2x",36),("A 2z",36)]:
    group = sgtbx.space_group(sym)
    for use_p1_algorithm in [False, True]:
      affine = sgtbx.find_affine(group, 2)
      cb_mx = affine.cb_mx()
      assert len(cb_mx) == n
      if (group.n_smx() == 1): break

def exercise_search_symmetry():
  f = sgtbx.search_symmetry_flags(use_space_group_symmetry=True)
  f = sgtbx.search_symmetry_flags(
    use_space_group_symmetry=False,
    use_space_group_ltr=1,
    use_seminvariants=False,
    use_normalizer_k2l=True,
    use_normalizer_l2n=False)
  assert not f.use_space_group_symmetry()
  assert f.use_space_group_ltr() > 0
  assert not f.use_seminvariants()
  assert f.use_normalizer_k2l()
  assert not f.use_normalizer_l2n()
  assert f == f
  assert not f != f
  assert f == sgtbx.search_symmetry_flags(False, 1, False, True, False)
  assert not f != sgtbx.search_symmetry_flags(False, 1, False, True, False)
  assert f != sgtbx.search_symmetry_flags(True, 0, True, False, True)
  assert not f == sgtbx.search_symmetry_flags(False, 0, True, False, True)
  for use_space_group_symmetry in [False, True]:
    for use_space_group_ltr in [-1,0,1]:
      for use_seminvariants in [False, True]:
        for use_normalizer_k2l in [False, True]:
          for use_normalizer_l2n in [False, True]:
            f = sgtbx.search_symmetry_flags(
              use_space_group_symmetry,
              use_space_group_ltr,
              use_seminvariants,
              use_normalizer_k2l,
              use_normalizer_l2n)
            p = pickle.dumps(f)
            l = pickle.loads(p)
            assert l == f
  sg143 = sgtbx.space_group_info("P 3")
  ss143 = sg143.structure_seminvariants()
  sg146 = sgtbx.space_group_info("R 3 h")
  ss146 = sg146.structure_seminvariants()
  sg149 = sgtbx.space_group_info("P 3 1 2")
  ss149 = sg149.structure_seminvariants()
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, False, False, False),
    space_group_type=sg149.type())
  assert s.subgroup().type().lookup_symbol() == "P 1"
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, False, False, False),
    space_group_type=sg149.type())
  assert s.subgroup().type().lookup_symbol() == "P 3 1 2"
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 1, False, False, False),
    space_group_type=sg149.type())
  assert s.subgroup().type().lookup_symbol() == "P 1"
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, True, False, False),
    space_group_type=sg149.type(),
    seminvariant=ss149)
  assert s.subgroup() \
      == sgtbx.space_group("P 1 (-1/3*y-1/3*z,1/3*y-2/3*z,1/2*x)")
  assert s.continuous_shifts() == ()
  assert s.continuous_shifts_are_principal()
  assert s.continuous_shift_flags(assert_principal=True) == (False,False,False)
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, False, True, False),
    space_group_type=sg149.type())
  assert s.subgroup() == sgtbx.space_group("-P 1")
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, False, False, True),
    space_group_type=sg149.type())
  assert s.subgroup() == sgtbx.space_group("P 2z")
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, False, True, False),
    space_group_type=sg149.type())
  assert s.subgroup().type().lookup_symbol() == "P -3 1 m"
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, False, True, True),
    space_group_type=sg149.type())
  assert s.flags() == sgtbx.search_symmetry_flags(True, 0, False, True, True)
  assert s.subgroup().type().lookup_symbol() == "P 6/m m m"
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, True, True, True),
    space_group_type=sg149.type(),
    seminvariant=ss149)
  assert s.subgroup() \
      == sgtbx.space_group("-P 6 2 (1/3*x+1/3*y,-1/3*x+2/3*y,1/2*z)")
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 1, False, False, False),
    space_group_type=sg146.type())
  assert s.subgroup() == sgtbx.space_group("R 1")
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, False, True, False),
    space_group_type=sg146.type())
  assert s.subgroup() == sgtbx.space_group("R -3")
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, False, True, True),
    space_group_type=sg146.type())
  assert s.subgroup() == sgtbx.space_group('-R 3 2"')
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, False, True, True),
    space_group_type=sg146.type())
  assert s.subgroup() == sgtbx.space_group('-R 3 2"')
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, False, False, False),
    space_group_type=sg146.type(),
    seminvariant=ss146)
  assert s.subgroup().type().lookup_symbol() == "P 1"
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, True, False, False),
    space_group_type=sg146.type(),
    seminvariant=ss146)
  assert s.subgroup() == sgtbx.space_group("R 1")
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, -1, True, False, False),
    space_group_type=sg146.type(),
    seminvariant=ss146)
  assert s.subgroup() == sgtbx.space_group("P 1")
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, True, False, False),
    space_group_type=sg143.type(),
    seminvariant=ss143)
  assert s.subgroup().type().hall_symbol() \
      == " P 1 (2/3*x-1/3*y,1/3*x+1/3*y,z)"
  assert s.continuous_shifts() == ((0,0,1),)
  assert s.continuous_shifts_are_principal()
  assert s.continuous_shift_flags() == (False,False,True)
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, True, True, True),
    space_group_type=sg143.type(),
    seminvariant=ss143)
  assert s.subgroup().type().hall_symbol() \
      == "-P 6 2 (1/3*x+1/3*y,-1/3*x+2/3*y,z)"
  assert s.continuous_shifts() == ((0,0,1),)
  sg1 = sgtbx.space_group_info(number=1)
  ss1 = sg1.structure_seminvariants()
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, True, True, True),
    space_group_type=sg1.type(),
    seminvariant=ss1)
  assert s.subgroup() == sgtbx.space_group("-P 1")
  assert s.continuous_shifts() == ((1, 0, 0), (0, 1, 0), (0, 0, 1))
  assert s.continuous_shifts_are_principal()
  assert s.continuous_shift_flags() == (True,True,True)
  sg6 = sgtbx.space_group_info(number=6)
  ss6 = sg6.structure_seminvariants()
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, True, True, True),
    space_group_type=sg6.type(),
    seminvariant=ss6)
  assert s.subgroup() == sgtbx.space_group("-P 2y (x,1/2*y,z)")
  assert s.continuous_shifts() == ((1, 0, 0), (0, 0, 1))
  assert s.continuous_shift_flags() == (True,False,True)
  sg146r = sgtbx.space_group_info("R 3 r")
  ss146r = sg146r.structure_seminvariants()
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(False, 0, True, False, False),
    space_group_type=sg146r.type(),
    seminvariant=ss146r)
  assert s.subgroup() == sgtbx.space_group("P 1")
  assert s.continuous_shifts() == ((1, 1, 1),)
  assert not s.continuous_shifts_are_principal()
  assert s.continuous_shift_flags(assert_principal=False) == (True,True,True)
  sg144 = sgtbx.space_group_info("P 31")
  ss144 = sg144.structure_seminvariants()
  s = sgtbx.search_symmetry(
    flags=sgtbx.search_symmetry_flags(True, 0, True, False, False),
    space_group_type=sg144.type(),
    seminvariant=ss144)
  assert s.subgroup() == sgtbx.space_group("P 31 (1/3*x+1/3*y,-1/3*x+2/3*y,z)")
  assert s.projected_subgroup() \
      == sgtbx.space_group("P 3 (1/3*x+1/3*y,-1/3*x+2/3*y,z)")

def exercise_tensor_rank_2_constraints():
  u = uctbx.unit_cell((12,12,15,90,90,120))
  g = sgtbx.space_group("P 6 2")
  ss = sgtbx.site_symmetry(u, g, (0,0,0))
  tr2c = sgtbx.tensor_rank_2_constraints
  for c in [tr2c(
              space_group=g,
              reciprocal_space=False),
            tr2c(
              symmetry_matrices=ss.matrices(),
              i_first_matrix_to_use=1,
              reciprocal_space=False)]:
    assert list(c.row_echelon_form()) \
        == [1,-1,0,0,0,0,0,1,0,2,0,0,0,0,0,0,1,1,0,0,0,0,0,1]
  for c in [tr2c(
              space_group=g,
              reciprocal_space=True),
            tr2c(
              symmetry_matrices=ss.matrices(),
              i_first_matrix_to_use=1,
              reciprocal_space=True),
            ss.adp_constraints()]:
    assert list(c.row_echelon_form()) \
        == [1,-1,0,0,0,0,0,1,0,-2,0,0,0,0,0,0,1,-1,0,0,0,0,0,1]
    assert list(c.independent_indices) == [2,3]
    assert approx_equal(c.gradient_sum_matrix(), [0,0,1,0,0,0, 2,2,0,1,0,0])
    assert c.n_independent_params() == 2
    assert c.n_dependent_params() == 4
    assert approx_equal(c.independent_params(all_params=[1,2,3,4,5,6]), [3,4])
    assert approx_equal(c.all_params(independent_params=[1,2]), [4,4,1,2,0,0])
    a = [1,2,3,4,5,6]
    s = [3+1/3.,3+1/3.,3,-3-1/3.,0,0]
    assert approx_equal(c.independent_gradients(all_gradients=a), [3,10])
    assert approx_equal(c.independent_gradients(all_gradients=s), [3,10])
    a = flex.double(xrange(1,22))
    assert approx_equal(c.independent_curvatures(all_curvatures=a),[12,35,116])

def exercise_hashing():
  # rt_mx
  sg = sgtbx.space_group('P 4 2 3')
  op_set = set([ op for op in sg ])
  assert len(op_set) == len(sg)
  for g in sg:
    assert g in op_set

  # rot_mx, tr_vec
  d = {}
  for op in sg:
    d.setdefault(op.r(), set())
    d[op.r()].add(op.t())
  for g in sg:
    assert g.r() in d
    assert g.t() in d[g.r()]

  # space_group
  symbol_set = set([ symbol.hall()
                     for symbol in sgtbx.space_group_symbol_iterator() ])
  assert len(symbol_set) == 527
  sg_list = [ sgtbx.space_group(symbol).make_tidy() for symbol in symbol_set ]
  sg_set = set(sg_list)
  assert len(sg_set) == len(sg_list)
  for sg in sg_list:
    assert sg in sg_set

  sg = sgtbx.space_group('P 2z')
  sg.expand_smx(sgtbx.rt_mx('-x, y, -z'))
  try:
    hash(sg)
    raise Exception_expected
  except RuntimeError, e:
    assert str(e).find("tidy") != -1

def exercise_fractional_mod():
  assert approx_equal(
    sgtbx.fractional_mod_positive((-0.1, 0.2, -0.3)),
    (0.9, 0.2, 0.7))
  assert approx_equal(
    sgtbx.fractional_mod_short((0.4, 0.9, 0.7)),
    (0.4, -0.1, -0.3))

def run():
  exercise_symbols()
  exercise_tr_vec()
  exercise_rot_mx()
  exercise_rt_mx()
  exercise_change_of_basis_op()
  exercise_space_group()
  exercise_space_group_type()
  exercise_phase_info()
  exercise_reciprocal_space_asu()
  exercise_brick()
  exercise_site_symmetry()
  exercise_wyckoff()
  exercise_sym_equiv_sites()
  exercise_seminvariant()
  exercise_lattice_symmetry()
  exercise_lattice_symmetry_options()
  exercise_find_affine()
  exercise_search_symmetry()
  exercise_tensor_rank_2_constraints()
  exercise_hashing()
  exercise_fractional_mod()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
