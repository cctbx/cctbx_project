import math
import weakref
import pickle
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx import uctbx
from scitbx.python_utils import complex_math
from scitbx.test_utils import approx_equal

def exercise_symbols():
  s = sgtbx.space_group_symbols("p 2")
  assert s.number() == 3
  assert s.schoenflies() == "C2^1"
  assert s.qualifier() == "b"
  assert s.hermann_mauguin() == "P 1 2 1"
  assert s.extension() == "\0"
  assert s.extended_hermann_mauguin() == "P 1 2 1"
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
  assert i.next().extended_hermann_mauguin() == "P 1"
  assert i.next().extended_hermann_mauguin() == "P -1"
  assert len(tuple(i)) == 528
  i = sgtbx.space_group_symbol_iterator()
  n = 0
  for s in i:
    n += 1
    hm = s.extended_hermann_mauguin()
    assert sgtbx.space_group_symbols(hm).extended_hermann_mauguin() == hm
  assert n == 530

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
  assert float(t) == t.as_double()
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
  assert a.as_string(0001) == "-.25,0,.5"
  assert a.as_string(00000, ";") == "-1/4;0;1/2"

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
  assert i.sense() == 0
  i = rot_mx().info()
  assert i.type() == 1
  assert i.ev() == (0,0,0)
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,0,0,0,-1,0,0,0,1)))
  assert i.type() == 2
  assert i.ev() == (0,0,1), i.ev()
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,0,0,0,1,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (0,1,0), i.ev()
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((1,0,0,0,-1,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (1,0,0), i.ev()
  i = rot_mx_info(rot_mx((1,0,0,0,1,0,0,0,-1)))
  assert i.type() == -2
  assert i.ev() == (0,0,1), i.ev()
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((1,0,0,0,-1,0,0,0,1)))
  assert i.type() == -2
  assert i.ev() == (0,1,0), i.ev()
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,0,0,0,1,0,0,0,1)))
  assert i.type() == -2
  assert i.ev() == (1,0,0), i.ev()
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((0,-1,0,1,-1,0,0,0,1)))
  assert i.type() == 3
  assert i.ev() == (0,0,1), i.ev()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,-1,0,1,-1,0,0,0,1)).inverse())
  assert i.type() == 3
  assert i.ev() == (0,0,1), i.ev()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,1,0,-1,1,0,0,0,-1)))
  assert i.type() == -3
  assert i.ev() == (0,0,1), i.ev()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,1,0,-1,1,0,0,0,-1)).inverse())
  assert i.type() == -3
  assert i.ev() == (0,0,1), i.ev()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,0,1,0,1,0,-1,0,0)))
  assert i.type() == 4
  assert i.ev() == (0,1,0), i.ev()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,0,1,0,1,0,-1,0,0)).inverse())
  assert i.type() == 4
  assert i.ev() == (0,1,0), i.ev()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,0,-1,0,-1,0,1,0,0)))
  assert i.type() == -4
  assert i.ev() == (0,1,0), i.ev()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,0,-1,0,-1,0,1,0,0)).inverse())
  assert i.type() == -4
  assert i.ev() == (0,1,0), i.ev()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((1,0,0,0,1,-1,0,1,0)))
  assert i.type() == 6
  assert i.ev() == (1,0,0), i.ev()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((1,0,0,0,1,-1,0,1,0)).inverse())
  assert i.type() == 6
  assert i.ev() == (1,0,0), i.ev()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((-1,0,0,0,-1,1,0,-1,0)))
  assert i.type() == -6
  assert i.ev() == (1,0,0), i.ev()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((-1,0,0,0,-1,1,0,-1,0)).inverse())
  assert i.type() == -6
  assert i.ev() == (1,0,0), i.ev()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,0,1,1,0,0,0,1,0)))
  assert i.type() == 3
  assert i.ev() == (1,1,1), i.ev()
  assert i.sense() == 1
  i = rot_mx_info(rot_mx((0,0,1,1,0,0,0,1,0)).inverse())
  assert i.type() == 3
  assert i.ev() == (1,1,1), i.ev()
  assert i.sense() == -1
  i = rot_mx_info(rot_mx((0,-1,0,-1,0,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (-1,1,0), i.ev()
  assert i.sense() == 0
  i = rot_mx_info(rot_mx((-1,1,0,0,1,0,0,0,-1)))
  assert i.type() == 2
  assert i.ev() == (1,2,0), i.ev()

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
  s = rt_mx("y,y-x,z+1/4", "", 3, 8)
  assert str(s) == "y,-x+y,z+1/4"
  assert s.as_xyz() == str(s)
  assert s.as_xyz(0) == str(s)
  assert s.as_xyz(1) == "y,-x+y,z+.25"
  assert s.as_xyz(0, 0) == str(s)
  assert s.as_xyz(0, 1) == "y,-x+y,1/4+z"
  assert s.as_xyz(0, 0, "xyz") == str(s)
  assert s.as_xyz(0, 0, "abc") == "b,-a+b,c+1/4"
  assert s.as_xyz(0, 0, "xyz", ",") == str(s)
  assert s.as_xyz(0, 0, "xyz", " ; ") == "y ; -x+y ; z+1/4"
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
  r = (-1,1,0,0,1,0,0,0,-1)
  t = (1/12.,2/12.,3/12.)
  assert str(rt_mx(r, t)) == "-x+y+1/12,y+1/6,-z+1/4"
  assert str(rt_mx(r, t, 2)) == "-x+y+1/12,y+1/6,-z+1/4"
  assert str(rt_mx(r, t, 2, 24)) == "-x+y+1/12,y+1/6,-z+1/4"
  assert str(rt_mx("+0*x-1*y+0*z,+0*x+0*y+1*z,-1*x+0*y-1*z")) == "-y,z,-x-z"

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
  assert c.select(0) == c.c()
  assert c.select(1) == c.c_inv()
  d = c.inverse()
  assert c.c() == d.c_inv()
  assert d.c() == c.c_inv()
  c = change_of_basis_op("x+3/2,y-5/4,z+1/3")
  c.mod_positive_in_place()
  assert str(c.c()) == "x+1/2,y+3/4,z+1/3"
  c.mod_short_in_place()
  assert str(c.c()) == "x+1/2,y-1/4,z+1/3"
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
  i = flex.miller_index(((1,2,3),(2,3,4)))
  assert tuple(c.apply(i)) == ((3,1,2),(4,2,3))
  s = pickle.dumps(c)
  l = pickle.loads(s)
  assert str(c.c()) == str(l.c())

def exercise_space_group():
  tr_vec = sgtbx.tr_vec
  rt_mx = sgtbx.rt_mx
  space_group = sgtbx.space_group
  g = space_group()
  p = sgtbx.parse_string("P 1")
  g = space_group(p)
  p = sgtbx.parse_string("P 1")
  g = space_group(p, 0)
  p = sgtbx.parse_string("P 1")
  g = space_group(p, 0, 0)
  p = sgtbx.parse_string("P 1")
  g = space_group(p, 0, 0, 0)
  assert g.t_den() == sgtbx.sg_t_den
  p = sgtbx.parse_string("P 1")
  g = space_group(p, 0, 0, 0, 2*sgtbx.sg_t_den)
  assert g.t_den() == 2*sgtbx.sg_t_den
  g = space_group("P 1")
  g = space_group("P 1", 1)
  g = space_group("1", 0, 1)
  g = space_group("1", 0, 1, 1)
  g = space_group("1", 0, 1, 1, 3*sgtbx.sg_t_den)
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
  g.expand_ltr(tr_vec())
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
  g.expand_inv(tr_vec((0,0,0)))
  assert g.f_inv() == 2
  assert g.is_centric()
  assert g.is_origin_centric()
  g.reset()
  g.expand_inv(tr_vec((1,2,3)))
  assert g.is_centric()
  assert not g.is_origin_centric()
  g.reset()
  assert g.n_smx() == 1
  g.expand_smx(rt_mx("-x,-y,z"))
  assert g.n_smx() == 2
  for z,n in (("P",1), ("a",2), ("B",2),("c",2), ("I",2), ("r",3), ("F",4)):
    g.reset()
    g.expand_conventional_centring_type(z)
    assert g.n_ltr() == n
  p = sgtbx.parse_string("P 4")
  g.reset()
  g.parse_hall_symbol(p)
  assert g.order_z() == 4
  p = sgtbx.parse_string("P 2x")
  g.parse_hall_symbol(p, 1)
  assert g.order_z() == 8
  p = sgtbx.parse_string("-1")
  g.parse_hall_symbol(p, 0, 1)
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
  assert tuple(g.is_sys_absent(m)) == (0001,00000,00000)
  assert g.is_centric((1,0,0))
  assert not g.is_centric((0,0,4))
  assert tuple(g.is_centric(m)) == (00000,00000,0001)
  p = g.phase_restriction((1,0,0))
  assert not p.sys_abs_was_tested()
  assert p.ht() == 2
  assert g.is_valid_phase((1,0,0), p.ht_angle())
  assert g.is_valid_phase((1,0,0), p.ht_angle(), 0)
  assert g.is_valid_phase((1,0,0), p.ht_angle(1), 1)
  assert not g.is_valid_phase((1,0,0), p.ht_angle()+math.pi/180, 0)
  assert g.is_valid_phase((1,0,0), p.ht_angle()+math.pi/180, 0, 1.e6)
  assert g.multiplicity((1,2,3), 0) == 8
  assert g.multiplicity((1,2,3), 1) == 4
  assert tuple(g.multiplicity(m, 0)) == (2,2,4)
  assert tuple(g.multiplicity(m, 1)) == (1,1,4)
  assert g.epsilon((1,2,3)) == 1
  assert tuple(g.epsilon(m)) == (4,4,1)
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
  g = space_group("C 2 -2c")
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
  g = space_group("P 31")
  assert g.gridding() == (1,1,3)
  assert g.refine_gridding((2,1,4)) == (2,2,12)
  mod_p = ["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3"]
  mod_s = ["x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z-1/3"]
  assert [str(s) for s in g.all_ops()] == mod_p
  assert [str(s) for s in g.all_ops(0)] == mod_p
  assert [str(s) for s in g.all_ops(1)] == mod_p
  assert [str(s) for s in g.all_ops(-1)] == mod_s
  assert [str(s) for s in g.all_ops(-1, 00000)] == mod_s
  assert [str(s) for s in g.all_ops(-1, 0001)] == mod_s
  assert g.all_ops(0, 00000)[1].t().den() == sgtbx.sg_t_den
  assert g.all_ops(0, 0001)[1].t().den() == 3
  g = space_group("-P 4 2")
  g = space_group("P 31")
  assert g.type().number() == 144
  p = pickle.dumps(g)
  l = pickle.loads(p)
  assert g == l

def exercise_space_group_type():
  space_group = sgtbx.space_group
  space_group_type = sgtbx.space_group_type
  t = space_group_type()
  assert t.group() == space_group()
  assert t.number() == 1
  assert t.cb_op().is_identity_op()
  t = space_group_type("P 2")
  assert t.hall_symbol() == " P 2y"
  t = space_group_type("P 2", "a1983")
  assert t.hall_symbol() == " P 2y"
  t = space_group_type("P 2", "i1952")
  assert t.hall_symbol() == " P 2y (z,x,y)"
  t = space_group_type("P 3")
  assert len(t.addl_generators_of_euclidean_normalizer(0, 0)) == 0
  assert len(t.addl_generators_of_euclidean_normalizer(1, 0)) == 1
  assert len(t.addl_generators_of_euclidean_normalizer(0, 1)) == 2
  assert len(t.addl_generators_of_euclidean_normalizer(1, 1)) == 3
  g = t.expand_addl_generators_of_euclidean_normalizer(0, 0)
  assert g == t.group()
  j = t.expand_addl_generators_of_euclidean_normalizer(1, 0).type()
  assert j.lookup_symbol() == "P -3"
  j = t.expand_addl_generators_of_euclidean_normalizer(0, 1).type()
  assert j.lookup_symbol() == "P 6 m m"
  j = t.expand_addl_generators_of_euclidean_normalizer(1, 1).type()
  assert j.lookup_symbol() == "P 6/m m m"
  assert not space_group_type("P 3").is_enantiomorphic()
  assert space_group_type("P 31").is_enantiomorphic()
  assert space_group_type("P -1").change_of_hand_op().is_identity_op()
  assert str(space_group_type("P 1").change_of_hand_op().c()) == "-x,-y,-z"
  assert str(space_group_type("P 31").change_of_hand_op().c()) == "-x,-y,-z"
  assert str(space_group_type("I41").change_of_hand_op().c()) == "-x+1/2,-y,-z"
  g = space_group(sgtbx.space_group_symbols("C c c a :1"))
  t = g.change_basis(g.z2p_op()).type()
  assert t.hall_symbol() == "-C 2a 2ac (x-y-1/4,x+y-3/4,z+1/4)"
  assert t.hall_symbol(1) == "-C 2a 2ac (x-y-1/4,x+y-3/4,z+1/4)"
  assert t.hall_symbol(0) == "-C 2a 2ac (x-y-1/4,x+y+1/4,z+1/4)"
  assert t.lookup_symbol() == "Hall: -C 2a 2ac (x-y-1/4,x+y-3/4,z+1/4)"
  assert g.type().lookup_symbol() == "C c c a :1"
  p = pickle.dumps(t)
  l = pickle.loads(p)
  assert t.group() == l.group()

def exercise_phase_info():
  phase_info = sgtbx.phase_info
  g = sgtbx.space_group("P 61 (1 2 0)")
  p = phase_info(g, (1,2,3))
  assert p.sys_abs_was_tested()
  p = phase_info(g, (1,2,3), 0)
  assert p.sys_abs_was_tested()
  p = phase_info(g, (1,2,3), 1)
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
  assert p.ht_angle(0) == p.ht_angle()
  assert approx_equal(p.ht_angle(1), p.ht_angle()*180/math.pi)
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
    assert p.is_valid_phase(phi*180/math.pi, 1)
    assert p.is_valid_phase(phi, 0, 1.e-6)
    phi = p.ht_angle() + i * math.pi + math.pi / 180
    assert not p.is_valid_phase(phi)
    assert not p.is_valid_phase(phi*180/math.pi, 1)
    assert not p.is_valid_phase(phi, 0, 1.e-6)
    assert p.is_valid_phase(phi, 0, 1.e6)
    assert p.is_valid_phase(p.nearest_valid_phase(phi))
    assert p.is_valid_phase(p.nearest_valid_phase(phi, 0), 0)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+1.e-15, 0), 0)
    assert p.is_valid_phase(p.nearest_valid_phase(phi-1.e-15, 0), 0)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+math.pi+1.e-15, 0), 0)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+math.pi-1.e-15, 0), 0)
    assert p.is_valid_phase(p.nearest_valid_phase(phi, 1), 1)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+1.e-15, 1), 1)
    assert p.is_valid_phase(p.nearest_valid_phase(phi-1.e-15, 1), 1)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+180+1.e-15, 1), 1)
    assert p.is_valid_phase(p.nearest_valid_phase(phi+180-1.e-15, 1), 1)
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
  s = site_symmetry(u, g, (0.05,0,0), 0.5, 1)
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
  u = uctbx.unit_cell((3,3,5,90,90,120))
  g = sgtbx.space_group("-P 3 2")
  s = site_symmetry(u, g, (0,0,0))
  a = (5,2,3,-1,2,-2)
  assert not s.is_compatible_u_star(a)
  assert s.is_compatible_u_star(a, 1.e6)
  a = s.average_u_star(a)
  assert s.is_compatible_u_star(a)

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
  assert str(m.sym_op()) in ("x+3,y-1,z+5", "x+3,-y+1,z+5")
  assert approx_equal(m.representative_site(), (0.5,0,0))
  assert approx_equal(m.exact_site(), (-2.5, 1.0, -5.0))
  assert approx_equal(m.distance_moved(), u.length((-0.02,0.123,-0.03)))
  assert str(m.special_op()) == "-5/2,1,-5"
  m = w.mapping(u, x, 0.2)
  assert m.position().letter() == "o"
  assert str(m.sym_op()) == "x+3,y-2,z+5"
  assert approx_equal(m.representative_site(), (0.5,1.123-2,0))
  assert approx_equal(m.exact_site(), (-2.5,1.123,-5))
  assert approx_equal(m.distance_moved(), u.length((-0.02,0,-0.03)))
  assert str(m.special_op()) == "-5/2,y,-5"

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
      for i in e.coordinates().indices():
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
    for i in xrange(3):
      if (i == 0): e = sym_equiv_sites(u, g, x)
      if (i == 1): e = sym_equiv_sites(u, g, x, 0.1)
      if (i == 2): e = sym_equiv_sites(u, g, x, 0.1, 0.001)
      assert e.unit_cell().is_similar_to(u)
      assert e.space_group() == g
      assert e.original_site() == x
      assert not e.special_op().is_valid()
      assert e.max_accepted_tolerance() >= 0
      assert e.coordinates().size() == mult
    e = sym_equiv_sites(u, g, x, 100, 50)
    assert e.coordinates().size() == 1
    for i in xrange(3):
      if (i == 0): e = sym_equiv_sites(ss)
      if (i == 1): e = sym_equiv_sites(wm)
      if (i == 2): e = sym_equiv_sites(u, g, x)
      d = sgtbx.min_sym_equiv_distance_info(e, x)
      assert str(d.sym_op()) == "x,y,z"
      assert d.continuous_shifts() == (0,0,0)
      assert approx_equal(d.diff(), (0,0,0))
      assert approx_equal(d.dist(), 0)
      assert approx_equal(d.sym_op() * x, x)
      a = d.apply(e.coordinates())
      assert a.size() == e.coordinates().size()
      for i in e.coordinates().indices():
        assert approx_equal(a[i], e.coordinates()[i])
      for i,y in e.coordinates().items():
        d = sgtbx.min_sym_equiv_distance_info(e, y)
        assert approx_equal(d.dist(), 0)
        assert approx_equal(d.sym_op() * y, x)
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
    for i in c.indices():
      assert approx_equal(g(i) * x, c[i])
  g = sgtbx.space_group()
  e = sym_equiv_sites(g, (0.016, 0.895, 0.111), uctbx.unit_cell(()))
  d = sgtbx.min_sym_equiv_distance_info(e, (0.939, 0.128, 0.178))
  assert approx_equal(d.dist()**2, 0.064707)
  for hall_symbol,shift_flags in (("P 1", (1,1,1)),
                                  ("P 2", (0,0,1)),
                                  ("P -2", (1,1,0)),
                                  ("P 3x", (1,0,0)),
                                  ("P 3 2", (0,0,0))):
    for x in ((1.732, -1.414, 2.236), (0.939, 0.128, 0.178)):
      g = sgtbx.space_group(hall_symbol)
      e = sym_equiv_sites(g, x, uctbx.unit_cell(()))
      d = sgtbx.min_sym_equiv_distance_info(e, x, shift_flags)
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
          assert not approx_equal(d.sym_op() * z, x)
        else:
          assert approx_equal(d.sym_op() * z, x)
        fz = flex.vec3_double(1, z)
        assert approx_equal(d.apply(fz)[0], x)
  g = sgtbx.space_group()
  e = sym_equiv_sites(g, (0.1,0.2,0.3), u)
  d = sgtbx.min_sym_equiv_distance_info(e, (0.2,0.4,0.1))
  assert approx_equal(d.diff(), (-0.1,-0.2,0.2))
  assert approx_equal(d.dist(), u.length((-0.1,-0.2,0.2)))

def exercise_seminvariant():
  space_group = sgtbx.space_group
  structure_seminvariant = sgtbx.structure_seminvariant
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
    ss = structure_seminvariant(g)
    assert ss.size() == len(ss.vectors_and_moduli())
    assert [(vm.v, vm.m) for vm in ss.vectors_and_moduli()] == vfy
  ss = structure_seminvariant(space_group("P 2 -2"))
  assert ss.is_ss((2,-2,0))
  assert ss.apply_mod((2,-2,0)) == (0,0,0)
  assert not ss.is_ss((2,1,0))
  assert ss.apply_mod((2,1,0)) == (0,1,0)
  ss = structure_seminvariant(space_group("C 2y"))
  assert ss.is_ss((2,0,4))
  assert ss.apply_mod((2,0,4)) == (0,0)
  assert not ss.is_ss((2,1,4))
  assert ss.apply_mod((2,1,4)) == (1,0)
  ss = structure_seminvariant(space_group("C 2 2"))
  assert ss.is_ss((4,2,6))
  assert ss.apply_mod((4,2,6)) == (0,0)
  assert not ss.is_ss((3,2,4))
  assert ss.apply_mod((3,2,4)) == (1,0)
  assert ss.apply_mod((2,3,5)) == (0,1)
  assert ss.apply_mod((3,3,5)) == (1,1)
  ss = structure_seminvariant(space_group("I 4"))
  assert ss.is_ss((1,3,0))
  assert ss.apply_mod((1,3,0)) == (0,)
  assert not ss.is_ss((1,3,4))
  assert ss.apply_mod((1,3,4)) == (4,)
  ss = structure_seminvariant(space_group("I 2 3"))
  assert ss.is_ss((1,3,4))
  assert ss.apply_mod((1,3,4)) == ()
  ss = structure_seminvariant(space_group("I -4"))
  assert ss.gridding() == (2,1,4)
  assert ss.refine_gridding((3,2,10)) == (6,2,20)
  ss = structure_seminvariant(space_group("P -2x"))
  a = ss.grid_adapted_moduli((3,5,7))
  assert [vm.m for vm in a] == [2,5,7]
  ss = structure_seminvariant(space_group("P 3*"))
  a = ss.grid_adapted_moduli((3,5,7))
  assert [vm.m for vm in a] == [3*5*7]
  a = ss.grid_adapted_moduli((3,5,12))
  assert [vm.m for vm in a] == [60]

def exercise_row_echelon():
  m = flex.int((1,1,1,1))
  m.resize(flex.grid(2,2))
  t = flex.int((2,3))
  t.resize(flex.grid(2,1))
  assert sgtbx.row_echelon_form_t(m, t) == 1
  assert m.focus() == (1,2)
  assert tuple(m) == (1,1)
  assert tuple(t) == (2,1)
  assert sgtbx.row_echelon_form(m) == 1
  assert m.focus() == (1,2)
  assert tuple(m) == (1,1)
  m = flex.int((0,-24,0,0,0,-24,24,0,24))
  m.resize(flex.grid(3,3))
  t = flex.int((-3, -6, 0))
  t.resize(flex.grid(3,1))
  assert sgtbx.row_echelon_form_t(m, t) == 3
  assert tuple(m) == (24,0,24,0,24,0,0,0,24)
  assert tuple(t) == (0,3,6)
  t.resize(flex.grid(3))
  sol = flex.int(3)
  assert sgtbx.row_echelon_back_substitution(m, t, sol) == 8
  assert tuple(sol) == (-2,1,2)
  indep = flex.bool((0001,0001,0001))
  assert sgtbx.row_echelon_back_substitution(m, indep=indep) == 1
  assert tuple(indep) == (00000,00000,00000)

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
  exercise_row_echelon()
  print "OK"

if (__name__ == "__main__"):
  run()
