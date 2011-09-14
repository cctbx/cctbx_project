from __future__ import division

test_cases = """\
P 1
-0.13463, 0.5133, 0.41428
x,y,z

P -1
-1, -0.5, 0.5
-1,-1/2,1/2

F d d 2 (-x-y+z,2*x,2*y)
0.51827, 0, 0
x,0,0

C 1 2/c 1 (-x+y,z,2*x-z)
-0.36537, 0.75, 0.25
x,3/4,1/4

C m c a (z,x+y,2*y)
0, -0.20681, 0
0,y,0

C 1 2/c 1
0.5, -0.13463, 0.25
1/2,y,1/4

F d d 2
0, 0, 0.51827
0,0,z

F d d 2
-0.25, 0.25, 0.5133
-1/4,1/4,z

C m c 21
0, 0.41428, 0.8342
0,y,z

P m c n
0.25, 0.209, 0.2009
1/4,y,z

C 1 m 1
0.2479, 0, 0.8201
x,0,z

P 1 21/m 1
0.65486, 0.75, 0.46088
x,3/4,z

C 1 m 1 (z,x-y,2*y)
0.8201, -0.2479, 0
x,y,0

P 1 1 21/m
-0.65486, -0.46088, -0.75
x,y,-3/4

A b a 2 (y+z,2*z,x)
0.21585, 0.4317, 0
x,2*x,0

R -3 m :H
0.38797, 0.193985, -0.00711
x,1/2*x,z

F d d 2 (x+y+z,x+y-z,2*y)
-0.0133, 1.0133, 0.5
x,-x+1,1/2

I 4/m c m
0.85106, 0.35106, 0
x,x-1/2,0

C 1 2/c 1 (x+y,z,-2*y)
-0.29421, -0.25, 0.58842
x,-1/4,-2*x

C 1 2/m 1 (z,x-y,2*y)
0, -0.3019, -0.3962
0,y,-2*y-1

C 1 2/m 1 (z,x-y,x+y)
-0.5, 1.38934, 0.61066
-1/2,y,-y+2

F -4 3 m
0.10515, 0.05925, 0.10515
x,y,x

C m c 21 (-x+y,z,2*y)
0.148327, 0.5362, -0.703346
x,y,2*x-1

R -3 m :R
0.38086, -0.201095, -0.201095
x,y,y

C 1 2/m 1 (-z,x+y,2*x)
0.74518, 1.4919, 1.9838
x,y,2*y-1

C m c m (x+y,-x+y,z)
0.4689, 0.4689, 0.10868
x,x,z

C m c a (x+y,2*y,z)
-0.347957, 0.304086, -0.1964
x,2*x+1,z

P -4 b 2
0.88761, 0.61239, 0.5
x,-x+3/2,1/2

R -3 m :H (x+z,-y+z,-3*z)
0.3102, -0.1548, -0.0006
x,y,-x-2*y NO_FVAR

C 1 2/m 1 (-x+y,z,2*x-z)
-0.17213, -0.09996, 0.44422
x,y,-2*x-y NO_FVAR

C 1 2/m 1 (x+y,-x+y+z,z)
-0.06876, -0.216, -0.28476
x,y,x+y NO_FVAR

C 1 2/m 1 (-x+y,z,2*x-z)
-1.34469, 0.83463, 0.85475
x,y,-2*x-y-1 NO_FVAR

C 1 2/m 1 (-x+y,z,2*x-z)
-1.70949, 0.67904, 0.73994
x,y,-2*x-y-2 NO_FVAR

F -4 3 m (-x+y+z,x-y+z,x+y-z)
-0.151837, -0.001949, 0.576893
x,y,-1/2*x-1/2*y+1/2 NO_FVAR

C 1 2/m 1 (x+y,-x+y+z,z)
-0.3761, -0.8643, -0.2404
x,y,x+y+1 NO_FVAR

F -4 3 m (-x+y+z,x-y+z,x+y-z)
0.324012, 0.226578, 0.125398
x,y,-2*x-y+1 NO_FVAR

F -4 3 m (-x+y+z,x-y+z,x+y-z)
-0.384945, -0.110849, 0.247897
x,y,-1/2*x-1/2*y NO_FVAR

R -3 c :R
0.0639, 0.0639, 0.0639
x,x,x

F -4 3 m
0.208141, 0.208141, 0.291859
x,x,-x+1/2

F -4 3 m
0.1718, 0.3282, -0.1718
x,-x+1/2,-x

R -3 :H (x+z,-y+z,-3*z)
0.206827, -0.793173, 1.37952
x,x-1,-3*x+2 NO_FVAR

P 21 3
0.25347, 0.24653, 0.75347
x,-x+1/2,x+1/2 NO_FVAR

"""

def test_case_generator():
  from cctbx.sgtbx import space_group_info
  lines = test_cases.splitlines()
  assert len(lines) % 4 == 0
  for i in xrange(0, len(lines), 4):
    assert len(lines[i+3]) == 0
    expr = lines[i+2]
    j = expr.find(" NO_FVAR")
    no_fvar = (j >= 0)
    if (no_fvar):
      expr = expr[:j]
    yield \
      space_group_info(symbol=lines[i]), \
      eval(lines[i+1]), \
      expr, \
      no_fvar

def exercise():
  from scitbx.array_family import flex
  from scitbx import matrix
  from libtbx.test_utils import approx_equal, show_diff
  import libtbx.load_env
  if (libtbx.env.has_module(name="iotbx")):
    from iotbx.shelx.parsers import decode_variables
    import iotbx.shelx.errors
  else:
    decode_variables = None
  mt = flex.mersenne_twister(seed=0)
  for sgi,site,expected_expr,no_fvar in test_case_generator():
    sps = sgi.any_compatible_crystal_symmetry(volume=1000) \
      .special_position_settings()
    ss = sps.site_symmetry(site)
    assert approx_equal(ss.exact_site(), site)
    sos = ss.special_op_simplified()
    expr = str(sos)
    assert not show_diff(expr, expected_expr)
    def check(site):
      ns = dict(zip("xyz", site))
      expr_site = eval(expr, ns, {})
      assert approx_equal(expr_site, site, 1e-4)
      #
      shifted_site = (
        matrix.col(site) + matrix.col(mt.random_double_point_on_sphere()))
      ns = dict(zip("xyz", shifted_site))
      expr_shifted_site = eval(expr, ns, {})
      expr_shifted_site_exact = ss.special_op() * expr_shifted_site
      assert approx_equal(expr_shifted_site, expr_shifted_site_exact)
    check(site)
    for i_trial in xrange(3):
      shifted_site = ss.special_op() * (
        matrix.col(site) + matrix.col(mt.random_double_point_on_sphere()))
      check(shifted_site)
    fvars = [None] # placeholder for scale factor
    if (decode_variables is not None):
      try:
        coded_variables = ss.shelx_fvar_encoding(site=site, fvars=fvars)
      except iotbx.shelx.errors.error:
        if (not no_fvar): raise
      else:
        assert not no_fvar
      if (not no_fvar):
        values, behaviors = decode_variables(
          free_variable=fvars, coded_variables=coded_variables)
        mismatch = sps.unit_cell().mod_short_distance(site, values)
        assert mismatch < 1e-10

def run(args):
  assert len(args) == 0
  exercise()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
