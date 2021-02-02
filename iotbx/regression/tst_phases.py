from __future__ import absolute_import, division, print_function
import iotbx.pdb.xray_structure
from cctbx.development import random_structure
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.test_utils import show_diff
from six.moves import cStringIO as StringIO
import math
import sys

def exercise_basic():
  a = flex.double((10,20))
  p = flex.double((-25,355))
  c = flex.polar(a, p, True)
  f = flex.double((0.3,0.9))
  s = miller.set(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(10,10,10,90,90,90),
      space_group_symbol="P1"),
    indices=flex.miller_index([(1,2,3),(-3,4,-6)]),
    anomalous_flag=False)
  out = StringIO()
  s.array(data=c).as_phases_phs(out=out)
  assert not show_diff(out.getvalue(), """\
   1   2   3 4999.99    1.00  -25.00
  -3   4  -6 9999.99    1.00   -5.00
""")
  out = StringIO()
  s.array(data=c).as_phases_phs(out=out, scale_amplitudes=False)
  assert not show_diff(out.getvalue(), """\
   1   2   3   10.00    1.00  -25.00
  -3   4  -6   20.00    1.00   -5.00
""")
  for phases in [s.array(data=p), p]:
    out = StringIO()
    s.array(data=c).amplitudes().as_phases_phs(
      out=out, phases=phases, phases_deg=True)
    assert not show_diff(out.getvalue(), """\
   1   2   3 4999.99    1.00  -25.00
  -3   4  -6 9999.99    1.00  355.00
""")
  for phases in [s.array(data=p*(math.pi/180)), p*(math.pi/180)]:
    out = StringIO()
    s.array(data=c).amplitudes().as_phases_phs(
      out=out, phases=phases, phases_deg=False)
    assert not show_diff(out.getvalue(), """\
   1   2   3 4999.99    1.00  -25.00
  -3   4  -6 9999.99    1.00  355.00
""")
  for figures_of_merit in [s.array(data=f), f]:
    out = StringIO()
    s.array(data=c).as_phases_phs(out=out, figures_of_merit=figures_of_merit)
    assert not show_diff(out.getvalue(), """\
   1   2   3 4999.99    0.30  -25.00
  -3   4  -6 9999.99    0.90   -5.00
""")

def generate_random_f_calc(
      space_group_info,
      n_elements=10,
      d_min=1.5,
      verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=1000,
    min_distance=3.,
    general_positions_only=False)
  if (0 or verbose):
    structure.show_summary()
    structure.show_scatterers()
    print()
  print("Writing tmp.pdb")
  s = structure.as_pdb_file(
    remark="random structure",
    resname="RND")
  with open("tmp.pdb", "w") as f:
    f.write(s)
  if (0 or verbose):
    print()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  if (0 or verbose):
    f_calc.show_summary()
    print()
  print("Writing: tmp.phs")
  with open("tmp.phs", "w") as f:
    f_calc.as_phases_phs(out=f)
  if (0 or verbose):
    print()

def run():
  exercise_basic()
  if (len(sys.argv) > 1):
    verbose = 0
    args = []
    for arg in sys.argv[1:]:
      if (arg.lower() == "--verbose"):
        verbose = 1
      else:
        args.append(arg)
    for arg in args:
      generate_random_f_calc(sgtbx.space_group_info(arg), verbose=verbose)
  print("OK")

if (__name__ == "__main__"):
  run()
