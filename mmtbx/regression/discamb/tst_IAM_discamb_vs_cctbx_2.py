from __future__ import absolute_import, division, print_function
import time, random
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
from cctbx.array_family import flex
#from cctbx.eltbx.discamb import tst_IAM_discamb_vs_cctbx as tst_disc
import pydiscamb

'''
Notes
-----
This test can only be executed if pyDiSCaMB is available
Right now, it needs to be installed manually, but it will be
included in bootstrap in the near future.
'''

def compare_structure_factors(x,y):
  x = flex.abs(x)
  y = flex.abs(y)
  scale = flex.sum(x*y)/flex.sum(y*y)
  num = flex.sum(flex.abs(x-scale*y))
  den = flex.sum(flex.abs(x+scale*y))
  diff = flex.abs(x-y)
  return num/den*2*100., flex.mean(diff), flex.max(diff)

def test_all_spacegroups():
  """
  Test structure factor calculations for all space groups.

  This function iterates over all 230 space groups, generates random crystal
  structures with a set of common elements, isotropic atomic displacement
  parameters, and occupancies, and compares structure factors calculated with
  cctbx and DiSCaMB.

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If the relative difference score or absolute differences exceed predefined
      thresholds.

  Notes
  -----
  For each space group, this function:
  - Creates a random crystal structure with certain elements.
  - Optionally applies random perturbations to fp or fdp.
  - Computes structure factors with cctbx and DiSCaMB
  - resolution limit is random between 2 and 5 A
  - Asserts that the difference between cctbx and DiSCaMB structure factors
    remains within acceptable limits.
  """
  for sg_number in range(1,231):
    sgi = space_group_info(sg_number)
    elements = ["C", "O", "N", "H", "S", "Na", "Fe", "Ca"] * 10
    random_occ = random.choice([True, False])
    #
    xrs = random_structure.xray_structure(
        space_group_info       = sgi,
        elements               = elements,
        general_positions_only = False,
        use_u_iso              = True,
        random_u_iso = True,
        random_occupancy=random_occ,
    )
    random_fp = random.choice([True, False])
    if random_fp:
      xrs.shake_fps()
    random_fdp = random.choice([True, False])
    if random_fdp:
      xrs.shake_fdps()
    xrs.scattering_type_registry(table = "electron")
    d_min = random.uniform(2, 5)
    # Calculate structure factors with cctbx
    fcalc_cctbx = xrs.structure_factors(
      d_min=d_min, algorithm="direct").f_calc().data()
    # Calculate structure factors with DiSCaMB
    fcalc_discamb = flex.complex_double(
      pydiscamb.calculate_structure_factors_IAM(xrs, d_min))
    assert fcalc_discamb.size() == fcalc_cctbx.size()
    # compute the difference between the structure factors
    score, mean_diff, max_diff = compare_structure_factors(
      x=fcalc_cctbx, y=fcalc_discamb)
    #print (sgi, score, mean_diff, max_diff)
    # Check if the differences are within the specified tolerance
    assert(score < 0.0003)
    assert(mean_diff < 0.0015)
    assert(max_diff < 0.011)

if (__name__ == "__main__"):
  """
  Run structure factor tests for all space groups.
  Takes 40-60 seconds.
  """
  t0 = time.time()
  test_all_spacegroups()
  print("OK. Time: %8.3f"%(time.time()-t0))
