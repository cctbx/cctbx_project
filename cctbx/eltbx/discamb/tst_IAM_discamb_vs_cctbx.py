from __future__ import absolute_import, division, print_function
import time
from cctbx.development import random_structure as cctbx_random_structure
from cctbx.sgtbx import space_group_info
from cctbx.array_family import flex
from cctbx.eltbx import e_scattering
from cctbx.eltbx import xray_scattering
import pydiscamb

'''
Notes
-----
This test can only be executed if pyDiSCaMB is available
Right now, it needs to be installed manually, but it will be
included in bootstrap in the near future.

The test currently fails due to missing values in pyDiSCaMB.
Got in touch with pyDiSCaMB developers.
'''

def compare_structure_factors(x,y):
  """
  Compute the relative mean and maximum differences between two sets of structure factors.

  Parameters
  ----------
  x : flex array
      The first set of structure factor data.
  y : flex array
      The second set of structure factor data.

  Returns
  -------
  score : float
      The relative difference score calculated as (num/den * 2 * 100).
  mean_diff : float
      The mean of the absolute differences between `x` and `y`.
  max_diff : float
      The maximum absolute difference between `x` and `y`.

  Notes
  -----
  The function computes a scaling factor between `x` and `y`, then calculates
  the relative difference using both absolute sums and element-wise absolute differences.
  """
  x = flex.abs(x)
  y = flex.abs(y)
  scale = flex.sum(x*y)/flex.sum(y*y)
  num = flex.sum(flex.abs(x-scale*y))
  den = flex.sum(flex.abs(x+scale*y))
  diff = flex.abs(x-y)
  return num/den*2*100., flex.mean(diff), flex.max(diff)

def test_all_elements_wk1995():
  """
  Test structure factor calculations for all elements using the WK1995 scattering table.

  This function iterates over all elements in the cctbx WK1995 scattering table,
  computes structure factors with cctbx and DiSCaMB, and checks the differences using
  the `compare_structure_factors` function.

  Returns
  -------
  None

  Notes
  -----
  Calls `check_element_in_cctbx_and_dyscamb` for each element, using 'wk1995' as the
  scattering table parameter.
  """
  for e in xray_scattering.wk1995_iterator():
    el = e.label()
    #if el in ['Hiso']: continue
    check_element_in_cctbx_and_dyscamb(
      element = el,
      table   = 'wk1995')


def test_all_elements_it1992():
  """
  Test structure factor calculations for all elements using the IT1992 scattering table.

  This function iterates over all elements in the cctbx IT1992 scattering table,
  computes structure factors with cctbx and DiSCaMB, and checks the differences using
  the `compare_structure_factors` function.

  Returns
  -------
  None

  Notes
  -----
  Calls `check_element_in_cctbx_and_dyscamb` for each element, using 'it1992' as the
  scattering table parameter.
  """
  for e in xray_scattering.it1992_iterator():
    el = e.label()
    #if el in ['Hiso']: continue
    check_element_in_cctbx_and_dyscamb(
      element = el,
      table   = 'it1992')

def test_all_elements_electron():
  """
  Test structure factor calculations for all elements using the electron scattering table.

  This function iterates over all elements in the cctbx electron scattering table,
  generates random crystal structures, calculates structure factors using both cctbx
  and DiSCaMB, and compares the results using the `compare_structure_factors` function.

  Returns
  -------
  None

  Notes
  -----
  Calls `check_element_in_cctbx_and_dyscamb` for each element, using 'electron' as the
  scattering table parameter.
  """
  for el in e_scattering.ito_vol_c_2011_table_4_3_2_2_elements():
    check_element_in_cctbx_and_dyscamb(
      element = el,
      table   = 'electron')


def check_element_in_cctbx_and_dyscamb(element, table=None):
  """
  Compute and compare structure factors for a given element using specified scattering tables.

  Parameters
  ----------
  element : str
      The element symbol for which structure factors are computed.
  table : str, optional
      The scattering table to use for structure factor calculation
      (e.g., 'wk1995', 'it1992', 'electron').

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If the relative difference score or absolute differences exceed predefined thresholds.

  Notes
  -----
  Computes structure factors using both cctbx and DiSCaMB and compares the results
  using `compare_structure_factors`. Assertions ensure differences remain within
  acceptable limits.
  """
  assert table is not None
  group = space_group_info(19)
  xrs = cctbx_random_structure.xray_structure(
      space_group_info=group,
      elements=[element] * 3,
      use_u_aniso = True)
  xrs.scattering_type_registry(table = table)
  d_min = 2
  # Calculate structure factors with cctbx
  fcalc_cctbx = xrs.structure_factors(
    d_min=d_min, algorithm="direct").f_calc().data()
  # Calculate structure factors with DiSCaMB
  fcalc_discamb = flex.complex_double(
    pydiscamb.calculate_structure_factors_IAM(xrs, d_min))

  assert fcalc_discamb.size() == fcalc_cctbx.size()
  # compute the difference of the results
  score, mean_diff, max_diff = compare_structure_factors(
    x=fcalc_cctbx, y=fcalc_discamb)
  print (element, score, mean_diff, max_diff)
  assert(score < 0.0001)
  assert(mean_diff < 0.00001)
  assert(max_diff < 0.0001)

if (__name__ == "__main__"):
  t0 = time.time()
  test_all_elements_electron()
  test_all_elements_it1992()
  test_all_elements_wk1995()
  print("OK. Time: %8.3f"%(time.time()-t0))
