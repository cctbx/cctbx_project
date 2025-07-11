from __future__ import absolute_import, division, print_function
import time
from libtbx.test_utils import approx_equal
from cctbx.eltbx import e_scattering
from cctbx.eltbx import xray_scattering
import pydiscamb

def test_all_coefficients_electron():
  """
  Test coefficient values for electron scattering factors across all elements.

  This function retrieves scattering factor coefficients for all elements from
  the cctbx electron scattering table and compares them with values from the
  DiSCaMB electron scattering table.

  Returns
  -------
  None

  Raises
  ------
  AssertionError
    If any coefficient (a, b, or c) differs between cctbx and DiSCaMB.
  """
  elements = e_scattering.ito_vol_c_2011_table_4_3_2_2_elements()
  table_discamb = pydiscamb.get_table('electron')
  for el in elements:
    entry = e_scattering.ito_vol_c_2011_table_4_3_2_2_entry_as_gaussian(label=el)
    entry_disc = table_discamb[el]
    #print(entry.array_of_a(), entry.array_of_b(), entry.c())
    #print(entry_disc.a, entry_disc.b, entry_disc.c)
    assert set(entry_disc.a) == set(entry.array_of_a())
    assert set(entry_disc.b) == set(entry.array_of_b())
    assert (entry.c() == entry_disc.c)

# ==============================================================================

def test_all_coefficients_it1992():
  """
  Test coefficient values for X-ray scattering factors based on the it1992 table.

  This function retrieves scattering factor coefficients for all elements from
  the cctbx it1992 scattering table and compares them with values from the
  DiSCaMB it1992 table.

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If any coefficient (a, b, or c) differs between cctbx and DiSCaMB.
  """
  table_discamb = pydiscamb.get_table('it1992')
  for it in xray_scattering.it1992_iterator():
    el = it.label()
    entry = it.fetch()
    entry_disc = table_discamb[el]
    #print(entry.array_of_a(), entry.array_of_b(), entry.c())
    #print(entry_disc.a, entry_disc.b, entry_disc.c)
    assert approx_equal(entry.array_of_a(), tuple(entry_disc.a))
    assert approx_equal(entry.array_of_b(), tuple(entry_disc.b))
    assert approx_equal(entry.c(), entry_disc.c)

# ==============================================================================

def test_all_coefficients_wk1995():
  """
  Test coefficient values for X-ray scattering factors based on the wk1995 table.

  This function retrieves scattering factor coefficients for all elements from
  the cctbx WK1995 scattering table and compares them with values from the
  DiSCaMB WK1995 table.

  Returns
  -------
  None

  Raises
  ------
  AssertionError
      If any coefficient (a, b, or c) differs between cctbx and DiSCaMB.
  """
  table_discamb = pydiscamb.get_table('wk')
  for it in xray_scattering.wk1995_iterator():
    el = it.label()
    if el in ['Hiso', 'O2-']: continue
    entry = it.fetch()
    entry_disc = table_discamb[el]
    #print(entry.array_of_a(), entry.array_of_b(), entry.c())
    #print(entry_disc.a, entry_disc.b, entry_disc.c)
    assert approx_equal(entry.array_of_a(), tuple(entry_disc.a))
    assert approx_equal(entry.array_of_b(), tuple(entry_disc.b))
    assert approx_equal(entry.c(), entry_disc.c)

# ==============================================================================

if (__name__ == "__main__"):
  """
  Run coefficient tests for electron, it1992, and wk1995 scattering tables.

  This main function executes tests for all supported scattering tables and.
  """
  t0 = time.time()
  test_all_coefficients_electron()
  test_all_coefficients_it1992()
  test_all_coefficients_wk1995()
  print("OK. Time: %8.3f"%(time.time()-t0))
