import iotbx.pdb
from libtbx.test_utils import approx_equal
import math
from cctbx import maptbx
from scitbx.array_family import flex
import inspect
from cctbx import crystal

def run(lines, R=5, is_periodic=True, target_cart=[0,0,0]):
  pdb_inp = iotbx.pdb.input(lines = lines, source_info = None)
  xrs = pdb_inp.xray_structure_simple()
  fc = xrs.structure_factors(d_min = 1).f_calc()
  fft_map = fc.fft_map(resolution_factor=0.2)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  for sf in xrs.sites_frac():
    mv = map_data.eight_point_interpolation(sf)
  O = maptbx.MapPeakLocator(
    map_data=map_data,
    unit_cell=fc.unit_cell(),
    is_periodic=is_periodic,
    threshold=1)
  nearby_peaks_cart, nearby_peaks_heights = O.get_peaks_within_radius(
    target_cart=target_cart, R=R)
  return nearby_peaks_cart

def call_test_01():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  5.00           O
"""
  for is_periodic in [True, False]:
    sites = run(lines = pdb_str, is_periodic=is_periodic)
    cntr = 0
    for i, site in enumerate(sites):
      site = [round(_,3) for _ in site]
      assert approx_equal(site, [0.0, 0.0, -0.0])
      cntr+=1
      print(i, site)
    assert cntr == 1, cntr

def call_test_02():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  5.00           O
HETATM    2  O   HOH A   2       3.000   3.000   3.000  1.00  5.00           O
"""
  for is_periodic in [True, False]:
    sites = run(lines = pdb_str, R=6, is_periodic=is_periodic)
    cntr = 0
    for i, site in enumerate(sites):
      site = [round(_,3) for _ in site]
      if i==0:
        assert approx_equal(site, [3.001, 3.001, 3.001])
        cntr+=1
      if i==1:
        assert approx_equal(site, [-0.001, -0.001, -0.001])
        cntr+=1
      print(i, site)
    assert cntr == 2, cntr

def call_test_03():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    2  O   HOH A   2      -3.000   0.000   0.000  1.00  5.00           O
"""
  sites = run(lines = pdb_str, is_periodic=True)
  cntr = 0
  for i, site in enumerate(sites):
    site = [round(_,3) for _ in site]
    assert approx_equal(site, [-3.0, 0.0, 0.0])
    cntr+=1
    print(i, site)
  assert cntr == 1, cntr

def call_test_04():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    2  O   HOH A   2      -3.000   0.000   0.000  1.00  5.00           O
"""
  sites = run(lines = pdb_str, is_periodic=False)
  assert len(sites) == 0

def call_test_05():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00  5.00           O
"""
  for is_periodic in [True, False]:
    sites = run(lines = pdb_str, is_periodic=is_periodic)
    cntr = 0
    for i, site in enumerate(sites):
      site = [round(_,3) for _ in site]
      assert approx_equal(site, [0.0, 0.0, -0.0])
      cntr+=1
      print(i, site)
    assert cntr == 1, cntr

def call_test_06():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1      10.100  10.200  10.300  1.00  5.00           O
"""
  for is_periodic in [True, False]:
    sites = run(lines = pdb_str, is_periodic=is_periodic)
    cntr = 0
    for i, site in enumerate(sites):
      site = [round(_,3) for _ in site]
      assert approx_equal(site, [0.105, 0.2, 0.305]) # peak is at [0.1, 0.2, 0.3]
      cntr+=1
      print(i, site)
    assert cntr == 1, cntr

def call_test_07():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1      10.100  10.200  10.300  1.00  5.00           O
"""
  for is_periodic in [True, False]:
    sites = run(lines = pdb_str, is_periodic=is_periodic, target_cart=[10,10,10])
    cntr = 0
    for i, site in enumerate(sites):
      site = [round(_,3) for _ in site]
      if is_periodic:
        assert approx_equal(site, [10.105, 10.2, 10.305])
      cntr+=1
      print(i, site)
    if is_periodic: assert cntr == 1, cntr
    else:           assert cntr == 0, cntr

def call_test_08():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       4.000   4.000   4.000  1.00  5.00           O
HETATM    1  O   HOH A   1       6.000   6.000   6.000  1.00  5.00           O
"""
  for is_periodic in [True, False]:
    sites = run(lines = pdb_str, is_periodic=is_periodic, target_cart=[5,5,5], R=5)
    cntr = 0
    for i, site in enumerate(sites):
      site = [round(_,3) for _ in site]
      if i==0:
        assert approx_equal(site, [6.0, 6.0, 6.0])
        cntr+=1
      if i==1:
        assert approx_equal(site, [4.0, 4.0, 4.0])
        cntr+=1
      print(i, site)
    assert cntr == 2, cntr

def call_test_09():
  print(inspect.currentframe().f_code.co_name)
  pdb_str="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       4.000   4.000   4.000  1.00  5.00           O
HETATM    1  O   HOH A   1       6.000   6.000   6.000  1.00  5.00           O
"""
  for is_periodic in [True, False]:
    sites = run(lines = pdb_str, is_periodic=is_periodic, target_cart=[5,5,5], R=1)
    assert len(sites) == 0

def call_test_10():
  print(inspect.currentframe().f_code.co_name)
  uc = sym = crystal.symmetry(
    unit_cell=(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
    space_group_symbol="P1").unit_cell()
  md = flex.double(flex.grid(10,10,10), 0)
  md[3,4,5] = 17
  O = maptbx.MapPeakLocator(
    map_data    = md,
    unit_cell   = uc,
    is_periodic = False,
    threshold   = 10)
  for R in [1, 5]:
    sites, h = O.get_peaks_within_radius(target_cart=[3,4,5], R=1)
    assert sites.size()==1
    assert approx_equal(sites[0], [3,4,5])
    assert approx_equal(h[0], 17)

def call_test_11():
  print(inspect.currentframe().f_code.co_name)
  uc = sym = crystal.symmetry(
    unit_cell=(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
    space_group_symbol="P1").unit_cell()
  md = flex.double(flex.grid(10,10,10), 0)
  md[0,0,0] = 17
  O = maptbx.MapPeakLocator(
    map_data    = md,
    unit_cell   = uc,
    is_periodic = False,
    threshold   = 10)
  for R in [1, 5]:
    sites, h = O.get_peaks_within_radius(target_cart=[0,0,0], R=1)
    assert sites.size()==1
    assert approx_equal(sites[0], [0,0,0])
    assert approx_equal(h[0], 17)

def call_test_12():
  print(inspect.currentframe().f_code.co_name)
  uc = sym = crystal.symmetry(
    unit_cell=(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
    space_group_symbol="P1").unit_cell()
  for point in [[0,0,0], [1,2,3]]:
    print("poit:", point)
    md = flex.double(flex.grid(10,10,10), 0)
    md[point] = 17
    O = maptbx.MapPeakLocator(
      map_data    = md,
      unit_cell   = uc,
      is_periodic = False,
      threshold   = 10)
    sites, h = O.get_peaks_within_radius(target_cart=[-0.1,-0.2,-0.3], R=7)
    print(list(sites))
    assert sites.size()==1
    assert approx_equal(sites[0], point)
    assert approx_equal(h[0], 17)

def call_test_13():
  """
  AI suggested to test its own code with the trieckiest test!
  """
  print(inspect.currentframe().f_code.co_name)
  # 1. Setup a standard 10x10x10 unit cell
  uc = crystal.symmetry(
      unit_cell=(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
      space_group_symbol="P1").unit_cell()
  # 2. Initialize map and place a massive peak exactly at the origin
  md = flex.double(flex.grid(10, 10, 10), 0)
  md[(0, 0, 0)] = 100.0
  # 3. Instantiate BOTH versions of the locator on the exact same map
  locator_cryo = maptbx.MapPeakLocator(md, uc, is_periodic=False, threshold=10)
  locator_xtal = maptbx.MapPeakLocator(md, uc, is_periodic=True, threshold=10)
  # 4. THE TRAP: Target a Cartesian point on the FAR opposite corner of the map.
  target_cart = [9.9, 9.9, 9.9]
  # Mathematical reality:
  # - Strict Cartesian distance from (9.9, 9.9, 9.9) to (0, 0, 0) is ~17.15 Å
  # - Periodic shortest distance (wrapped) is distance to (10, 10, 10), which
  # is ~0.17 Å
  dist_cartesian = math.sqrt(9.9**2 + 9.9**2 + 9.9**2)
  dist_periodic  = math.sqrt(0.1**2 + 0.1**2 + 0.1**2)
  # 5. Query both locators with a small radius (R=2.0)
  R = 2.0
  sites_cryo, h_cryo = locator_cryo.get_peaks_within_radius(target_cart, R)
  sites_xtal, h_xtal = locator_xtal.get_peaks_within_radius(target_cart, R)
  # 6. THE TRICKY ASSERTIONS
  # A) For Crystallography (Periodic), the peak SHOULD be found!
  # Because 0.17 Å < 2.0 Å, the target dynamically wraps around the unit cell.
  assert sites_xtal.size() == 1, \
      "Xtal Test Failed: Periodic locator missed the wrapped peak!"
  # B) For CryoEM (Non-periodic), the peak MUST NOT be found!
  # Because 17.15 Å > 2.0 Å, strict Euclidean distance must drop it.
  # If sites_cryo.size() == 1, it means cctbx C++ wrapping bled into the final
  # results!
  assert sites_cryo.size() == 0, \
      "CryoEM Test Failed: Locator found an aliased ghost peak across the box!"
  # 7. Final verification that the periodic peak found was translated correctly
  # The periodic engine should translate the (0,0,0) peak to sit near the target
  assert approx_equal(sites_xtal[0], (10.0, 10.0, 10.0), eps=1e-2)

if(__name__ == "__main__"):
  call_test_01()
  call_test_02()
  call_test_03()
  call_test_04()
  call_test_05()
  call_test_07()
  call_test_08()
  call_test_09()
  call_test_10()
  call_test_11()
  call_test_12()
  call_test_13()
