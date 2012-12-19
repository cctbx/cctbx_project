from __future__ import division
from scitbx.array_family import flex
import iotbx.pdb
from cctbx import maptbx
from cctbx import adptbx, sgtbx
import scitbx
import time
from libtbx.test_utils import approx_equal

def ru(crystal_symmetry, u_scale=1, u_min=0.1):
  from cctbx import sgtbx
  symbol = crystal_symmetry.space_group().type().lookup_symbol()
  point_group = sgtbx.space_group_info(
    symbol=symbol).group().build_derived_point_group()
  adp_constraints = sgtbx.tensor_rank_2_constraints(
    space_group=point_group,
    reciprocal_space=True)
  u_star = adptbx.u_cart_as_u_star(crystal_symmetry.unit_cell(),
    adptbx.random_u_cart(u_scale=u_scale,u_min=u_min))
  u_indep = adp_constraints.independent_params(all_params=u_star)
  u_star = adp_constraints.all_params(independent_params=u_indep)
  r = flex.sym_mat3_double()
  r.append(adptbx.u_star_as_u_cart(crystal_symmetry.unit_cell(), u_star))
  return r

def cartesian_space_version(site_cart, map_data, unit_cell, radius=2):
  sites = flex.vec3_double()
  values = flex.double()
  center_frac = unit_cell.fractionalize(site_cart)
  cut = map_data.eight_point_interpolation(center_frac)/3.
  ix=-radius
  step=0.5
  while ix<=radius:
    iy=-radius
    while iy<=radius:
      iz=-radius
      while iz<=radius:
        site_cart_ = [ix+site_cart[0],iy+site_cart[1],iz+site_cart[2]]
        site_frac = unit_cell.fractionalize(site_cart_)
        mv = map_data.eight_point_interpolation(site_frac)
        if(mv>cut):
          sites.append([ix,iy,iz])
          values.append(mv)
        iz+=step
      iy+=step
    ix+=step
  #
  v = scitbx.math.principal_axes_of_inertia(
    points=sites,
    weights=values).eigensystem().values()
  #
  if 0: # for debugging
    of = open("tmp.pdb", "w")
    print >> of, "CRYST1   10.000   10.000   10.000  70.00 100.00 120.00 P 1"
    for i,site in enumerate(sites):
      fmt = "HETATM%5d  O   HOH %5d    %8.3f%8.3f%8.3f  1.00%6.2f           O"
      print >> of, fmt%(i,i,site[0],site[1],site[2], values[i])
    of.close()
  #
  v_max = flex.max(v)
  if(v_max!=0):
    return flex.min(v)/v_max
  else: return 0

def run():
  pdb_str1="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       5.000   5.000   5.000  1.00  0.00           C
ANISOU    1  C    C      1     1000  10000   5000      0      0      0       C
END
"""
  pdb_str2="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       5.000   7.000  -9.000  1.00  0.00           C
ANISOU    1  C    C      1    10000  10000   5000      0      0      0       C
END
"""
  pdb_str3="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       5.000   7.000  -9.000  1.00  0.00           C
ANISOU    1  C    C      1    10000  10000  10000      0      0      0       C
END
"""
  pdb_str4="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       5.000   7.000  -9.000  1.00  0.00           C
ANISOU    1  C    C      1        0      0      0      0      0      0       C
END
"""
  pdb_str5="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       5.000   7.000  -9.000  1.00  0.00           C
ANISOU    1  C    C      1        1      1      0      0      0      0       C
END
"""
  pdb_str6="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       5.000   7.000  -9.000  1.00 20.00           C
END
"""
  expected_a = [0.100, 0.500, 1.000, 1.000, 0.000, 1.000]
  expected_s = [0.161, 0.557, 0.999, 0.904, 0.907, 0.997]
  result_a = []
  result_s = []
  for i, pdb_str in enumerate([pdb_str1, pdb_str2, pdb_str3, pdb_str4, pdb_str5,
                               pdb_str6]):
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    ph.write_pdb_file(file_name = "m.pdb")
    # To try more random adp values:
    #u_cart = ru(crystal_symmetry=xrs.crystal_symmetry(), u_scale=1, u_min=0.01)
    #xrs.set_u_cart(u_cart = u_cart)
    #
    fc = xrs.structure_factors(d_min = 0.26, algorithm = "direct").f_calc()
    fft_map = fc.fft_map(resolution_factor = 0.25)
    fft_map.apply_volume_scaling()
    map_data = fft_map.real_map_unpadded()
    # Implementation that intuitively should work, but it doesn't for unclear to
    # me reasons
    alt = cartesian_space_version(site_cart=xrs.sites_cart()[0],
      map_data=map_data, unit_cell=xrs.unit_cell(), radius=2)
    # Calculated anisotropy from ADP: this is exact answer
    a = list(xrs.scatterers().anisotropy(xrs.unit_cell()))[0]
    result_a.append(a)
    # Sphericity: main working implementation
    t0=time.time()
    s = maptbx.sphericity(map_data = map_data, unit_cell = xrs.unit_cell(),
      radius = 2., sites_frac = xrs.sites_frac())[0]
    result_s.append(s)
    t1=time.time()
    # Alternative implementation that doesn't seem to work correctly
    pa = maptbx.principal_axes_of_inertia(real_map=map_data,
      site_cart=xrs.sites_cart()[0], unit_cell=xrs.unit_cell(), radius=2.0)
    v = pa.eigensystem().values()
    n = flex.min(v)/flex.max(v)
    t2=time.time()
    #
    time_gain = (t2-t1)/(t1-t0)
    print "#%d: answer: %5.3f sphericity: %5.3f other: %5.3f time_gain: %6.2f"%(
      i+1, a, s, n, time_gain), alt
  #
  assert approx_equal(expected_a, result_a, 1.e-3)
  assert approx_equal(expected_s, result_s, 1.e-3)

if __name__ == "__main__" :
  run()
