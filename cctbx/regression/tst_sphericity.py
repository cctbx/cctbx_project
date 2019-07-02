from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
from cctbx import maptbx
from cctbx import adptbx, sgtbx
from libtbx.test_utils import approx_equal
import time

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
HETATM    1  C    C      1       5.000   7.000  -9.000  1.00 20.00           C
END
"""

  pdb_str6="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       5.000   7.000   8.000  1.00  0.00           C
ANISOU    1  C    C      1    10000  10000   5000      0      0      0       C
END
"""
  pdb_str7="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1       0.000   0.000   0.000  1.00  0.00           C
ANISOU    1  C    C      1    10000  10000   5000      0      0      0       C
END
"""
  pdb_str8="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       0.000   0.000   0.000  1.00  0.00           C
ANISOU    1  C    C      1     1000  10000   5000      0      0      0       C
END
"""
  pdb_str9="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1      -9.000   0.000  50.000  1.00  0.00           C
ANISOU    1  C    C      1     1000  10000   5000      0      0      0       C
END
"""
  pdb_str10="""
CRYST1   10.000   15.000   17.000  70.00 110.00 120.00 P 1
HETATM    1  C    C      1      -9.197  10.398  53.561  0.13  0.00           C
ANISOU    1  C    C      1     1000  10000   5000      0      0      0       C
END
"""
  pdb_str11="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       0.000   0.000   0.000  0.10  0.00           C
ANISOU    1  C    C      1     1000  10000   5000      0      0      0       C
HETATM    2  C    C      2       5.000   5.000   5.000  1.00  9.00           C
END
"""
  pdb_str12="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       5.000   5.000   5.000  1.00 20.00           C
END
"""

  found_two = 0
  for i, pdb_str in enumerate([pdb_str1, pdb_str2, pdb_str3, pdb_str4,
                               pdb_str5, pdb_str6, pdb_str7, pdb_str8,
                               pdb_str9, pdb_str10, pdb_str11, pdb_str12]):
    pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    ph.write_pdb_file(file_name = "m%d.pdb"%i,
      crystal_symmetry=xrs.crystal_symmetry())
    # To try more random adp values:
    #u_cart = ru(crystal_symmetry=xrs.crystal_symmetry(), u_scale=1, u_min=0.01)
    #xrs.set_u_cart(u_cart = u_cart)
    #
    fc = xrs.structure_factors(d_min = 0.5, algorithm = "direct").f_calc()
    fft_map = fc.fft_map(resolution_factor = 0.2)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    map_data = map_data.set_selected(map_data<4.0, 0)
    # Calculated anisotropy from ADP: this is exact answer
    a = list(xrs.scatterers().anisotropy(xrs.unit_cell()))[0]
    # Sphericity: main working implementation
    t0=time.time()
    s = maptbx.sphericity(map_data = map_data, unit_cell = xrs.unit_cell(),
      radius = 2.5, sites_frac = xrs.sites_frac())
    if(s.size()==2):
      assert approx_equal(s[1], 1.0, 0.001)
      found_two += 1
    assert approx_equal(a, s[0], 0.1)
    t1=time.time()
    print("#%2d: answer: %5.3f sphericity: %5.3f" % (i+1, a, s[0]))
  #
  assert found_two == 1

if __name__ == "__main__" :
  run()
