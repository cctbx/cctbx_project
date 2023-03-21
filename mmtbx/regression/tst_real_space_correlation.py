
from __future__ import absolute_import, division, print_function
from mmtbx import real_space_correlation
import mmtbx.utils
from iotbx import file_reader
import iotbx.pdb
from scitbx.array_family import flex
import libtbx.load_env
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from libtbx import group_args
import os

def exercise_1():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1yjp_h.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1yjp.mtz",
    test=os.path.isfile)
  if (None in [pdb_file, mtz_file]):
    print("phenix_regression not found, skipping test")
    return False
  pdb_in = iotbx.pdb.input(pdb_file)
  hierarchy = pdb_in.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  xrs = pdb_in.xray_structure_simple()
  mtz_in = file_reader.any_file(mtz_file)
  f_obs = mtz_in.file_server.miller_arrays[0]
  r_free = mtz_in.file_server.miller_arrays[1]
  r_free = r_free.customized_copy(data=(r_free.data()==1))
  fmodel = mmtbx.utils.fmodel_simple(
    f_obs=f_obs,
    r_free_flags=r_free,
    xray_structures=[xrs],
    scattering_table="n_gaussian")
  map_stats = real_space_correlation.map_statistics_for_fragment(
    fragment=hierarchy,
    fmodel=fmodel)
  assert approx_equal(map_stats.cc, 0.960, eps=0.01)
  edm = fmodel.electron_density_map()
  map1_coeffs = edm.map_coefficients("2mFo-DFc")
  map1 = map1_coeffs.fft_map(
    resolution_factor=0.25).apply_sigma_scaling().real_map()
  map2_coeffs = edm.map_coefficients("Fmodel")
  map2 = map2_coeffs.fft_map(
    resolution_factor=0.25).apply_sigma_scaling().real_map()
  xray_structure = fmodel.xray_structure
  map_stats2 = real_space_correlation.map_statistics_for_atom_selection(
    atom_selection=flex.bool(xrs.sites_cart().size(), True),
    map1=map1,
    map2=map2,
    xray_structure=xrs)
  assert approx_equal(map_stats2.cc, map_stats.cc, 0.01)
  # XXX other code outside cctbx depends on the current API - do not simply
  # change the test if this breaks!
  overall_cc, results = real_space_correlation.simple(
    fmodel=fmodel,
    pdb_hierarchy=hierarchy,
    log=null_out())
  assert isinstance(overall_cc, float)
  assert isinstance(results, list)
  assert isinstance(results[0], group_args)
  assert (results[0].n_atoms == 1)
  assert (results[0].id_str == " A   GLY    1    N  ")
  return True

# test with synthetic twinned data
def exercise_2():
  pdb_str = """\
CRYST1   12.000    8.000   12.000  90.02  89.96  90.05 P 1           1
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
TER
ATOM      1  N   ASN B   6       1.414   5.113   6.019  1.00 12.99           N
ATOM      2  CA  ASN B   6       2.720   4.776   5.445  1.00 13.30           C
ATOM      3  C   ASN B   6       3.763   5.209   6.438  1.00 14.40           C
ATOM      4  O   ASN B   6       4.125   6.391   6.507  1.00 14.92           O
ATOM      5  CB  ASN B   6       2.922   5.513   4.131  1.00 13.13           C
ATOM      6  CG  ASN B   6       1.798   5.250   3.160  1.00 13.77           C
ATOM      7  OD1 ASN B   6       1.629   4.129   2.686  1.00 15.27           O
ATOM      8  ND2 ASN B   6       1.022   6.266   2.875  1.00 11.07           N
ATOM      9  N   TYR B   7       4.222   4.248   7.230  1.00 15.70           N
ATOM     10  CA  TYR B   7       5.113   4.552   8.370  1.00 16.18           C
ATOM     11  C   TYR B   7       6.547   4.754   7.929  1.00 16.91           C
ATOM     12  O   TYR B   7       6.964   4.259   6.878  1.00 16.76           O
ATOM     13  CB  TYR B   7       5.042   3.449   9.417  1.00 16.35           C
ATOM     14  CG  TYR B   7       3.659   3.296   9.977  1.00 15.45           C
ATOM     15  CD1 TYR B   7       2.756   2.398   9.402  1.00 16.68           C
ATOM     16  CD2 TYR B   7       3.224   4.098  11.023  1.00 15.80           C
ATOM     17  CE1 TYR B   7       1.476   2.267   9.896  1.00 14.46           C
ATOM     18  CE2 TYR B   7       1.929   3.975  11.545  1.00 15.33           C
ATOM     19  CZ  TYR B   7       1.063   3.065  10.959  1.00 16.09           C
ATOM     20  OH  TYR B   7      -0.207   2.910  11.443  1.00 15.39           O
ATOM     21  OXT TYR B   7       7.316   5.408   8.654  1.00 18.49           O
END
"""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  fc = fc.set_observation_type_xray_amplitude()
  sigf = flex.double(fc.size(), 0.1) + (fc.data() * 0.03)
  fc = fc.customized_copy(sigmas=sigf)
  # and now add twinning
  fc_twin = fc.twin_data(twin_law='-l,-k,-h', alpha=0.4)
  flags = fc_twin.generate_r_free_flags(use_lattice_symmetry=True)
  fmodel = mmtbx.utils.fmodel_simple(
    f_obs=fc_twin,
    r_free_flags=flags,
    xray_structures=[xrs],
    scattering_table="n_gaussian")
  assert (fmodel.twin_law is not None)
  map_stats = real_space_correlation.extract_map_stats_for_single_atoms(
    xray_structure=xrs,
    pdb_atoms=hierarchy.atoms(),
    fmodel=fmodel)
  sel_cache = hierarchy.atom_selection_cache()
  sel = sel_cache.selection("chain B")
  map_stats_2 = real_space_correlation.extract_map_stats_for_single_atoms(
    xray_structure=xrs,
    pdb_atoms=hierarchy.atoms(),
    fmodel=fmodel,
    selection=sel)
  map_stats_3 = real_space_correlation.map_statistics_for_atom_selection(
    atom_selection=sel,
    fmodel=fmodel)

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  print("OK")
