from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import time
from libtbx.test_utils import approx_equal
import iotbx.pdb
from iotbx import reflection_file_reader
from cctbx import miller
from cctbx import maptbx
from scitbx.array_family import flex
from six.moves import range

pdb_str = """\
CRYST1   28.992   28.409   27.440  90.00  90.00  90.00 P 1
ATOM      1  N   ALA E   1       9.731  23.364   9.222  0.20 20.00           N
ATOM      2  CA  ALA E   1      10.928  22.678   9.693  0.20 20.00           C
ATOM      3  C   ALA E   1      10.619  21.229  10.055  0.20 20.00           C
ATOM      4  O   ALA E   1      11.301  20.629  10.886  0.20 20.00           O
ATOM      5  CB  ALA E   1      11.522  23.409  10.887  0.20 20.00           C
ATOM      6  N   HIS E   2       9.586  20.672   9.419  1.00 20.00           N
ATOM      7  CA  HIS E   2       9.202  19.291   9.695  1.00 20.00           C
ATOM      8  C   HIS E   2      10.295  18.321   9.264  1.00 20.00           C
ATOM      9  O   HIS E   2      10.653  17.402  10.010  1.00 20.00           O
ATOM     10  CB  HIS E   2       7.882  18.965   8.997  1.00 20.00           C
ATOM     11  CG  HIS E   2       6.738  19.828   9.431  1.00 20.00           C
ATOM     12  ND1 HIS E   2       5.915  19.498  10.485  1.00 20.00           N
ATOM     13  CD2 HIS E   2       6.282  21.010   8.953  1.00 20.00           C
ATOM     14  CE1 HIS E   2       5.000  20.438  10.638  1.00 20.00           C
ATOM     15  NE2 HIS E   2       5.200  21.367   9.721  1.00 20.00           N
ATOM     16  N   CYS E   3      10.837  18.510   8.058  1.00 20.00           N
ATOM     17  CA  CYS E   3      11.937  17.667   7.605  1.00 20.00           C
ATOM     18  C   CYS E   3      13.176  17.854   8.470  1.00 20.00           C
ATOM     19  O   CYS E   3      13.943  16.905   8.669  1.00 20.00           O
ATOM     20  CB  CYS E   3      12.260  17.966   6.141  1.00 20.00           C
ATOM     21  SG  CYS E   3      10.887  17.680   5.000  1.00 20.00           S
ATOM     22  N   ALA E   4      13.386  19.064   8.995  1.00 20.00           N
ATOM     23  CA  ALA E   4      14.520  19.295   9.885  1.00 20.00           C
ATOM     24  C   ALA E   4      14.353  18.539  11.196  1.00 20.00           C
ATOM     25  O   ALA E   4      15.308  17.933  11.696  1.00 20.00           O
ATOM     26  CB  ALA E   4      14.689  20.792  10.142  1.00 20.00           C
TER
ATOM     98  N   PHE A   1       9.174   9.310  14.969  0.20 20.00           N
ATOM     99  CA  PHE A   1      10.235   8.421  14.508  0.20 20.00           C
ATOM    100  C   PHE A   1      11.171   8.048  15.652  0.20 20.00           C
ATOM    101  O   PHE A   1      12.391   8.032  15.489  0.20 20.00           O
ATOM    102  CB  PHE A   1      11.025   9.073  13.371  0.20 20.00           C
ATOM    103  CG  PHE A   1      10.196   9.391  12.160  0.20 20.00           C
ATOM    104  CD1 PHE A   1      10.002   8.444  11.168  0.20 20.00           C
ATOM    105  CD2 PHE A   1       9.611  10.638  12.012  0.20 20.00           C
ATOM    106  CE1 PHE A   1       9.240   8.734  10.052  0.20 20.00           C
ATOM    107  CE2 PHE A   1       8.847  10.934  10.898  0.20 20.00           C
ATOM    108  CZ  PHE A   1       8.662   9.980   9.917  0.20 20.00           C
TER
END
"""

def get_map(cg, mc):
  fft_map = miller.fft_map(
    crystal_gridding     = cg,
    fourier_coefficients = mc)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def get_map_stats(map, sites_frac):
  map_values = flex.double()
  for sf in sites_frac:
    map_values.append(map.eight_point_interpolation(sf))
  return map_values

def format_map_stat(m):
  return m.min_max_mean().as_tuple(), (m>flex.mean(m)).count(True)

# Test reading mtz with anomalous arrays for F and Rfree
def exercise_01(fobs_1, flags_1, prefix):
  mtz = fobs_1.as_mtz_dataset(column_root_label="FP1")
  mtz.add_miller_array(flags_1, column_root_label="R-free-flags")
  mtz.mtz_object().write(prefix+".mtz")
  #
  selection = "chain A"
  cmd = " ".join([
    "phenix.polder",
    "tst_polder_2.pdb",
    "%s.mtz" % prefix,
    "sphere_radius=3",
    "output_file_name_prefix=%s" % prefix,
    'solvent_exclusion_mask_selection="%s"' % selection,
    "&> %s.log" % prefix
  ])
  print(cmd)
  easy_run.call(cmd)
  # exact values: 11.129  16.152  13.006
  check(
    tuple_calc = [11.1, 16.2, 13.0],
    selection  = selection,
    prefix     = prefix+'_')

# Test reading mtz with anomalous array for F (no Rfree present)
def exercise_02(fobs_1, prefix):
  mtz = fobs_1.as_mtz_dataset(column_root_label="FP1")
  mtz.mtz_object().write(prefix+".mtz")
  #
  selection = "chain E and resseq 1"
  cmd = " ".join([
    "phenix.polder",
    "tst_polder_2.pdb",
    "%s.mtz" % prefix,
    "sphere_radius=3",
    "output_file_name_prefix=%s" % prefix,
    'solvent_exclusion_mask_selection="%s"' % selection,
    "&> %s.log" % prefix
  ])
  print(cmd)
  easy_run.call(cmd)
  check(
    tuple_calc = [12.7, 22.7, 17.7],
    selection  = selection,
    prefix     = prefix+'_')

# Mixture: mtz with anomalous F and "normal" Rfree
def exercise_03(fobs_1, flags_1, prefix):
  mtz = fobs_1.as_mtz_dataset(column_root_label="FP1")
  mtz.add_miller_array(flags_1, column_root_label="R-free-flags")
  mtz.mtz_object().write(prefix+".mtz")
  #
  selection = "chain A"
  cmd = " ".join([
    "phenix.polder",
    "tst_polder_2.pdb",
    "%s.mtz" % prefix,
    "sphere_radius=3",
    "output_file_name_prefix=%s" % prefix,
    'solvent_exclusion_mask_selection="%s"' % selection,
    "&> %s.log" % prefix
  ])
  print(cmd)
  easy_run.call(cmd)
  # exact values: 11.129  16.152  13.006
  check(
    tuple_calc = [11.1, 16.2, 13.0],
    selection  = selection,
    prefix     = prefix+'_')

def exercise():
  """
  Test for phenix.polder: accepting anomalous data labels.
  """
  f = open("tst_polder_2.pdb", "w")
  f.write(pdb_str)
  f.close()
  #
  pdb_in = iotbx.pdb.input(source_info=None,lines=pdb_str)
  xrs = pdb_in.xray_structure_simple()
  pdb_hierarchy = iotbx.pdb.input(
    source_info=None, lines=pdb_str).construct_hierarchy()
  # map values are not reproducible if not systematic Rfree!
  def generate_r_free_flags_systematic(miller_array):
    result = flex.bool()
    for i in range(miller_array.indices().size()):
      if(i%10==0): result.append(True)
      else: result.append(False)
    return miller_array.array(data = result)
  # generate anomalous and non-anomalous data arrays
  f_anom = abs(xrs.structure_factors(d_min=3.0).f_calc().generate_bijvoet_mates())
  f_obs = abs(xrs.structure_factors(d_min=3.0).f_calc())
  flags_obs = generate_r_free_flags_systematic(miller_array=f_obs)
  flags_anom = flags_obs.generate_bijvoet_mates()
  #
  print('*'*79)
  print('Test reading one mtz with F and Rfree as anomalous arrays')
  exercise_01(
    fobs_1  = f_anom,
    flags_1 = flags_anom,
    prefix  = "tst_polder_2_1")
  print("OK")
  print('*'*79)
  print('Test reading one mtz with F as anomalous array (no Rfree present)')
  exercise_02(
    fobs_1 = f_anom,
    prefix = "tst_polder_2_2")
  print("OK")
  print('*'*79)
  print('Test reading mtz with F as anomalous array and Rfree is usual array')
  exercise_03(
    fobs_1  = f_anom,
    flags_1 = flags_obs,
    prefix  = "tst_polder_2_3")
  print("OK")

def check(tuple_calc, selection, prefix):
  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    prefix+"polder_map_coeffs.mtz").as_miller_arrays()
  mc_polder = None
  for ma in miller_arrays:
    lbl = ma.info().label_string()
    if(lbl == "mFo-DFc_polder,PHImFo-DFc_polder"):
      mc_polder = ma.deep_copy()
  assert (mc_polder is not None)
  cg = maptbx.crystal_gridding(
    unit_cell         = mc_polder.unit_cell(),
    d_min             = mc_polder.d_min(),
    resolution_factor = 0.25,
    space_group_info  = mc_polder.space_group_info())
  map_polder   = get_map(cg=cg, mc=mc_polder)
  pdb_hierarchy = iotbx.pdb.input(
    source_info=None, lines=pdb_str).construct_hierarchy()
  sel = pdb_hierarchy.atom_selection_cache().selection(string = selection)
  sites_cart_lig = pdb_hierarchy.atoms().extract_xyz().select(sel)
  sites_frac_lig = mc_polder.unit_cell().fractionalize(sites_cart_lig)
  mp  = get_map_stats(
    map        = map_polder,
    sites_frac = sites_frac_lig)
  #
  mmm_mp = mp.min_max_mean().as_tuple()
  print("Polder map : %7.3f %7.3f %7.3f" % mmm_mp)
  assert approx_equal(mmm_mp, tuple_calc, eps=1.0), "\
   calculated is %s and expected is %s" % (mmm_mp, tuple_calc)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f"%(time.time()-t0))
