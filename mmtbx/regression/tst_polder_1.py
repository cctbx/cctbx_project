from __future__ import absolute_import, division, print_function
import time, os
from libtbx.test_utils import approx_equal
from iotbx.data_manager import DataManager
from iotbx.cli_parser import run_program
import iotbx.pdb
from mmtbx.programs import polder
from iotbx import reflection_file_reader
from cctbx import maptbx
from scitbx.array_family import flex
from six.moves import range
from libtbx.utils import null_out
from mmtbx.regression.tst_polder import get_map, get_map_stats


def generate_r_free_flags_systematic(miller_array):
  result = flex.bool()
  for i in range(miller_array.indices().size()):
    if(i%10==0): result.append(True)
    else: result.append(False)
  return miller_array.array(data = result)

# ---------------------------------------------------------------------------

def exercise_01(fobs_1, fobs_2, flags_1, flags_2):
  '''
  Test reading one mtz with different choices for F and Rfree
  '''
  mtz = fobs_1.as_mtz_dataset(column_root_label="FP1")
  mtz.add_miller_array(fobs_2, column_root_label="FP2")
  mtz.add_miller_array(flags_1, column_root_label="R-free-flags-1")
  mtz.add_miller_array(flags_2, column_root_label="R-free-flags-2")
  mtz.mtz_object().write("tst_polder_1_1.mtz")
  #
  selection = "chain A"

  args = [
    "tst_polder_1.pdb",
    "tst_polder_1_1.mtz",
    "sphere_radius=3",
    'solvent_exclusion_mask_selection="%s"' % selection,
    "user_selected_labels=FP1",
    "user_selected_labels=R-free-flags-1",
    "output_file_name_prefix=tst_polder_1_1"]

  r = run_program(program_class=polder.Program, args=args, logger=null_out())

  check(
    tuple_calc=[10.6, 16.0, 13.5],
    selection=selection,
    filename = r.output_file)

  os.remove("tst_polder_1_1.mtz")
  os.remove(r.output_file)

# ---------------------------------------------------------------------------

def exercise_02(fobs_1, flags_1):
  '''
  Test reading two cns files: one with Fobs, one with Rfree;
  without specifying labels --> test automatic machinery
  '''
  with open("tst_polder_1_2.cns", "w") as file_object:
    fobs_1.export_as_cns_hkl(
      file_object=file_object,
      file_name="tst_polder_1_2.cns",
      info=["calculated structure factors FP1"])
  #
  with open("tst_polder_1_2_free.cns", "w") as file_object:
    flags_1.export_as_cns_hkl(
      file_object=file_object,
      file_name="tst_polder_1_2_free.cns",
      info=["R-free-flags for FP1"])
  #
  selection = "chain E and resseq 14"
  args = [
    "tst_polder_1.pdb",
    "tst_polder_1_2.cns",
    "tst_polder_1_2_free.cns",
    "sphere_radius=5",
    'solvent_exclusion_mask_selection="%s"' % selection,
    "output_file_name_prefix=tst_polder_1_2"]

  r = run_program(program_class=polder.Program, args=args, logger=null_out())

  check(
    tuple_calc=[10.7, 15.1, 13.4],
    selection=selection,
    filename = r.output_file)

  os.remove("tst_polder_1_2.cns")
  os.remove("tst_polder_1_2_free.cns")
  os.remove(r.output_file)

# ---------------------------------------------------------------------------

def exercise_03(fobs_2, flags_2):
  '''
  Test reading model as cif, Fobs from mtz and Rfree from cns file
  '''
  mtz = fobs_2.as_mtz_dataset(column_root_label="FP1")
  mtz.mtz_object().write("tst_polder_1_3.mtz")
  with open("tst_polder_1_3.cns", "w") as file_object:
    flags_2.export_as_cns_hkl(
      file_object=file_object,
      file_name="tst_polder_1_3.cns",
      info=["R-free-flags for FP1"])
  #
  selection = "chain E and resseq 7"
  args = [
    "tst_polder_1.cif",
    "tst_polder_1_3.mtz",
    "tst_polder_1_3.cns",
    "sphere_radius=5",
    'solvent_exclusion_mask_selection="%s"' % selection,
    "output_file_name_prefix=tst_polder_1_3"]

  r = run_program(program_class=polder.Program, args=args, logger=null_out())

  check(
    tuple_calc=[13.7, 20.3, 18.4],
    selection=selection,
    filename = r.output_file)

  # Clean up files
  os.remove('tst_polder_1_3.mtz')
  os.remove('tst_polder_1_3.cns')
  os.remove(r.output_file)

# ---------------------------------------------------------------------------

def exercise():
  """
  Test for phenix.polder: Reading files.
  """
  dm = DataManager(['model'])
  dm.process_model_str("tst_polder_1.pdb", pdb_str)
  model = dm.get_model()
  dm.write_model_file(model.model_as_pdb(),
                      filename  = "tst_polder_1.pdb",
                      overwrite = True)
  dm.write_model_file(model.model_as_mmcif(),
                      filename  = "tst_polder_1.cif",
                      overwrite = True)
  xrs = model.get_xray_structure()

  fobs_1 = abs(xrs.structure_factors(d_min=3.0).f_calc())
  fobs_2 = abs(xrs.structure_factors(d_min=2.5).f_calc())
  flags_1 = generate_r_free_flags_systematic(miller_array=fobs_1)
  flags_2 = generate_r_free_flags_systematic(miller_array=fobs_2)

  exercise_01(fobs_1, fobs_2, flags_1, flags_2)
  exercise_02(fobs_1, flags_1)
  exercise_03(fobs_2, flags_2)

  # Clean up files
  os.remove('tst_polder_1.pdb')
  os.remove('tst_polder_1.cif')

# ---------------------------------------------------------------------------

def check(tuple_calc, selection, filename):
  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    filename).as_miller_arrays()
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
  mp  = get_map_stats(map=map_polder,   sites_frac=sites_frac_lig)
  #
  mmm_mp = mp.min_max_mean().as_tuple()
  #print("Polder map : %7.3f %7.3f %7.3f" % mmm_mp)
  assert approx_equal(mmm_mp, tuple_calc, eps=1.0), "\
   calculated is %s and expected is %s" % (mmm_mp, tuple_calc)

# ---------------------------------------------------------------------------

pdb_str = """
CRYST1   28.992   28.409   27.440  90.00  90.00  90.00 P 1
ATOM      1  N   ALA E   1       9.731  23.364   9.222  1.00 20.00           N
ATOM      2  CA  ALA E   1      10.928  22.678   9.693  1.00 20.00           C
ATOM      3  C   ALA E   1      10.619  21.229  10.055  1.00 20.00           C
ATOM      4  O   ALA E   1      11.301  20.629  10.886  1.00 20.00           O
ATOM      5  CB  ALA E   1      11.522  23.409  10.887  1.00 20.00           C
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
ATOM     27  N   ILE E   5      13.146  18.559  11.766  1.00 20.00           N
ATOM     28  CA  ILE E   5      12.893  17.811  12.993  1.00 20.00           C
ATOM     29  C   ILE E   5      13.003  16.313  12.735  1.00 20.00           C
ATOM     30  O   ILE E   5      13.481  15.555  13.588  1.00 20.00           O
ATOM     31  CB  ILE E   5      11.519  18.193  13.576  1.00 20.00           C
ATOM     32  CG1 ILE E   5      11.493  19.677  13.949  1.00 20.00           C
ATOM     33  CG2 ILE E   5      11.187  17.338  14.789  1.00 20.00           C
ATOM     34  CD1 ILE E   5      10.171  20.139  14.522  1.00 20.00           C
ATOM     35  N   TYR E   6      12.577  15.864  11.552  1.00 20.00           N
ATOM     36  CA  TYR E   6      12.688  14.447  11.218  1.00 20.00           C
ATOM     37  C   TYR E   6      14.146  14.025  11.072  1.00 20.00           C
ATOM     38  O   TYR E   6      14.524  12.925  11.492  1.00 20.00           O
ATOM     39  CB  TYR E   6      11.910  14.145   9.937  1.00 20.00           C
ATOM     40  CG  TYR E   6      10.412  14.302  10.074  1.00 20.00           C
ATOM     41  CD1 TYR E   6       9.796  14.221  11.316  1.00 20.00           C
ATOM     42  CD2 TYR E   6       9.614  14.530   8.961  1.00 20.00           C
ATOM     43  CE1 TYR E   6       8.427  14.363  11.446  1.00 20.00           C
ATOM     44  CE2 TYR E   6       8.244  14.674   9.080  1.00 20.00           C
ATOM     45  CZ  TYR E   6       7.656  14.589  10.325  1.00 20.00           C
ATOM     46  OH  TYR E   6       6.293  14.731  10.449  1.00 20.00           O
ATOM     47  N   THR E   7      14.979  14.884  10.479  0.20 20.00           N
ATOM     48  CA  THR E   7      16.391  14.547  10.326  0.20 20.00           C
ATOM     49  C   THR E   7      17.120  14.587  11.664  0.20 20.00           C
ATOM     50  O   THR E   7      17.960  13.724  11.944  0.20 20.00           O
ATOM     51  CB  THR E   7      17.058  15.492   9.326  0.20 20.00           C
ATOM     52  OG1 THR E   7      16.757  16.850   9.672  0.20 20.00           O
ATOM     53  CG2 THR E   7      16.571  15.207   7.912  0.20 20.00           C
ATOM     54  N   ILE E   8      16.814  15.580  12.502  1.00 20.00           N
ATOM     55  CA  ILE E   8      17.450  15.665  13.815  1.00 20.00           C
ATOM     56  C   ILE E   8      17.050  14.475  14.679  1.00 20.00           C
ATOM     57  O   ILE E   8      17.893  13.855  15.340  1.00 20.00           O
ATOM     58  CB  ILE E   8      17.103  17.004  14.492  1.00 20.00           C
ATOM     59  CG1 ILE E   8      17.763  18.165  13.745  1.00 20.00           C
ATOM     60  CG2 ILE E   8      17.530  16.996  15.951  1.00 20.00           C
ATOM     61  CD1 ILE E   8      17.464  19.524  14.342  1.00 20.00           C
ATOM     62  N   HIS E   9      15.760  14.128  14.678  1.00 20.00           N
ATOM     63  CA  HIS E   9      15.302  12.994  15.472  1.00 20.00           C
ATOM     64  C   HIS E   9      15.787  11.669  14.897  1.00 20.00           C
ATOM     65  O   HIS E   9      15.926  10.689  15.637  1.00 20.00           O
ATOM     66  CB  HIS E   9      13.776  13.002  15.572  1.00 20.00           C
ATOM     67  CG  HIS E   9      13.232  14.116  16.411  1.00 20.00           C
ATOM     68  ND1 HIS E   9      14.038  15.058  17.012  1.00 20.00           N
ATOM     69  CD2 HIS E   9      11.961  14.436  16.753  1.00 20.00           C
ATOM     70  CE1 HIS E   9      13.288  15.912  17.685  1.00 20.00           C
ATOM     71  NE2 HIS E   9      12.024  15.557  17.545  1.00 20.00           N
ATOM     72  N   SER E  10      16.046  11.615  13.589  1.00 20.00           N
ATOM     73  CA  SER E  10      16.560  10.388  12.991  1.00 20.00           C
ATOM     74  C   SER E  10      18.038  10.191  13.305  1.00 20.00           C
ATOM     75  O   SER E  10      18.458   9.080  13.651  1.00 20.00           O
ATOM     76  CB  SER E  10      16.334  10.403  11.479  1.00 20.00           C
ATOM     77  OG  SER E  10      17.100  11.420  10.858  1.00 20.00           O
ATOM     78  N   VAL E  11      18.839  11.252  13.191  1.00 20.00           N
ATOM     79  CA  VAL E  11      20.262  11.143  13.498  1.00 20.00           C
ATOM     80  C   VAL E  11      20.471  10.920  14.991  1.00 20.00           C
ATOM     81  O   VAL E  11      21.273  10.071  15.398  1.00 20.00           O
ATOM     82  CB  VAL E  11      21.014  12.389  12.998  1.00 20.00           C
ATOM     83  CG1 VAL E  11      22.469  12.345  13.440  1.00 20.00           C
ATOM     84  CG2 VAL E  11      20.920  12.491  11.483  1.00 20.00           C
ATOM     85  N   ASP E  12      19.755  11.672  15.830  1.00 20.00           N
ATOM     86  CA  ASP E  12      19.881  11.484  17.272  1.00 20.00           C
ATOM     87  C   ASP E  12      19.303  10.148  17.719  1.00 20.00           C
ATOM     88  O   ASP E  12      19.822   9.534  18.658  1.00 20.00           O
ATOM     89  CB  ASP E  12      19.200  12.633  18.016  1.00 20.00           C
ATOM     90  CG  ASP E  12      19.888  13.964  17.788  1.00 20.00           C
ATOM     91  OD1 ASP E  12      21.099  13.964  17.481  1.00 20.00           O
ATOM     92  OD2 ASP E  12      19.220  15.010  17.921  1.00 20.00           O
ATOM     93  N   ALA E  13      18.236   9.683  17.066  1.00 20.00           N
ATOM     94  CA  ALA E  13      17.640   8.405  17.442  1.00 20.00           C
ATOM     95  C   ALA E  13      18.531   7.237  17.035  1.00 20.00           C
ATOM     96  O   ALA E  13      18.671   6.266  17.787  1.00 20.00           O
ATOM     97  CB  ALA E  13      16.253   8.268  16.816  1.00 20.00           C
ATOM     98  N   PHE E  14      19.140   7.312  15.850  0.20 20.00           N
ATOM     99  CA  PHE E  14      20.024   6.243  15.402  0.20 20.00           C
ATOM    100  C   PHE E  14      21.374   6.283  16.108  0.20 20.00           C
ATOM    101  O   PHE E  14      22.022   5.241  16.254  0.20 20.00           O
ATOM    102  CB  PHE E  14      20.223   6.320  13.887  0.20 20.00           C
ATOM    103  CG  PHE E  14      18.980   6.023  13.096  0.20 20.00           C
ATOM    104  CD1 PHE E  14      17.953   5.275  13.646  0.20 20.00           C
ATOM    105  CD2 PHE E  14      18.841   6.492  11.800  0.20 20.00           C
ATOM    106  CE1 PHE E  14      16.809   5.000  12.919  0.20 20.00           C
ATOM    107  CE2 PHE E  14      17.700   6.221  11.068  0.20 20.00           C
ATOM    108  CZ  PHE E  14      16.683   5.474  11.629  0.20 20.00           C
ATOM    109  N   ALA E  15      21.810   7.463  16.548  1.00 20.00           N
ATOM    110  CA  ALA E  15      23.084   7.598  17.240  1.00 20.00           C
ATOM    111  C   ALA E  15      22.974   7.379  18.742  1.00 20.00           C
ATOM    112  O   ALA E  15      23.992   7.108  19.390  1.00 20.00           O
ATOM    113  CB  ALA E  15      23.689   8.980  16.973  1.00 20.00           C
ATOM    114  N   GLU E  16      21.776   7.489  19.307  1.00 20.00           N
ATOM    115  CA  GLU E  16      21.584   7.294  20.739  1.00 20.00           C
ATOM    116  C   GLU E  16      20.666   6.107  21.011  1.00 20.00           C
ATOM    117  O   GLU E  16      20.973   5.248  21.838  1.00 20.00           O
ATOM    118  CB  GLU E  16      21.012   8.559  21.382  1.00 20.00           C
ATOM    119  CG  GLU E  16      21.901   9.784  21.238  1.00 20.00           C
ATOM    120  CD  GLU E  16      21.296  11.020  21.875  1.00 20.00           C
ATOM    121  OE1 GLU E  16      20.188  10.917  22.440  1.00 20.00           O
ATOM    122  OE2 GLU E  16      21.930  12.095  21.809  1.00 20.00           O
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

# ---------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("OK. Time: %8.3f"%(time.time()-t0))
