from __future__ import absolute_import, division, print_function
import time, os
from libtbx.test_utils import approx_equal
from iotbx import reflection_file_reader
from cctbx import maptbx
from mmtbx.regression.tst_polder import get_map, get_map_stats
from iotbx.data_manager import DataManager
from iotbx.cli_parser import run_program
from mmtbx.programs import fmodel, polder
from libtbx.utils import null_out

def exercise_00(prefix="tst_polder_box"):
  """
  Test for phenix.polder using the compute_box=True option.
  """
  dm = DataManager(['model'])
  dm.process_model_str('tst_polder_box.pdb', pdb_str)
  dm.process_model_str('tst_polder_box_ligand.pdb', pdb_str + pdb_str_ligand)
  dm.write_model_file(dm.get_model('tst_polder_box.pdb').model_as_pdb(),
                      filename  = 'tst_polder_box.pdb',
                      overwrite = True)
  dm.write_model_file(dm.get_model('tst_polder_box_ligand.pdb').model_as_pdb(),
                      filename  = 'tst_polder_box_ligand.pdb',
                      overwrite = True)

  # create test data with phenix.fmodel
  mtz_fn = "tst_polder_box.mtz"
  args = [
    'tst_polder_box_ligand.pdb',
    "high_res=2.0",
    "type=real",
    "label=f-obs",
    "k_sol=0.4",
    "b_sol=50",
    "output.file_name=%s" % mtz_fn]
  run_program(program_class=fmodel.Program, args=args, logger=null_out())

  selection = '((resseq 65 and name CG) or (resseq 7 and name OH) or \
    resseq 183 or resseq 95 or (resseq 148 and name CG))'

  args = [
    "tst_polder_box.pdb",
    "tst_polder_box.mtz",
    "compute_box=True",
    'solvent_exclusion_mask_selection="%s"' % selection]
    #'mask_output=True']

  r = run_program(program_class=polder.Program, args=args, logger=null_out())

  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    r.output_file).as_miller_arrays()
  mc_polder, mc_bias_omit, mc_omit = [None,]*3
  for ma in miller_arrays:
    lbl = ma.info().label_string()
    if(lbl == "mFo-DFc_polder,PHImFo-DFc_polder"):
      mc_polder = ma.deep_copy()
    if(lbl == "mFo-DFc_bias_omit,PHImFo-DFc_bias_omit"):
      mc_bias_omit = ma.deep_copy()
    if(lbl == "mFo-DFc_omit,PHImFo-DFc_omit"):
      mc_omit = ma.deep_copy()
  assert [mc_polder, mc_omit].count(None)==0
  cg = maptbx.crystal_gridding(
    unit_cell         = mc_polder.unit_cell(),
    d_min             = mc_polder.d_min(),
    resolution_factor = 0.25,
    space_group_info  = mc_polder.space_group_info())
  map_polder   = get_map(cg=cg, mc=mc_polder)
  map_omit     = get_map(cg=cg, mc=mc_omit)
  pdb_hierarchy = dm.get_model('tst_polder_box.pdb').get_hierarchy()
  sel = pdb_hierarchy.atom_selection_cache().selection(string = selection)
  sites_cart_lig = pdb_hierarchy.atoms().extract_xyz().select(sel)
  sites_frac_lig = mc_polder.unit_cell().fractionalize(sites_cart_lig)
  mp  = get_map_stats(map=map_polder,   sites_frac=sites_frac_lig)
  mo  = get_map_stats(map=map_omit,     sites_frac=sites_frac_lig)
  #
  mmm_mp = mp.min_max_mean().as_tuple()
  mmm_o = mo.min_max_mean().as_tuple()
  #print("Polder map : %7.3f %7.3f %7.3f"%mmm_mp)
  #print("Omit       : %7.3f %7.3f %7.3f"%mmm_o)
  #
  assert approx_equal(mmm_mp, [-2.051, 1.255, -0.162], eps=0.15)
  assert approx_equal(mmm_o,  [-0.558, 0.138, -0.068], eps=0.15)

  os.remove('tst_polder_box.pdb')
  os.remove('tst_polder_box_ligand.pdb')
  os.remove('tst_polder_box.mtz')
  os.remove(r.output_file)

# ---------------------------------------------------------------------------

pdb_str = """\
CRYST1   21.830   27.276   27.424  90.00  90.00  90.00 P 1
SCALE1      0.045809  0.000000  0.000000        0.00000
SCALE2      0.000000  0.036662  0.000000        0.00000
SCALE3      0.000000  0.000000  0.036464        0.00000
ATOM      1  N   TYR A   7      13.337  15.752  19.645  1.00 10.24           N
ATOM      2  CA  TYR A   7      12.308  16.751  19.902  1.00 11.36           C
ATOM      3  C   TYR A   7      12.565  17.364  21.244  1.00 12.55           C
ATOM      4  O   TYR A   7      12.979  16.662  22.161  1.00 14.36           O
ATOM      5  CB  TYR A   7      10.900  16.172  20.015  1.00 12.24           C
ATOM      6  CG  TYR A   7      10.374  15.482  18.784  1.00 13.20           C
ATOM      7  CD1 TYR A   7       9.910  16.256  17.734  1.00 15.61           C
ATOM      8  CD2 TYR A   7      10.326  14.099  18.719  1.00 13.25           C
ATOM      9  CE1 TYR A   7       9.393  15.632  16.618  1.00 16.74           C
ATOM     10  CE2 TYR A   7       9.809  13.467  17.600  1.00 14.24           C
ATOM     11  CZ  TYR A   7       9.346  14.248  16.562  1.00 16.67           C
ATOM     12  OH  TYR A   7       8.789  13.651  15.447  1.00 18.39           O
ATOM     13  N   HIS A  12       7.131  18.758  21.927  1.00 17.15           N
ATOM     14  CA  HIS A  12       6.572  17.487  21.517  1.00 15.83           C
ATOM     15  C   HIS A  12       7.214  16.460  22.424  1.00 16.05           C
ATOM     16  O   HIS A  12       8.426  16.251  22.372  1.00 15.95           O
ATOM     17  CB  HIS A  12       6.933  17.253  20.079  1.00 14.81           C
ATOM     18  CG  HIS A  12       6.260  16.046  19.452  1.00 16.05           C
ATOM     19  ND1 HIS A  12       5.199  16.054  18.655  1.00 17.40           N
ATOM     20  CD2 HIS A  12       6.719  14.757  19.509  1.00 16.62           C
ATOM     21  CE1 HIS A  12       5.000  14.839  18.206  1.00 17.35           C
ATOM     22  NE2 HIS A  12       5.919  14.078  18.726  1.00 19.34           N
ATOM     23  N   PHE A  43       7.914  21.178  11.433  1.00 16.32           N
ATOM     24  CA  PHE A  43       9.201  20.521  11.286  1.00 16.57           C
ATOM     25  C   PHE A  43      10.206  21.423  10.616  1.00 16.47           C
ATOM     26  O   PHE A  43       9.888  22.276   9.769  1.00 17.79           O
ATOM     27  CB  PHE A  43       9.088  19.272  10.437  1.00 17.43           C
ATOM     28  CG  PHE A  43       8.447  18.153  11.230  1.00 20.51           C
ATOM     29  CD1 PHE A  43       9.221  17.393  12.096  1.00 20.64           C
ATOM     30  CD2 PHE A  43       7.102  17.882  11.076  1.00 22.06           C
ATOM     31  CE1 PHE A  43       8.649  16.354  12.805  1.00 21.49           C
ATOM     32  CE2 PHE A  43       6.536  16.837  11.792  1.00 23.71           C
ATOM     33  CZ  PHE A  43       7.310  16.073  12.653  1.00 22.33           C
ATOM     34  N   ILE A  48      15.849  17.789   8.930  1.00 12.46           N
ATOM     35  CA  ILE A  48      15.469  16.420   8.550  1.00 13.43           C
ATOM     36  C   ILE A  48      16.295  16.013   7.331  1.00 13.28           C
ATOM     37  O   ILE A  48      16.830  14.898   7.337  1.00 12.47           O
ATOM     38  CB  ILE A  48      13.931  16.337   8.250  1.00 14.44           C
ATOM     39  CG1 ILE A  48      13.181  16.678   9.560  1.00 17.20           C
ATOM     40  CG2 ILE A  48      13.499  14.930   7.806  1.00 13.84           C
ATOM     41  CD1 ILE A  48      11.668  16.906   9.422  1.00 20.34           C
ATOM     42  N   GLN A  60      11.818   7.832   5.000  1.00 18.95           N
ATOM     43  CA  GLN A  60      11.755   8.086   6.430  1.00 20.21           C
ATOM     44  C   GLN A  60      11.086   6.941   7.175  1.00 19.68           C
ATOM     45  O   GLN A  60      10.875   7.000   8.389  1.00 19.28           O
ATOM     46  CB  GLN A  60      10.978   9.383   6.675  1.00 21.96           C
ATOM     47  CG  GLN A  60      11.647  10.627   6.126  1.00 25.34           C
ATOM     48  CD  GLN A  60      10.695  11.816   6.163  1.00 27.40           C
ATOM     49  OE1 GLN A  60       9.790  11.966   5.335  1.00 29.45           O
ATOM     50  NE2 GLN A  60      10.841  12.673   7.155  1.00 28.52           N
ATOM     51  N   LEU A  63      11.754   5.581  10.915  1.00 17.62           N
ATOM     52  CA  LEU A  63      12.627   6.435  11.685  1.00 17.71           C
ATOM     53  C   LEU A  63      11.929   7.021  12.904  1.00 15.54           C
ATOM     54  O   LEU A  63      10.864   7.652  12.785  1.00 17.03           O
ATOM     55  CB  LEU A  63      13.111   7.577  10.812  1.00 21.79           C
ATOM     56  CG ALEU A  63      14.302   8.243  10.167  0.50 25.94           C
ATOM     57  CD1ALEU A  63      13.988   9.692   9.787  0.50 26.52           C
ATOM     58  CD2ALEU A  63      15.437   8.172  11.170  0.50 26.38           C
ATOM     59  CG BLEU A  63      14.448   7.458  10.122  0.50 25.94           C
ATOM     60  CD1BLEU A  63      14.585   8.487   8.998  0.50 26.52           C
ATOM     61  CD2BLEU A  63      15.508   7.663  11.187  0.50 26.38           C
ATOM     62  N   THR A  64      12.504   6.810  14.068  1.00 12.15           N
ATOM     63  CA  THR A  64      11.993   7.397  15.285  1.00 12.07           C
ATOM     64  C   THR A  64      13.092   8.265  15.911  1.00 10.59           C
ATOM     65  O   THR A  64      14.263   8.236  15.515  1.00 10.38           O
ATOM     66  CB  THR A  64      11.573   6.306  16.277  1.00 12.33           C
ATOM     67  OG1 THR A  64      12.700   5.444  16.461  1.00 13.86           O
ATOM     68  CG2 THR A  64      10.370   5.516  15.788  1.00 12.38           C
ATOM     69  N   MET A  65      12.656   9.167  16.778  1.00 10.56           N
ATOM     70  CA  MET A  65      13.564  10.027  17.504  1.00  9.30           C
ATOM     71  C   MET A  65      13.794   9.445  18.893  1.00  9.44           C
ATOM     72  O   MET A  65      12.845   8.932  19.489  1.00 10.75           O
ATOM     73  CB  MET A  65      12.966  11.435  17.629  1.00 10.13           C
ATOM     74  CG  MET A  65      12.862  12.139  16.284  1.00 11.59           C
ATOM     75  SD  MET A  65      14.472  12.368  15.465  1.00 14.86           S
ATOM     76  CE  MET A  65      14.089  11.535  13.957  1.00 13.71           C
HETATM   89  O   HOH A  91       8.366  11.587  20.482  1.00 26.84           O
HETATM   90  O   HOH A  92       6.401  11.429  18.139  1.00 25.10           O
HETATM   91  O   HOH A  94      10.085   9.313  19.814  1.00 19.48           O
HETATM   92  O   HOH A  95       5.707   8.844  18.491  1.00 67.81           O
HETATM   93  O   HOH A 161       6.064  13.749  14.882  1.00 46.91           O
HETATM   94  O   HOH A 183       7.620   9.005   8.058  1.00 72.39           O
HETATM   95  O   HOH A 185       8.376   5.000  12.590  1.00 64.39           O
HETATM   96  O   HOH A 197       8.894   7.872  10.214  1.00 53.20           O
"""

pdb_str_ligand = """
HETATM   77  C2  MES A  88      11.178  10.844  10.521  1.00 66.12           C
HETATM   78  C3  MES A  88      10.035  10.700  11.264  1.00 65.44           C
HETATM   79  C5  MES A  88       9.898  12.901  11.698  1.00 65.50           C
HETATM   80  C6  MES A  88      11.005  13.173  10.925  1.00 66.33           C
HETATM   81  C7  MES A  88       8.949  11.458  13.133  1.00 63.09           C
HETATM   82  C8  MES A  88       8.983  10.302  13.874  1.00 61.16           C
HETATM   83  N4  MES A  88      10.025  11.657  12.279  1.00 64.72           N
HETATM   84  O1  MES A  88      11.134  12.150   9.913  1.00 66.57           O
HETATM   85  O1S MES A  88       7.392  11.377  15.767  1.00 59.34           O
HETATM   86  O2S MES A  88       8.043   8.588  15.664  1.00 59.97           O
HETATM   87  O3S MES A  88      10.061  10.480  16.373  1.00 59.16           O
HETATM   88  S   MES A  88       8.674  10.225  15.261  1.00 59.12           S
"""

# ---------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  print("OK. Time: %8.3f"%(time.time()-t0))
