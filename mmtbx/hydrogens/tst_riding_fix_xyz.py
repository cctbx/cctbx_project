from __future__ import absolute_import, division, print_function

import time
import mmtbx.model
import iotbx.pdb
from libtbx import easy_run
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
import iotbx.pdb
from scitbx.array_family import flex


def exercise_00():
  """
  Test
  """
  # shake coordinates --> test does not depend on initial xyz
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    build_grm   = False,
    log         = null_out())
  xrs = model.get_xray_structure()
  xrs_shaken = xrs.deep_copy_scatterers()
  xrs_shaken.shake_sites_in_place(rms_difference=0.3)
  model.set_xray_structure(xrs_shaken)
  # save pdb with perturped coordinates
  f = open("tst_riding_fix_xyz_shaken.pdb", "w")
  f.write(model.model_as_pdb())
  f.close()
  # save original
  prefix = 'tst_riding_fix_xyz'
  f = open("%s.pdb" % prefix, "w")
  f.write(pdb_str)
  f.close()
  #
  # always use the same Rfree flags
  def generate_r_free_flags_systematic(miller_array):
    result = flex.bool()
    for i in range(miller_array.indices().size()):
      if(i%10==0): result.append(True)
      else: result.append(False)
    return miller_array.array(data = result)
  fobs_1 = abs(xrs.structure_factors(d_min=2.5).f_calc())
  flags_1 = generate_r_free_flags_systematic(miller_array=fobs_1)
  # Save mtz for refinement
  mtz = fobs_1.as_mtz_dataset(column_root_label="FP")
  mtz.add_miller_array(flags_1, column_root_label="R-free-flags")
  mtz.mtz_object().write(prefix+".mtz")
  # run phenix.refine
  selection_bool = "chain D"
  cmd = " ".join([
    "phenix.refine",
    "tst_riding_fix_xyz_shaken.pdb",
    "%s.mtz"%prefix,
    "refinement.main.number_of_macro_cycles=1",
    "refinement.refine.strategy=individual_sites",
    "refinement.refine.sites.individual='%s'" % selection_bool,
    "write_eff_file=False",
    "write_geo_file=False",
    "write_def_file=False",
    "write_model_cif_file=False",
    "write_map_coefficients=False",
    "--overwrite"
  ])
  print(cmd)
  easy_run.call(cmd)

  pdb_inp_refined = iotbx.pdb.input(file_name='tst_riding_fix_xyz_shaken_refine_001.pdb', source_info=None)
  model_refined = mmtbx.model.manager(
    model_input = pdb_inp_refined,
    log         = null_out())

  ph = model.get_hierarchy()
  ph_refined = model_refined.get_hierarchy()
  sele_bool = ph.atom_selection_cache().selection(
    string = selection_bool)
  sele_bool_refined = ph_refined.atom_selection_cache().selection(
    string = selection_bool)
  assert (sele_bool == sele_bool_refined)

  for a, ar, rbool in zip(ph.atoms(), ph_refined.atoms(), sele_bool):
    if not rbool:
      for ax, arx in zip(a.xyz, ar.xyz):
        #print(ax, arx)
        assert(approx_equal(ax, arx, eps=0.005))


pdb_str = """
CRYST1   18.717   24.019   18.868  90.00  90.00  90.00 P 1
SCALE1      0.053427  0.000000  0.000000        0.00000
SCALE2      0.000000  0.041634  0.000000        0.00000
SCALE3      0.000000  0.000000  0.053000        0.00000
ATOM      1  N   ARG A   2      14.830  19.124  14.055  1.00  0.00           N
ATOM      2  CA  ARG A   2      15.539  19.908  13.051  1.00  0.00           C
ATOM      3  C   ARG A   2      17.033  19.947  13.354  1.00  0.00           C
ATOM      4  O   ARG A   2      17.439  20.159  14.497  1.00  0.00           O
ATOM      5  CB  ARG A   2      14.977  21.330  12.986  1.00  0.00           C
ATOM      6  CG  ARG A   2      15.580  22.185  11.883  1.00  0.00           C
ATOM      7  CD  ARG A   2      14.908  23.547  11.806  1.00  0.00           C
ATOM      8  NE  ARG A   2      15.501  24.394  10.775  1.00  0.00           N
ATOM      9  CZ  ARG A   2      15.091  25.626  10.489  1.00  0.00           C
ATOM     10  NH1 ARG A   2      14.081  26.168  11.156  1.00  0.00           N
ATOM     11  NH2 ARG A   2      15.694  26.319   9.533  1.00  0.00           N
ATOM     12  HA  ARG A   2      15.418  19.495  12.182  1.00  0.00           H
ATOM     13  HB2 ARG A   2      14.020  21.280  12.832  1.00  0.00           H
ATOM     14  HB3 ARG A   2      15.150  21.772  13.832  1.00  0.00           H
ATOM     15  HG2 ARG A   2      16.523  22.323  12.063  1.00  0.00           H
ATOM     16  HG3 ARG A   2      15.462  21.739  11.030  1.00  0.00           H
ATOM     17  HD2 ARG A   2      13.969  23.426  11.595  1.00  0.00           H
ATOM     18  HD3 ARG A   2      15.004  23.998  12.659  1.00  0.00           H
ATOM     19  HE  ARG A   2      16.160  24.075  10.323  1.00  0.00           H
ATOM     20 HH11 ARG A   2      13.685  25.725  11.778  1.00  0.00           H
ATOM     21 HH12 ARG A   2      13.821  26.965  10.966  1.00  0.00           H
ATOM     22 HH21 ARG A   2      16.349  25.972   9.097  1.00  0.00           H
ATOM     23 HH22 ARG A   2      15.429  27.116   9.348  1.00  0.00           H
ATOM    283  N   TYR D   8      20.146  13.338   9.431  1.00  0.00           N
ATOM    284  CA  TYR D   8      20.988  14.478   9.066  1.00  0.00           C
ATOM    285  C   TYR D   8      22.308  14.337   9.818  1.00  0.00           C
ATOM    286  O   TYR D   8      22.402  14.703  10.993  1.00  0.00           O
ATOM    287  CB  TYR D   8      20.302  15.800   9.392  1.00  0.00           C
ATOM    288  CG  TYR D   8      19.076  16.075   8.551  1.00  0.00           C
ATOM    289  CD1 TYR D   8      19.192  16.579   7.263  1.00  0.00           C
ATOM    290  CD2 TYR D   8      17.801  15.832   9.046  1.00  0.00           C
ATOM    291  CE1 TYR D   8      18.075  16.832   6.490  1.00  0.00           C
ATOM    292  CE2 TYR D   8      16.678  16.082   8.282  1.00  0.00           C
ATOM    293  CZ  TYR D   8      16.820  16.582   7.005  1.00  0.00           C
ATOM    294  OH  TYR D   8      15.704  16.833   6.240  1.00  0.00           O
ATOM    295  H   TYR D   8      20.228  13.097  10.252  1.00  0.00           H
ATOM    296  HA  TYR D   8      21.173  14.454   8.115  1.00  0.00           H
ATOM    297  HB2 TYR D   8      20.027  15.788  10.322  1.00  0.00           H
ATOM    298  HB3 TYR D   8      20.931  16.523   9.244  1.00  0.00           H
ATOM    299  HD1 TYR D   8      20.037  16.749   6.913  1.00  0.00           H
ATOM    300  HD2 TYR D   8      17.703  15.494   9.907  1.00  0.00           H
ATOM    301  HE1 TYR D   8      18.168  17.170   5.629  1.00  0.00           H
ATOM    302  HE2 TYR D   8      15.831  15.914   8.626  1.00  0.00           H
ATOM    303  HH  TYR D   8      15.010  16.638   6.671  1.00  0.00           H
END
"""


#  cmd = " ".join([
#    "phenix.fmodel",
#    "%s.pdb"%prefix,
#    "high_res=2.5",
#    "type=real",
#    "label=f-obs",
#    "k_sol=0.4",
#    "b_sol=50",
#    "output.file_name=%s.mtz"%prefix,
#    "&> %s.log"%prefix
#  ])
#  print(cmd)
#  easy_run.call(cmd)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  print("OK. Time: %8.3f"%(time.time()-t0))
