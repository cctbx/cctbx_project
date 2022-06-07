from __future__ import absolute_import, division, print_function
from mmtbx import utils
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import format_cpu_times, null_out, Sorry
import sys
import mmtbx.model

def exercise_d_data_target_d_atomic_params():
  import iotbx.pdb
  import mmtbx.f_model
  from scitbx.array_family import flex
  good = """
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 212121
HETATM  115  O   HOH A  18       3.000   5.000   5.000  1.00 10.00           O
HETATM  115  O   HOH A  18       5.000   5.000   8.000  1.00 10.00           O
TER
END
  """
  bad = """
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 212121
HETATM  115  O   HOH A  18       3.000   5.000   5.000  0.50 10.00           O
HETATM  115  O   HOH A  19       5.000   5.000   8.000  0.90 10.00           O
TER
END
  """
  for update_f_part1 in [True, False]:
    pdb_inp = iotbx.pdb.input(source_info=None, lines=bad)
    xrs = pdb_inp.xray_structure_simple()
    xrs.scattering_type_registry(table = "wk1995")
    #
    xb = iotbx.pdb.input(source_info=None, lines=good).xray_structure_simple()
    xb.scattering_type_registry(table = "wk1995")
    f_obs = abs(xb.structure_factors(d_min=1,
      algorithm="direct").f_calc()).set_observation_type_xray_amplitude()
    #
    sfp = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    sfp.algorithm="direct"
    for target_name in ["ls_wunit_kunit", "ls_wunit_k1", "ml"]:
      print("%s"%target_name, "-"*50)
      fmodel = mmtbx.f_model.manager(
        f_obs                        = f_obs,
        xray_structure               = xrs.deep_copy_scatterers(),
        target_name                  = target_name,
        sf_and_grads_accuracy_params = sfp)
      fmodel.update_all_scales(update_f_part1=update_f_part1)
      alpha_beta = fmodel.alpha_beta()
      print("R-work: %6.4f"%fmodel.r_work())
      #
      tg = mmtbx.utils.experimental_data_target_and_gradients(fmodel = fmodel,
        alpha_beta=alpha_beta)
      tg.show()
      eps = 1.e-6
      # occupancies
      go = tg.grad_occ()
      for i in [0,1]:
        xrs1 = fmodel.xray_structure.deep_copy_scatterers()
        xrs2 = fmodel.xray_structure.deep_copy_scatterers()
        #
        xrs1.scatterers()[i].occupancy+=+eps
        tg.update_xray_structure(xray_structure=xrs1, alpha_beta=alpha_beta)
        t1 = tg.target()
        #
        xrs2.scatterers()[i].occupancy+=-eps
        tg.update_xray_structure(xray_structure=xrs2, alpha_beta=alpha_beta)
        t2 = tg.target()
        #
        gfd = (t1-t2)/(2*eps)
        print(gfd, go[i])
        assert approx_equal(go[i], gfd, 1.e-5)
      # sites_cart
      gc = tg.grad_sites_cart()
      uc = fmodel.xray_structure.unit_cell()
      for i in [0,1]:
        gfd = []
        for e in [(eps,0,0),(0,eps,0),(0,0,eps)]:
          xrs1 = fmodel.xray_structure.deep_copy_scatterers()
          xrs2 = fmodel.xray_structure.deep_copy_scatterers()
          #
          s1 = flex.vec3_double([uc.orthogonalize(xrs1.scatterers()[i].site)])
          s1 = s1+flex.vec3_double([e])
          xrs1.scatterers()[i].site = uc.fractionalize(s1[0])
          tg.update_xray_structure(xray_structure=xrs1, alpha_beta=alpha_beta)
          t1 = tg.target()
          #
          s2 = flex.vec3_double([uc.orthogonalize(xrs2.scatterers()[i].site)])
          s2 = s2-flex.vec3_double([e])
          xrs2.scatterers()[i].site = uc.fractionalize(s2[0])
          tg.update_xray_structure(xray_structure=xrs2, alpha_beta=alpha_beta)
          t2 = tg.target()
          #
          gfd.append( (t1-t2)/(2*eps) )
        print(gfd, list(gc[i]))
        assert approx_equal(gc[i], gfd, 1.e-5)

def exercise_d_data_target_d_atomic_params2():
  import iotbx.pdb
  import mmtbx.f_model
  from scitbx.array_family import flex
  good = """
CRYST1   26.771   27.605   16.145  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1      21.771  22.605   5.434  1.00 10.00           N
ATOM      2  CA  ASP A   1      20.399  22.215   5.000  1.00 10.00           C
ATOM      3  C   ASP A   1      19.546  21.744   6.168  1.00 10.00           C
ATOM      4  O   ASP A   1      19.997  20.958   7.001  1.00 10.00           O
ATOM      5  N   VAL A   2      18.313  22.233   6.229  1.00 10.00           N
ATOM      6  CA  VAL A   2      17.402  21.833   7.287  1.00 10.00           C
ATOM      7  C   VAL A   2      16.896  20.436   6.954  1.00 10.00           C
ATOM      8  O   VAL A   2      16.413  20.188   5.850  1.00 10.00           O
ATOM      9  N   GLN A   3      17.009  19.531   7.918  1.00 10.00           N
ATOM     10  CA  GLN A   3      16.590  18.147   7.740  1.00 10.00           C
ATOM     11  C   GLN A   3      15.154  17.904   8.207  1.00 10.00           C
ATOM     12  O   GLN A   3      14.813  18.175   9.356  1.00 10.00           O
ATOM     13  N   MET A   4      14.318  17.383   7.310  1.00 10.00           N
ATOM     14  CA  MET A   4      12.921  17.091   7.630  1.00 10.00           C
ATOM     15  C   MET A   4      12.750  15.585   7.849  1.00 10.00           C
ATOM     16  O   MET A   4      13.057  14.784   6.965  1.00 10.00           O
TER
ATOM     17  N   THR B   5      14.260  17.201  11.025  1.00 10.00           N
ATOM     18  CA  THR B   5      14.076  15.783  11.355  1.00 10.00           C
ATOM     19  C   THR B   5      12.612  15.405  11.550  1.00 10.00           C
ATOM     20  O   THR B   5      11.958  15.902  12.463  1.00 10.00           O
ATOM     21  N   GLN B   6      12.104  14.518  10.697  1.00 10.00           N
ATOM     22  CA  GLN B   6      10.712  14.084  10.797  1.00 10.00           C
ATOM     23  C   GLN B   6      10.571  12.704  11.424  1.00 10.00           C
ATOM     24  O   GLN B   6      11.336  11.787  11.118  1.00 10.00           O
ATOM     25  N   THR B   7       9.565  12.570  12.283  1.00 10.00           N
ATOM     26  CA  THR B   7       9.266  11.322  12.973  1.00 10.00           C
ATOM     27  C   THR B   7       7.756  11.127  12.915  1.00 10.00           C
ATOM     28  O   THR B   7       7.000  12.070  13.145  1.00 10.00           O
ATOM     29  N   PRO B   8       7.291   9.907  12.609  1.00 10.00           N
ATOM     30  CA  PRO B   8       8.060   8.698  12.310  1.00 10.00           C
ATOM     31  C   PRO B   8       8.406   8.619  10.827  1.00 10.00           C
ATOM     32  O   PRO B   8       8.058   9.514  10.054  1.00 10.00           O
ATOM     33  N   LEU B   9       9.093   7.553  10.427  1.00 10.00           N
ATOM     34  CA  LEU B   9       9.459   7.385   9.023  1.00 10.00           C
ATOM     35  C   LEU B   9       8.232   7.000   8.213  1.00 10.00           C
ATOM     36  O   LEU B   9       8.026   7.506   7.113  1.00 10.00           O
TER
  """
  bad = """
CRYST1   26.771   27.605   16.145  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1      21.771  22.605   5.434  1.00 10.00           N
ATOM      2  CA  ASP A   1      20.399  22.215   5.000  1.00 10.00           C
ATOM      3  C   ASP A   1      19.546  21.744   6.168  1.00 10.00           C
ATOM      4  O   ASP A   1      19.997  20.958   7.001  1.00 10.00           O
ATOM      5  N   VAL A   2      18.313  22.233   6.229  1.00 10.00           N
ATOM      6  CA  VAL A   2      17.402  21.833   7.287  1.00 10.00           C
ATOM      7  C   VAL A   2      16.896  20.436   6.954  1.00 10.00           C
ATOM      8  O   VAL A   2      16.413  20.188   5.850  1.00 10.00           O
ATOM      9  N   GLN A   3      17.009  19.531   7.918  1.00 10.00           N
ATOM     10  CA  GLN A   3      16.590  18.147   7.740  1.00 10.00           C
ATOM     11  C   GLN A   3      15.154  17.904   8.207  1.00 10.00           C
ATOM     12  O   GLN A   3      14.813  18.175   9.356  1.00 10.00           O
ATOM     13  N   MET A   4      14.318  17.383   7.310  1.00 10.00           N
ATOM     14  CA  MET A   4      12.921  17.091   7.630  1.00 10.00           C
ATOM     15  C   MET A   4      12.750  15.585   7.849  1.00 10.00           C
ATOM     16  O   MET A   4      13.057  14.784   6.965  1.00 10.00           O
TER
ATOM     17  N   THR B   5      14.260  17.201  11.025  1.00 10.00           N
ATOM     18  CA  THR B   5      14.076  15.783  11.355  1.00 10.00           C
ATOM     19  C   THR B   5      12.612  15.405  11.550  1.00 10.00           C
ATOM     20  O   THR B   5      11.958  15.902  12.463  1.00 10.00           O
ATOM     21  N   GLN B   6      12.104  14.518  10.697  0.90 10.00           N
ATOM     22  CA  GLN B   6      10.712  14.084  10.797  0.90 10.00           C
ATOM     23  C   GLN B   6      10.571  12.704  11.424  0.90 10.00           C
ATOM     24  O   GLN B   6      11.336  11.787  11.118  0.90 10.00           O
ATOM     25  N   THR B   7       9.565  12.570  12.283  1.00 10.00           N
ATOM     26  CA  THR B   7       9.266  11.322  12.973  1.00 10.00           C
ATOM     27  C   THR B   7       7.756  11.127  12.915  1.00 10.00           C
ATOM     28  O   THR B   7       7.000  12.070  13.145  1.00 10.00           O
ATOM     29  N   PRO B   8       7.291   9.907  12.609  1.00 10.00           N
ATOM     30  CA  PRO B   8       8.060   8.698  12.310  1.00 10.00           C
ATOM     31  C   PRO B   8       8.406   8.619  10.827  1.00 10.00           C
ATOM     32  O   PRO B   8       8.058   9.514  10.054  1.00 10.00           O
ATOM     33  N   LEU B   9       9.093   7.553  10.427  1.10 10.00           N
ATOM     34  CA  LEU B   9       9.459   7.385   9.023  1.10 10.00           C
ATOM     35  C   LEU B   9       8.232   7.000   8.213  1.10 10.00           C
ATOM     36  O   LEU B   9       8.026   7.506   7.113  1.10 10.00           O
TER
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=good)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = "wk1995")
  f_obs = abs(xrs.structure_factors(d_min=1,
    algorithm="direct").f_calc()).set_observation_type_xray_amplitude()
  sfp = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfp.algorithm="direct"
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=bad)
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = "wk1995")
  #
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    xray_structure               = xrs,
    target_name                  = "ml",
    sf_and_grads_accuracy_params = sfp)
  alpha_beta = fmodel.alpha_beta()
  tg = mmtbx.utils.experimental_data_target_and_gradients(
    fmodel = fmodel,
    alpha_beta=alpha_beta)
  result = tg.group_occupancy_grads(
    pdb_hierarchy       = hierarchy,
    residues_per_window = 1)
  for r in result:
    print("chainID_resseqs: %s occupancy_grad: %-15.6f"%tuple(r))
  # get group gradients using custom selections
  selections = [flex.size_t([0,1,2,3]), flex.size_t([4,5,6,7,8])]
  result = tg.group_occupancy_grads(selections = selections)
  for i, r in enumerate(result):
    print("selection#: %s occupancy_grad: %-15.6f"%(str(i), r[1]))

def exercise_get_atom_selections(verbose=False):
  pdb_in = """\
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 212121
HETATM  115  O   HOH A  18       3.000   5.000   5.000  1.00 10.00           O
HETATM  115  O   HOH A  19       5.000   5.000   8.000  1.00 10.00           O
HETATM  115  O   HOH A  20       5.000   5.000   8.000  1.00 10.00           O
END"""
  log = null_out()
  if (verbose):
    log = sys.stdout
  model = mmtbx.model.manager(
      model_input = iotbx.pdb.input(lines=pdb_in.splitlines(), source_info=None))
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records=pdb_in.splitlines())
  selections1 = utils.get_atom_selections(
    model = model,
    selection_strings=["resseq 18", "resseq 19", "resseq 20"],
    parameter_name="refine.occupancy")
  try :
    selections2 = utils.get_atom_selections(
      model = model,
      selection_strings=["resseq 18:19", "resseq 19:20"],
      parameter_name="refine.occupancy")
  except Sorry as s :
    assert (str(s) == """\
One or more overlapping selections for refine.occupancy:
resseq 18:19
resseq 19:20""")
  else :
    raise Exception_expected

def exercise_f_000():
  pdb_str="""
CRYST1    5.000    5.000    5.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       0.000   0.000   0.000  1.00  5.00           C
END
"""
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  #
  f = utils.f_000(xray_structure=xrs, mean_solvent_density=0)
  assert approx_equal(f.f_000, 6, 1.e-3)
  #
  f = utils.f_000(mean_solvent_density=1, unit_cell_volume=125)
  assert approx_equal(f.f_000, 125)
  #
  f = utils.f_000(mean_solvent_density=0.25, unit_cell_volume=125,
                  solvent_fraction=0.3, xray_structure = xrs)
  assert approx_equal(f.f_000, 0.25*125*0.3+6, 1.e-3)
  assert approx_equal(f.solvent_fraction, 0.3)
  #
  f = utils.f_000(mean_solvent_density=0.25, unit_cell_volume=125,
                  xray_structure = xrs)
  assert approx_equal(f.f_000, 0.25*125*0.687355324074+6, 1.e-3)

def exercise_detect_link_problems():
  with open("tmp_mmtbx_utils_asn_nag.pdb", "w") as f:
    f.write("""\
CRYST1  124.702  124.702   71.573  90.00  90.00  90.00 P 4 21 2
ATOM   3196  N   ASN A 284      36.622 -19.654  35.782  1.00 19.63           N
ATOM   3197  CA  ASN A 284      36.491 -18.279  35.327  1.00 19.79           C
ATOM   3198  C   ASN A 284      35.037 -17.835  35.322  1.00 20.05           C
ATOM   3199  O   ASN A 284      34.751 -16.634  35.341  1.00 28.11           O
ATOM   3200  CB  ASN A 284      37.063 -18.130  33.915  1.00 21.07           C
ATOM   3201  CG  ASN A 284      38.549 -17.829  33.904  1.00 20.78           C
ATOM   3202  OD1 ASN A 284      39.040 -17.053  34.702  1.00 19.74           O
ATOM   3203  ND2 ASN A 284      39.263 -18.449  32.968  1.00 21.82           N
ATOM   3204  H   ASN A 284      36.875 -20.208  35.174  1.00 23.56           H
ATOM   3205 HD21 ASN A 284      38.893 -18.679  32.227  1.00 26.19           H
ATOM   3206  HA  ASN A 284      36.987 -17.719  35.944  1.00 19.79           H
ATOM   3207  HB2 ASN A 284      36.900 -18.947  33.418  1.00 21.07           H
ATOM   3208  HB3 ASN A 284      36.591 -17.419  33.454  1.00 21.07           H
ATOM   3209 HD22 ASN A 284      38.878 -18.991  32.422  1.00 21.82           H
HETATM 5988  C1  NAG A 467      40.601 -17.959  32.799  1.00 27.22           C
HETATM 5989  C2  NAG A 467      41.289 -19.314  32.714  1.00 22.16           C
HETATM 5990  C3  NAG A 467      42.783 -19.123  32.507  1.00 54.68           C
HETATM 5991  C4  NAG A 467      43.034 -18.265  31.278  1.00 23.55           C
HETATM 5992  C5  NAG A 467      42.261 -16.957  31.391  1.00 33.78           C
HETATM 5993  C6  NAG A 467      42.388 -16.125  30.141  1.00 24.49           C
HETATM 5994  C7  NAG A 467      41.114 -21.444  33.906  1.00 21.47           C
HETATM 5995  C8  NAG A 467      40.844 -22.107  35.214  1.00 20.34           C
HETATM 5996  N2  NAG A 467      41.041 -20.110  33.902  1.00 24.73           N
HETATM 5997  O3  NAG A 467      43.399 -20.391  32.338  1.00 54.77           O
HETATM 5998  O4  NAG A 467      44.417 -17.960  31.155  1.00 53.51           O
HETATM 5999  O5  NAG A 467      40.861 -17.215  31.600  1.00 31.39           O
HETATM 6000  O6  NAG A 467      41.470 -15.043  30.154  1.00 46.51           O
HETATM 6001  O7  NAG A 467      41.392 -22.081  32.897  1.00 22.76           O
HETATM 6002  H1  NAG A 467      40.952 -17.467  33.566  1.00 32.66           H
HETATM 6003  H2  NAG A 467      40.934 -19.790  31.940  1.00 26.59           H
HETATM 6004  H3  NAG A 467      43.163 -18.682  33.290  1.00 65.62           H
HETATM 6005  H4  NAG A 467      42.738 -18.746  30.482  1.00 28.26           H
HETATM 6006  H5  NAG A 467      42.608 -16.449  32.148  1.00 40.54           H
HETATM 6007  H61 NAG A 467      42.210 -16.687  29.363  1.00 29.39           H
HETATM 6008  H62 NAG A 467      43.296 -15.773  30.082  1.00 29.39           H
HETATM 6009  H81 NAG A 467      40.882 -23.076  35.101  1.00 24.41           H
HETATM 6010  H82 NAG A 467      39.958 -21.851  35.532  1.00 24.41           H
HETATM 6011  H83 NAG A 467      41.516 -21.829  35.865  1.00 24.41           H
HETATM 6012  HN2 NAG A 467      40.836 -19.681  34.680  1.00 29.67           H
HETATM 6013  HO3 NAG A 467      42.779 -20.998  32.145  1.00 65.72           H
HETATM 6014  HO4 NAG A 467      44.884 -18.471  31.711  1.00 64.21           H
HETATM 6015  HO6 NAG A 467      40.829 -15.206  30.746  1.00 55.81           H
""")
  result = utils.detect_hydrogen_nomenclature_problem(
    pdb_file="tmp_mmtbx_utils_asn_nag.pdb")
  assert (result.n_asn_hd22 == 0)

def exercise_corrupt_cryst1():
  """
  Consistency check for cs derived in different scenarious (inspired by PDB code
  2y9k).
  """
  pdb_str = """
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1          15
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
ATOM      1  N   GLY A  34     -74.292  12.386 -24.083  1.00  0.00           N
ATOM      2  CA  GLY A  34     -74.465  11.623 -25.332  1.00  0.00           C
ATOM      3  C   GLY A  34     -73.114  11.306 -25.993  1.00  0.00           C
ATOM      4  O   GLY A  34     -72.196  12.128 -25.997  1.00  0.00           O
"""
  of = open("tmp_exercise_corrupt_cryst1.pdb", "w")
  print(pdb_str, file=of)
  of.close()
  base_arg = ["tmp_exercise_corrupt_cryst1.pdb"]
  for extra_arg in [[], ["bla=bla"], ["bla"]]:
    o = utils.process_command_line_args(args=base_arg+extra_arg)
    assert o.crystal_symmetry is None, o.crystal_symmetry
    o = utils.process_command_line_args(args=extra_arg+base_arg)
    assert o.crystal_symmetry is None

def exercise_rotatable_bonds():
  pdb_str = """
CRYST1  124.792  235.089   82.032  90.00  90.00  90.00 C 2 2 21
ATOM      1  N   LYS C   1     -31.148 -19.390 -11.207  1.00  0.00           N
ATOM      2  CA  LYS C   1     -31.030 -18.256 -10.299  1.00  0.00           C
ATOM      3  C   LYS C   1     -31.972 -18.403  -9.110  1.00  0.00           C
ATOM      4  O   LYS C   1     -31.567 -18.231  -7.960  1.00  0.00           O
ATOM      5  CB  LYS C   1     -31.312 -16.946 -11.035  1.00  0.00           C
ATOM      6  CG  LYS C   1     -31.183 -15.698 -10.172  1.00  0.00           C
ATOM      7  CD  LYS C   1     -31.387 -14.435 -10.995  1.00  0.00           C
ATOM      8  CE  LYS C   1     -31.277 -13.186 -10.131  1.00  0.00           C
ATOM      9  NZ  LYS C   1     -31.474 -11.943 -10.922  1.00  0.00           N
ATOM     10  HA  LYS C   1     -30.010 -18.227  -9.914  1.00  0.00           H
ATOM     11  HB1 LYS C   1     -30.622 -16.845 -11.874  1.00  0.00           H
ATOM     12  HB2 LYS C   1     -32.323 -16.968 -11.444  1.00  0.00           H
ATOM     13  HG1 LYS C   1     -31.929 -15.728  -9.377  1.00  0.00           H
ATOM     14  HG2 LYS C   1     -30.193 -15.670  -9.719  1.00  0.00           H
ATOM     15  HD1 LYS C   1     -30.634 -14.388 -11.783  1.00  0.00           H
ATOM     16  HD2 LYS C   1     -32.373 -14.458 -11.458  1.00  0.00           H
ATOM     17  HE1 LYS C   1     -32.025 -13.224  -9.341  1.00  0.00           H
ATOM     18  HE2 LYS C   1     -30.291 -13.153  -9.665  1.00  0.00           H
ATOM     19  HZ1 LYS C   1     -31.393 -11.140 -10.314  1.00  0.00           H
ATOM     20  HZ2 LYS C   1     -30.771 -11.889 -11.646  1.00  0.00           H
ATOM     21  HZ3 LYS C   1     -32.391 -11.953 -11.344  1.00  0.00           H
ATOM     22  H 1 LYS C   1     -30.521 -19.268 -11.976  1.00  0.00           H
ATOM     23  H 2 LYS C   1     -30.920 -20.233 -10.720  1.00  0.00           H
ATOM     24  H 3 LYS C   1     -32.087 -19.447 -11.549  1.00  0.00           H
TER
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input = pdb_inp, log=null_out())
  model.process(make_restraints=True)
  residue = model.get_hierarchy().only_residue()
  r = mmtbx.utils.rotatable_bonds.axes_and_atoms_aa_specific(
      residue = residue, mon_lib_srv = model.get_mon_lib_srv(), log=null_out())
  assert r is None

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_corrupt_cryst1()
  exercise_d_data_target_d_atomic_params()
  exercise_d_data_target_d_atomic_params2()
  exercise_get_atom_selections(verbose=verbose)
  exercise_f_000()
  exercise_detect_link_problems()
  exercise_rotatable_bonds()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
