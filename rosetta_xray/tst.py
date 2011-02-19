
import os
import sys

def exercise () :
  from rosetta_xray.xray_target import xray_target, master_phil
  import iotbx.pdb
  import scitbx.lbfgs
  import cStringIO
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
CRYST1   16.000   16.000   16.000  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ALA A  10      -0.961  12.543   2.657  1.00 19.69           N
ATOM      2  CA  ALA A  10       0.210  13.408   2.561  1.00 19.69           C
ATOM      3  C   ALA A  10      -0.191  14.849   2.274  1.00 19.69           C
ATOM      4  O   ALA A  10      -0.938  15.120   1.332  1.00 19.69           O
ATOM      5  CB  ALA A  10       1.180  12.913   1.505  1.00 19.69           C
ATOM      6  N   ALA A  11       0.307  15.771   3.091  1.00 19.69           N
ATOM      7  CA  ALA A  11      -0.020  17.184   2.942  1.00 16.70           C
ATOM      8  C   ALA A  11       0.030  17.609   1.480  1.00 16.70           C
ATOM      9  O   ALA A  11      -0.840  18.340   1.007  1.00 16.70           O
ATOM     10  CB  ALA A  11       0.936  18.074   3.758  1.00 16.70           C
ATOM     11  N   ALA A  12       1.054  17.149   0.770  1.00 16.70           N
ATOM     12  CA  ALA A  12       1.250  17.510  -0.624  1.00 16.70           C
ATOM     13  C   ALA A  12       0.062  17.094  -1.475  1.00 19.69           C
ATOM     14  O   ALA A  12      -0.254  17.737  -2.453  1.00 13.88           O
ATOM     15  CB  ALA A  12       2.532  16.876  -1.167  1.00 35.67           C
END""")
  xrs = pdb_in.xray_structure_simple()
  f_calc = xrs.structure_factors(d_min=1.5).f_calc().as_amplitude_array()
  sites_orig = xrs.sites_cart().deep_copy()
  xrs.shake_sites_in_place(rms_difference=0.5)
  sites_shaken = xrs.sites_cart().deep_copy()
  rmsd1 = xrs.sites_cart().rms_difference(sites_orig)
  params = master_phil.fetch().extract()
  params.options.target_name = "lsq"
  params.options.bulk_solvent_and_scale = False
  params.output.ccp4_map.write_map = True
  tf = xray_target(
    params=params,
    xray_structure=xrs,
    f_obs=f_calc,
    r_free_flags=f_calc.generate_r_free_flags(),
    out=cStringIO.StringIO(),
    ignore_bulk_solvent=True)
  r_free1 = tf.fmodel.r_free()
  tf.prepare_for_minimization()
  # test replacing sites
  tf.compute_target(compute_gradients=False)
  target1 = tf.target()
  tf.update_sites_1d(sites_orig.as_double())
  tf.compute_target(compute_gradients=False)
  target2 = tf.target()
  tf.update_sites_1d(list(sites_shaken.as_double()))
  tf.compute_target(compute_gradients=False)
  target3 = tf.target()
  assert (target3 == target1) and (target2 < 0.0001)
  # test actual LBFGS minimization
  tf.flag_apply_shifts = True
  minimizer = scitbx.lbfgs.run(tf)
  tf.clean_up_after_minimization()
  rmsd2 = tf.xray_structure.sites_cart().rms_difference(sites_orig)
  r_free2 = tf.fmodel.r_free()
  print "starting model: rmsd=%.4f r_free=%.4f" % (rmsd1, r_free1)
  print "final model:    rmsd=%.4f r_free=%.4f" % (rmsd2, r_free2)
  assert (r_free2 < 0.01) and (rmsd2 < 0.01)
  hierarchy = pdb_in.construct_hierarchy()
  tf.update_pdb_hierarchy(hierarchy)
  tf.write_files()
  assert os.path.isfile("2mFo-DFc.ccp4")
  print "OK"

if __name__ == "__main__" :
  exercise()
