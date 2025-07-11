
from __future__ import absolute_import, division, print_function
from mmtbx.command_line import molprobity
import mmtbx.model
import mmtbx.validation.molprobity
from libtbx.easy_pickle import loads, dumps, dump
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import libtbx.load_env
from six.moves import cStringIO as StringIO
import os.path as op

# test on protein - we need real model/data for this
def exercise_protein():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3ifk.pdb",
    test=op.isfile)
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/3ifk.mtz",
    test=op.isfile)
  if (pdb_file is None):
    print("phenix_regression not available, skipping.")
    return
  args1 = [
    pdb_file,
    "outliers_only=True",
    "output.prefix=tst_molprobity",
    "--pickle",
    "flags.xtriage=True",
  ]
  result = molprobity.run(args=args1, out=null_out()).validation
  out1 = StringIO()
  result.show(out=out1)
  result = loads(dumps(result))
  out2 = StringIO()
  result.show(out=out2)
  assert (result.nqh_flips.n_outliers == 6)
  assert (not "RNA validation" in out2.getvalue())
  assert (out2.getvalue() == out1.getvalue())
  dump("tst_molprobity.pkl", result)
  mc = result.as_multi_criterion_view()
  assert (result.neutron_stats is None)
  mpscore = result.molprobity_score()
  # percentiles
  out4 = StringIO()
  result.show_summary(out=out4, show_percentiles=True)
  assert ("""  Clashscore            =  49.96 (percentile: 0.2)""" in
          out4.getvalue())
  # misc
  assert approx_equal(result.r_work(), 0.237) # from PDB header
  assert approx_equal(result.r_free(), 0.293) # from PDB header
  assert approx_equal(result.d_min(), 2.03)   # from PDB header
  assert (result.d_max_min() is None)
  assert approx_equal(result.rms_bonds(), 0.02586, 1e-5)
  assert approx_equal(result.rms_angles(), 2.35285, 1e-5)
  assert approx_equal(result.rama_favored(), 96.47059)
  assert (result.cbeta_outliers() == 10)
  assert approx_equal(result.molprobity_score(), 3.39, eps=0.01)
  summary = result.summarize()
  gui_fields = list(summary.iter_molprobity_gui_fields())
  assert (len(gui_fields) == 6)
  #result.show()
  assert (str(mc.data()[2]) == ' A   5  THR  rota,cb,clash')
  import mmtbx.validation.molprobity
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(pdb_file)
  model = mmtbx.model.manager(pdb_in)
  result = mmtbx.validation.molprobity.molprobity(model)
  out3 = StringIO()
  result.show_summary(out=out3)
  assert  """\
  Ramachandran outliers =   1.76 %
                favored =  96.47 %
  Rotamer outliers      =  20.00 %
""" in out3.getvalue()
  # now with data
  args2 = args1 + [ hkl_file, "--maps" ]
  result, cmdline = molprobity.run(args=args2,
    out=null_out(),
    return_input_objects=True)
  out = StringIO()
  result.show(out=out)
  stats = result.get_statistics_for_phenix_gui()
  #print stats
  stats = result.get_polygon_statistics(["r_work","r_free","adp_mean_all",
    "angle_rmsd", "bond_rmsd", "clashscore"])
  #print stats
  assert approx_equal(result.r_work(), 0.2291, eps=0.001)
  assert approx_equal(result.r_free(), 0.2804, eps=0.001)
  assert approx_equal(result.d_min(), 2.0302, eps=0.0001)
  assert approx_equal(result.d_max_min(), [34.546125, 2.0302], eps=0.0001)
  assert approx_equal(result.rms_bonds(), 0.02586, 1e-5)
  assert approx_equal(result.rms_angles(), 2.35285, 1e-5)
  assert approx_equal(result.rama_favored(), 96.47059)
  assert (result.cbeta_outliers() == 10)
  assert approx_equal(result.unit_cell().parameters(),
          (55.285, 58.851, 67.115,90,90,90))
  assert (str(result.space_group_info()) == "P 21 21 21")
  bins = result.fmodel_statistics_by_resolution()
  assert (len(bins) == 10)
  assert approx_equal(result.atoms_to_observations_ratio(), 0.09755,
    eps=0.0001)
  assert approx_equal(result.b_iso_mean(), 31.11739)
  assert op.isfile("tst_molprobity_maps.mtz")
  bins = result.fmodel_statistics_by_resolution()
  #bins.show()
  bin_plot = result.fmodel_statistics_graph_data()
  lg = bin_plot.format_loggraph()
  # fake fmodel_neutron
  fmodel_neutron = cmdline.fmodel.deep_copy()
  result2 = mmtbx.validation.molprobity.molprobity(
    cmdline.model,
    fmodel=cmdline.fmodel,
    fmodel_neutron=fmodel_neutron,
    nuclear=True,
    keep_hydrogens=True)
  stats = result2.get_statistics_for_phenix_gui()
  assert ('R-work (neutron)' in [ label for (label, stat) in stats ])

def exercise_rna():
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb2goz_refmac_tls.ent",
    test=op.isfile)
  if (regression_pdb is None):
    print("Skipping exercise_regression(): input pdb (pdb2goz_refmac_tls.ent) not available")
    return
  result = molprobity.run(args=[regression_pdb], out=null_out()).validation
  assert (result.rna is not None)
  out = StringIO()
  result.show(out=out)
  assert ("2/58 pucker outliers present" in out.getvalue())
  result = loads(dumps(result))
  out2 = StringIO()
  result.show(out=out2)
  assert (out2.getvalue() == out.getvalue())

if (__name__ == "__main__"):
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping tests: probe not configured")
  else :
    exercise_protein()
    exercise_rna()
    print("OK")
