
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
  # GUI multi-criterion plot path: no real-space data here, so every residue
  # (and every bin gap) is padded with numpy.nan
  import math
  y_limits = mc.get_y_limits()
  assert set(y_limits.keys()) == {"rho", "b", "cc"}
  assert all(math.isnan(v) for v in y_limits["rho"]), y_limits
  binner = mc.binned_data()
  assert len(binner.bins) > 0
  for res_bin in binner.bins:
    rs_values = res_bin.get_real_space_plot_values()
    assert rs_values.shape == (4, res_bin.n_res())
    assert all(math.isnan(v) for v in rs_values.flat)
    outlier_values = res_bin.get_outlier_plot_values()
    assert outlier_values.shape == (4, res_bin.n_res())
  # gaps between residues are padded with None and also plot as numpy.nan
  from mmtbx.validation import graphics
  res_bin = graphics.residue_bin()
  res_bin.add_residue(mc.data()[0])
  res_bin.add_empty(2)
  assert res_bin.n_res() == 3
  for values in (res_bin.get_real_space_plot_values(),
                 res_bin.get_outlier_plot_values()):
    assert values.shape == (4, 3)
    assert all(math.isnan(v) for v in values[:, 1:].flat)
  assert (result.neutron_stats is None)
  mpscore = result.molprobity_score()
  # percentiles
  out4 = StringIO()
  result.show_summary(out=out4, show_percentiles=True)
  # Clashscore (and its percentile) depends on the H-placement engine: reduce2
  # places H differently than reduce1, so gate the expected value on the switch
  # for now (until reduce2 becomes the default).
  from mmtbx import hydrogens as reduce_switch
  expected_clashscore = (
    "  Clashscore            =  49.96 (percentile: 0.2)"
    if reduce_switch.use_old_reduce() else
    "  Clashscore            =  49.17 (percentile: 0.3)")
  assert expected_clashscore in out4.getvalue(), out4.getvalue()
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

def exercise_nucleic_acid_multi_criterion():
  # a model without protein, with data: the multi-criterion plot must get
  # real-space values for the nucleic-acid residues too (real_space.results
  # only lists protein residues); with all values NaN the GUI plot fails with
  # "Axis limits cannot be NaN or Inf"
  import math
  import iotbx.pdb
  import mmtbx.f_model
  from mmtbx.regression import model_1zew_dna
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(source_info=None, lines=model_1zew_dna))
  model.process(make_restraints=True)
  xrs = model.get_xray_structure()
  f_obs = xrs.structure_factors(d_min=2.0).f_calc().as_amplitude_array()
  fmodel = mmtbx.f_model.manager(
    f_obs=f_obs,
    r_free_flags=f_obs.generate_r_free_flags(fraction=0.1),
    xray_structure=xrs)
  result = mmtbx.validation.molprobity.molprobity(
    model=model, fmodel=fmodel, outliers_only=False)
  assert (result.real_space is not None)
  assert (result.nqh_flips is None)
  mc = result.as_multi_criterion_view()
  residues = mc.data()
  assert (len(residues) == 3)
  for residue in residues:
    b_iso, cc, two_fofc, fmodel_value = residue.get_real_space_plot_values()
    assert (not math.isnan(cc)), residue.id_str()
    assert (not math.isnan(b_iso)), residue.id_str()
  y_limits = mc.get_y_limits()
  for key in ("rho", "b", "cc"):
    assert (not any(math.isnan(v) for v in y_limits[key])), y_limits
  # views matching the residue-type menu of the GUI (Protein, Other, Water,
  # Everything): same residues as the corresponding real-space table
  assert (result.waters is not None) and (len(result.waters.results) == 1)
  for residue_type, expected_ids in [
      ("other", [" A   1 ", " A   2 ", " A   3 "]),
      ("water", [" A  43 "]),
      ("everything", [" A   1 ", " A   2 ", " A   3 ", " A  43 "])]:
    mc = result.as_multi_criterion_view(residue_type=residue_type)
    residues = mc.data()
    assert ([r.residue_group_id_str() for r in residues] == expected_ids), \
      (residue_type, [r.residue_group_id_str() for r in residues])
    for residue in residues:
      values = residue.get_real_space_plot_values()
      assert (not any(math.isnan(v) for v in values)), (residue_type,
        residue.id_str())
    y_limits = mc.get_y_limits()
    for key in ("rho", "b", "cc"):
      assert (not any(math.isnan(v) for v in y_limits[key])), y_limits
      # usable axis limits even for a single residue (the water view)
      assert (y_limits[key][0] < y_limits[key][1]), (residue_type, y_limits)
    assert (len(mc.binned_data().bins) == 1)
  # no residues of the selected type: an empty view that can still be
  # plotted (no exception, no bins, no limits)
  mc = result.as_multi_criterion_view(residue_type="protein")
  assert (mc.data() == [])
  assert (mc.binned_data().bins == [])
  assert (mc.get_y_limits() == {"rho": (None, None), "b": (None, None),
                                "cc": (None, None)})
  # the default view is unchanged (protein and nucleic acid chains, no water)
  mc = result.as_multi_criterion_view()
  assert (len(mc.data()) == 3)

def exercise_multi_criterion_plot():
  # draw the plot without a GUI (Agg backend), including the empty case
  try:
    import matplotlib
  except ImportError:
    print("matplotlib not available, skipping exercise_multi_criterion_plot")
    return
  matplotlib.use("Agg")
  from matplotlib.figure import Figure
  from matplotlib.backends.backend_agg import FigureCanvasAgg
  from matplotlib.font_manager import FontProperties
  import matplotlib.ticker
  from mmtbx.validation import graphics
  import iotbx.pdb
  import mmtbx.f_model
  from mmtbx.regression import model_1zew_dna

  class headless_plot(graphics.multi_criterion_plot_mixin):
    def __init__(self, binner, y_limits):
      graphics.multi_criterion_plot_mixin.__init__(self, binner, y_limits)
      self.figure = Figure(figsize=(8, 5))
      self.canvas = FigureCanvasAgg(self.figure)
      self.null_fmt = matplotlib.ticker.NullFormatter()
    def get_font(self, font_type):
      return FontProperties(size=8)

  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(source_info=None, lines=model_1zew_dna))
  model.process(make_restraints=True)
  xrs = model.get_xray_structure()
  f_obs = xrs.structure_factors(d_min=2.0).f_calc().as_amplitude_array()
  fmodel = mmtbx.f_model.manager(
    f_obs=f_obs,
    r_free_flags=f_obs.generate_r_free_flags(fraction=0.1),
    xray_structure=xrs)
  result = mmtbx.validation.molprobity.molprobity(
    model=model, fmodel=fmodel, outliers_only=False)
  for residue_type in ["other", "water", "everything"]:
    mc = result.as_multi_criterion_view(residue_type=residue_type)
    plot = headless_plot(mc.binned_data(), mc.get_y_limits())
    plot.plot_range(0)
    assert (plot._current_bin is not None)
    assert (len(plot.figure.axes) == 3), residue_type # density, CC, B
  # empty selection: an empty plot instead of an exception
  mc = result.as_multi_criterion_view(residue_type="protein")
  plot = headless_plot(mc.binned_data(), mc.get_y_limits())
  plot.plot_range(0)
  assert (plot._current_bin is None)
  assert (len(plot.figure.axes) == 1)
  assert (len(plot.figure.axes[0].lines) == 0)

if (__name__ == "__main__"):
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping tests: probe not configured")
  else :
    exercise_protein()
    exercise_rna()
    exercise_nucleic_acid_multi_criterion()
    exercise_multi_criterion_plot()
    print("OK")
