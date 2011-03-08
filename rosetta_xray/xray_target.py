
from __future__ import division
import mmtbx.maps.utils
import mmtbx.f_model
import mmtbx.utils
import iotbx.phil
from cctbx import xray
from scitbx.array_family import flex
from libtbx.str_utils import make_header
from libtbx.utils import Usage
import time
import os
import sys

master_phil = iotbx.phil.parse("""
include scope mmtbx.utils.cmdline_input_phil_str
output {
  write_grads = True
    .type = bool
  write_eff = True
    .type = bool
  prefix = gradients
    .type = str
  serial = 1
    .type = int
  output_dir = None
    .type = path
  ccp4_map {
    write_map = False
      .type = bool
    map_type = 2mFo-DFc
      .type = str
    map_file_name = 2mFo-DFc.ccp4
      .type = path
    fill_missing_f_obs = False
      .type = bool
    grid_resolution_factor = 1/3
      .type = float
  }
}
options {
  bulk_solvent_and_scale = True
    .type = bool
  force_bss = False
    .type = bool
  remove_outliers = True
    .type = bool
  target_name = *ml lsq
    .type = choice(multi=False)
  compute_gradients = True
    .type = bool
}
include scope mmtbx.command_line.fmodel.fmodel_from_xray_structure_params
bss {
  include scope mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params
}
""", process_includes=True)

class xray_target (object) :
  def __init__ (self,
                params,
                xray_structure,
                f_obs,
                r_free_flags,
                out=sys.stdout,
                ignore_bulk_solvent=False) :
    self.params = params
    self.out = out
    t1 = time.time()
    params.bss.apply_back_trace_of_b_cart=False
    if ((not params.options.bulk_solvent_and_scale) and
        (not ignore_bulk_solvent)) :
      assert (not None in [params.fmodel.k_sol,params.fmodel.b_sol,
              params.fmodel.b_cart])
    target_name = params.options.target_name
    if (target_name == "lsq") :
      target_name = "ls_wunit_k1"
    fmodel = mmtbx.f_model.manager(
      xray_structure               = xray_structure,
      sf_and_grads_accuracy_params = params.structure_factors_accuracy,
      r_free_flags                 = r_free_flags,
      mask_params                  = params.mask,
      f_obs                        = f_obs,
      target_name                  = target_name,
      k_sol                        = params.fmodel.k_sol,
      b_sol                        = params.fmodel.b_sol,
      b_cart                       = params.fmodel.b_cart)
    #fmodel_info = fmodel.info()
    #fmodel_info.show_rfactors_targets_scales_overall(out=out)
    if params.options.remove_outliers :
      print >> out, ""
      fmodel = fmodel.remove_outliers(show=True, log=out)
    if params.options.bulk_solvent_and_scale :
      make_header("Bulk solvent and scaling", out=out)
      fmodel.update_solvent_and_scale(params=params.bss)
      fmodel_info = fmodel.info()
      fmodel_info.show_rfactors_targets_scales_overall(out=out)
    else :
      make_header("Update X-ray structure", out=out)
      fmodel.apply_back_b_iso()
      fmodel.update_xray_structure(update_f_mask=True)
    self.fmodel = fmodel
    self.xray_structure = fmodel.xray_structure
    self.flag_apply_shifts = False

  def prepare_for_minimization (self) :
    make_header("Preparing for minimization", out=self.out)
    xrs = self.xray_structure
    xrs.scatterers().flags_set_grads(state=False)
    selection = flex.bool(xrs.sites_cart().size(), True)
    #print selection.iselection().size()
    xrs.scatterers().flags_set_grad_site(iselection=selection.iselection())
    self.x = flex.double(xrs.n_parameters(), 0)
    self._scatterers_start = xrs.scatterers()
    self.fmodel.update_xray_structure(xray_structure=self.xray_structure)
    self.target_functor = self.fmodel.target_functor()
    self.target_functor.prepare_for_minimization()
    self.fmodel.info().show_targets(out=self.out, text="Target values")

  def compute_target (self, compute_gradients=True, use_vec3_array=False) :
    if (getattr(self, "target_functor", None) is None) :
      raise RuntimeError("Please call prepare_for_minimization() before "+
        "computing target and gradients.")
    if (self.flag_apply_shifts) :
      self.apply_shifts()
    out = self.out
    fmodel = self.fmodel
    fmodel.update_xray_structure(xray_structure=self.xray_structure,
      update_f_calc=True)
    t1 = time.time()
    tfx_r = self.target_functor(compute_gradients=compute_gradients)
    self.f = tfx_r.target_work()
    if compute_gradients :
      grads = tfx_r.gradients_wrt_atomic_parameters().packed()
      assert (grads.size() == self.x.size())
      if (use_vec3_array) :
        self.g = flex.vec3_double(grads)
      else :
        self.g = grads
    else :
      if (use_vec3_array) :
        self.g = flex.vec3_double(int(self.x.size() / 3), (0.0,0.0,0.0))
      else :
        self.g = flex.double(self.x.size(), 0.0)
    t2 = time.time()
    return self.f, self.g

  def compute_functional_and_gradients (self) :
    self.apply_shifts()
    return self.compute_target()

  def compute_functional_and_gradients_rosetta (self) :
    self.apply_shifts()
    return self.compute_target(use_vec3_array=True)

  def target (self) :
    return self.f

  def gradients (self) :
    return self.g

  def apply_shifts (self) :
    xrs = self.xray_structure
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = xrs.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    site_symmetry_table = xrs.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      scatterers_shifted[i_seq].site = crystal.correct_special_position(
        crystal_symmetry = self.xray_structure,
        special_op       = site_symmetry_table.get(i_seq).special_op(),
        site_frac        = scatterers_shifted[i_seq].site,
        site_label       = scatterers_shifted[i_seq].label,
        tolerance        = self.correct_special_position_tolerance)
    xrs.replace_scatterers(scatterers = scatterers_shifted)
    return None

  def update_sites (self, sites_cart) :
    self.xray_structure.set_sites_cart(sites_cart)
    self.fmodel.update_xray_structure(
      xray_structure=self.xray_structure,
      update_f_calc=False,
      update_f_mask=True)

  def update_sites_1d (self, sites_cart_as_1d) :
    if isinstance(sites_cart_as_1d, list) :
      assert (len(sites_cart_as_1d) == self.x.size())
      sites_cart_as_1d = flex.double(sites_cart_as_1d)
    else :
      assert (sites_cart_as_1d.size() == self.x.size())
    sites_cart = flex.vec3_double(sites_cart_as_1d)
    self.update_sites(sites_cart)

  def set_shifts (self, site_shifts) :
    assert (site_shifts.size() == self.x.size())
    self.x = site_shifts

  def clean_up_after_minimization (self) :
    self.apply_shifts()
    info = self.fmodel.info()
    info.show_targets(out=self.out, text="Target values after minimization")

  def update_pdb_hierarchy (self, pdb_hierarchy) :
    pdb_hierarchy.atoms().set_xyz(self.xray_structure.sites_cart())

  def compute_ccp4_map (self,
                        map_type="2mFo-DFc",
                        file_name="2mFo-DFc.ccp4",
                        fill_missing_f_obs=False,
                        grid_resolution_factor=1/3) :
    map_manager = self.fmodel.electron_density_map(
      fill_missing_f_obs=fill_missing_f_obs,
      fill_mode="dfmodel")
    map_coeffs = map_manager.map_coefficients(map_type=map_type)
    if map_coeffs.anomalous_flag() :
      map_coeffs = map_coeffs.average_bijvoet_mates()
    fft_map = map_coeffs.fft_map(resolution_factor=grid_resolution_factor)
    fft_map.apply_sigma_scaling()
    mmtbx.maps.utils.write_ccp4_map(
      sites_cart=self.xray_structure.sites_cart(),
      unit_cell=fft_map.unit_cell(),
      map_data=fft_map.real_map(),
      n_real=fft_map.n_real(),
      file_name=file_name)

  def write_files (self) :
    out = self.out
    params = self.params
    fmodel = self.fmodel
    print >> out, ""
    params.fmodel.b_cart = fmodel.b_cart()
    params.fmodel.k_sol = fmodel.k_sol()
    params.fmodel.b_sol = fmodel.b_sol()
    if (params.output.output_dir is None) :
      params.output.output_dir = os.getcwd()
    file_base = os.path.join(params.output.output_dir,
      "%s_%d" % (params.output.prefix, params.output.serial))
    params.output.serial += 1
    if (not params.options.force_bss) :
      params.options.bulk_solvent_and_scale = False
    if params.output.write_eff :
      f = open("%s.eff" % file_base, "w")
      final_phil = master_phil.format(python_object=params)
      final_phil.show(out=f)
      f.close()
      print >> out, "Parameters for next run:"
      print >> out, "  %s.eff" % file_base
    if params.output.write_grads :
      grads_out = open("%s.xyz" % file_base, "w")
      print >> grads_out, "# TARGET = %.12g" % self.f
      if isinstance(self.g, flex.double) :
        grads = flex.vec3_double(self.g)
      else :
        grads = self.g
        for site_grads in grads :
          print >> grads_out, "%e %e %e" % site_grads
      grads_out.close()
      print >> out, "Gradients:"
      print >> out, "  %s.xyz" % file_base
    if params.output.ccp4_map.write_map :
      self.compute_ccp4_map(
        map_type=params.output.ccp4_map.map_type,
        file_name=params.output.ccp4_map.map_file_name,
        fill_missing_f_obs=params.output.ccp4_map.fill_missing_f_obs,
        grid_resolution_factor=params.output.ccp4_map.grid_resolution_factor)
      print >> out, "CCP4 map:"
      print >> out, "  %s" % params.output.ccp4_map.map_file_name

def run (args, out=sys.stdout) :
  if (len(args) == 0) :
    raise Usage("""
phenix.compute_xray_gradients [model.pdb] [data.mtz] [params.eff] [options...]

This program will write out a parameter (.eff) file when it is finished,
which can be used for the next round without bulk-solvent correction.
""")
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    args=args,
    master_phil=master_phil,
    out=out,
    process_pdb_file=False,
    create_fmodel=False)
  params = cmdline.params
  f_obs = cmdline.f_obs
  r_free_flags = cmdline.r_free_flags
  xray_structure = cmdline.xray_structure
  target_evaluator = xray_target(
    params=cmdline.params,
    xray_structure=cmdline.xray_structure,
    f_obs=cmdline.f_obs,
    r_free_flags=cmdline.r_free_flags,
    out=out)
  target_evaluator.prepare_for_minimization()
  target_evaluator.compute_target(
    compute_gradients=params.options.compute_gradients,
    use_vec3_array=True)
  target_evaluator.write_files()
  return target_evaluator

########################################################################
# REGRESSION TESTING
def exercise () :
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
  if ("--test" in sys.argv[1:]) :
    exercise()
  else :
    run(sys.argv[1:])
