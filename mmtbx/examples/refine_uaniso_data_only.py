from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from libtbx import group_args
from cctbx import xray
import scitbx.lbfgs
from libtbx import adopt_init_args
import cctbx
from six.moves import range

pdb_str = """
CRYST1   12.000   11.000   13.000  80.00  70.00 100.00 P 1
SCALE1      0.083333  0.014694 -0.035164        0.00000
SCALE2      0.000000  0.092312 -0.024020        0.00000
SCALE3      0.000000  0.000000  0.084586        0.00000
ATOM      1  CB  PHE A   1       7.353   5.743   7.446  1.00 11.07           C
ANISOU    1  CB  PHE A   1     1417   1711   1077   -802   -534    562       C
ATOM      2  CG  PHE A   1       6.587   5.028   8.521  1.00 12.41           C
ATOM      3  CD1 PHE A   1       5.463   4.281   8.210  1.00 15.10           C
ANISOU    3  CD1 PHE A   1     2242   1692   1805   -865   -520    173       C
ATOM      4  CD2 PHE A   1       6.993   5.104   9.843  1.00 11.77           C
ATOM      5  CE1 PHE A   1       4.758   3.623   9.198  1.00 15.96           C
ATOM      6  CE2 PHE A   1       6.292   4.449  10.836  1.00 12.44           C
ANISOU    6  CE2 PHE A   1     1794   1178   1756   -466   -772     83       C
ATOM      7  CZ  PHE A   1       5.173   3.707  10.513  1.00 14.49           C
ANISOU    7  CZ  PHE A   1     2230   1388   1889   -462   -737     32       C
ATOM      8  C   PHE A   1       7.886   7.946   6.389  1.00 15.51           C
ANISOU    8  C   PHE A   1     1740   2635   1517   -904   -600    967       C
ATOM      9  O   PHE A   1       8.151   7.695   5.214  1.00 16.93           O
ANISOU    9  O   PHE A   1     1943   2817   1671  -1003   -687   1048       O
ATOM     10  OXT PHE A   1       8.501   8.858   6.941  1.00 19.45           O
ATOM     11  N   PHE A   1       5.580   7.078   6.395  1.00 13.11           N
ANISOU   11  N   PHE A   1     1400   1945   1635   -826   -589    838       N
ATOM     12  CA  PHE A   1       6.829   7.148   7.143  1.00 12.44           C
TER
END
"""

class minimizer(object):
  def __init__(self,
        fmodel,
        max_iterations=25):
    self.fmodel = fmodel
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    self.fmodel.xray_structure.scatterers().flags_set_grad_u_aniso(
      iselection = self.fmodel.xray_structure.use_u_aniso().iselection())
    self.fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
      iselection = self.fmodel.xray_structure.use_u_iso().iselection())
    self.x = flex.double(self.fmodel.xray_structure.n_parameters(), 0)
    self._scatterers_start = self.fmodel.xray_structure.scatterers()
    self.call_count = 0
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    self.fmodel.xray_structure.tidy_us()
    self.apply_shifts()
    del self._scatterers_start
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

  def apply_shifts(self):
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = self.fmodel.xray_structure.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    self.fmodel.xray_structure.replace_scatterers(
      scatterers = scatterers_shifted)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

  def compute_functional_and_gradients(self):
    self.apply_shifts()
    tgx = func(fmodel=self.fmodel)
    f = tgx.target_work()
    g = tgx.grads_u_anisos()
    xray.minimization.add_gradients(
      scatterers        = self.fmodel.xray_structure.scatterers(),
      xray_gradients    = g)
    self.call_count += 1
    return f, g

class func(object):
  def __init__(self, fmodel):
    adopt_init_args(self, locals())
    weights = flex.double(self.fmodel.f_obs().data().size(), 1.0)
    self.core = xray.target_functors.least_squares(
      compute_scale_using_all_data = False,
      f_obs                        = self.fmodel.f_obs(),
      r_free_flags                 = self.fmodel.r_free_flags(),
      weights                      = weights,
      scale_factor                 = 1)
    self.r = self.core(f_calc = self.fmodel.f_model(), compute_gradients=True)
    self.d_target_d_f_calc = self.r.gradients_work() # XXX needs scales
    self.ge = cctbx.xray.structure_factors.gradients(
      miller_set = self.fmodel.f_obs())

  def target_work(self):
    return self.r.target_work()

  def grads_u_anisos(self):
    return self.ge(
      u_iso_refinable_params = None,
      d_target_d_f_calc      = self.d_target_d_f_calc,
      xray_structure         = self.fmodel.xray_structure,
      n_parameters           = self.fmodel.xray_structure.n_parameters(),
      miller_set             = self.fmodel.f_obs_work(),
      algorithm              = "direct").packed()

def get_inputs(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  ph = pdb_inp.construct_hierarchy()
  xrs = ph.extract_xray_structure(crystal_symmetry =
    pdb_inp.crystal_symmetry())
  return group_args(pdb_hierarchy = ph, xray_structure = xrs)

def run():
  # get xray_structure from PDB file
  inp = get_inputs(pdb_str = pdb_str)
  if(1):
    inp.pdb_hierarchy.adopt_xray_structure(inp.xray_structure)
    inp.pdb_hierarchy.write_pdb_file(file_name="start.pdb")
  # simulate poor starting model
  xrs_poor = inp.xray_structure.deep_copy_scatterers()
  xrs_poor.shake_adp(aniso_spread=1.5, random_u_cart_scale=10.0)
  if(1):
    inp.pdb_hierarchy.adopt_xray_structure(xrs_poor)
    inp.pdb_hierarchy.write_pdb_file(file_name="poor.pdb")
  # simulate Fobs
  f_obs = abs(inp.xray_structure.structure_factors(
    d_min = 1.0,
    algorithm="direct").f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  # get fmodel
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = xrs_poor,
    sf_and_grads_accuracy_params = params,
    target_name                  = "ls_wunit_kunit")
  # refinement loop
  print("start r_factor: %6.4f" % fmodel.r_work())
  total_call_count = 0
  for macro_cycle in range(10):
    minimized = minimizer(fmodel = fmodel)
    total_call_count += minimized.call_count
    print("  macro_cycle %3d (adp)   r_factor: %6.4f  call_count=%d" % \
    (macro_cycle, fmodel.r_work(), minimized.call_count))
  print('total_call_count =', total_call_count)
  if(1):
    inp.pdb_hierarchy.adopt_xray_structure(fmodel.xray_structure)
    inp.pdb_hierarchy.write_pdb_file(file_name="refined.pdb")

if (__name__ == "__main__"):
  run()
