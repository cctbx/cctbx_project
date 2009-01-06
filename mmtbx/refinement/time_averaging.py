from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys, copy
from libtbx.test_utils import approx_equal
from cctbx import xray
from mmtbx.dynamics import cartesian_dynamics

def run(fmodels, model, target_weights, log):
  tx = 1.
  f_calc_average = None
  n_cycles = 5
  if 0: fmodels.fmodel_xray().xray_structure.set_b_iso(value=1.0)
  for i_cycle in xrange(n_cycles):
    print >> log, "Cycle number: ", i_cycle
    print >> log, "start R=%6.4f" % fmodels.fmodel_xray().r_work()
    xrs_start = fmodels.fmodel_xray().xray_structure.deep_copy_scatterers()
    cd_manager = cartesian_dynamics.cartesian_dynamics(
      structure                   = fmodels.fmodel_xray().xray_structure,
      restraints_manager          = model.restraints_manager,
      temperature                 = 300.,
      n_steps                     = 200,
      time_step                   = 0.0005,
      fmodel                      = fmodels.fmodel_xray(),
      shift_update                = None,
      xray_target_weight          = target_weights.xyz_weights_result.wx * target_weights.xyz_weights_result.wx_scale,
      chem_target_weight          = 1,
      xray_structure_last_updated = model.xray_structure,
      log = log)
    xrs_final = cd_manager.xray_structure_last_updated
    result = xrs_start.distances(other = xrs_final)
    fmodels.update_xray_structure(
      xray_structure = cd_manager.xray_structure_last_updated,
      update_f_calc  = True)
    print >> log, "final R=%6.4f" % fmodels.fmodel_xray().r_work(), \
      xrs_start.distances(other = xrs_final).min_max_mean().as_tuple()
    f_calc = fmodels.fmodel_xray().f_calc()
    if(f_calc_average is None):
      f_calc_average = f_calc
    else:
      a_prime = math.exp(-cd_manager.n_steps*n_cycles/tx)
      f_calc_data = f_calc.data()
      f_calc_average_data = f_calc_average.data()
      f_calc_average_data = a_prime * f_calc_average_data + (1.-a_prime) * f_calc_data
      f_calc_average = f_calc_average.array(data = f_calc_average_data)
    fmodels.fmodel_xray().update(f_calc = f_calc_average)
    print >> log, "final R=%6.4f" % fmodels.fmodel_xray().r_work()
