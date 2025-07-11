from __future__ import division, print_function
import time, os
import numpy as np
import iotbx.pdb
import mmtbx.model
import libtbx.load_env
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from iotbx.data_manager import DataManager
import pydiscamb

#def set_gradient_flags(scatterers, scf):
#  scatterers.flags_set_grads(state=False)
#  scf(iselection = iselection)

phil_str = '''
data_manager {
    model {
      file = %s
      type = x_ray neutron *electron reference
    }
    miller_array {
      file = %s
      labels {
        name = %s
        type = x_ray neutron *electron
      }
      labels {
        name = %s
        type = x_ray neutron *electron
      }
  }
}'''

def run():
  pdb_fn = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/7n2l.pdb",
    test=os.path.isfile)
  mtz_fn = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/mtz/7n2l.mtz",
    test=os.path.isfile)
  label_fobs = 'FOBS,SIGFOBS'
  label_free = 'R-free-flags'
  dm_phil_str = phil_str % (pdb_fn, mtz_fn, label_fobs, label_free)
  dm_phil = iotbx.phil.parse(dm_phil_str)
  dm = DataManager(['model', 'miller_array'])
  dm.load_phil_scope(dm_phil)
  # get fmodel from data manager
  fmodel = dm.get_fmodel(scattering_table="electron")
  fmodel.update_all_scales()

  m = dm.get_model()
  isel_iso = m.selection('element H').iselection()
  isel_aniso = m.selection('not element H').iselection()

  #print("r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(), fmodel.r_free()))
  #fmodel.show_short(show_k_mask=True, log=None, prefix="")

  xrs = fmodel.xray_structure
  # Make sure we use the correct scattering table
  assert(xrs.get_scattering_table() == 'electron')
  #
  scatterers = xrs.scatterers()
  fmodel.sfg_params.algorithm = 'direct' # don't forget this!
  target = fmodel.target_functor()(compute_gradients=True)

  # Get target derivatives
  d_target_d_fcalc = target.d_target_d_f_calc_work()

  # Calculate derivatives with discamb
  w = pydiscamb.DiscambWrapper(xrs)
  w.set_indices(d_target_d_fcalc.indices())
  discamb_result = w.d_target_d_params(list(d_target_d_fcalc.data()))

  iselection = flex.bool(scatterers.size(), True).iselection()

  # coordinates
  scatterers.flags_set_grads(state=False)
  scatterers.flags_set_grad_site(iselection=iselection)
  expected = flex.vec3_double(target.gradients_wrt_atomic_parameters().packed())
  actual = flex.vec3_double([res.site_derivatives for res in discamb_result])
  assert(approx_equal(expected,actual, eps=1e-6))

  # isotropic ADP
  scatterers.flags_set_grads(state=False)
  scatterers.flags_set_grad_u_iso(iselection=isel_iso)
  expected = target.gradients_wrt_atomic_parameters().packed()
  #actual = np.array([res.adp_derivatives for res in discamb_result]).flatten()
  _actual = [res.adp_derivatives for res in discamb_result]
  actual = np.array([_actual[i] for i in isel_iso]).flatten()
  assert(approx_equal(expected,actual, eps=1e-6))

  # anisotropic ADP
  scatterers.flags_set_grads(state=False)
  scatterers.flags_set_grad_u_aniso(iselection=isel_aniso)
  expected = target.gradients_wrt_atomic_parameters().packed()
  _actual = [res.adp_derivatives for res in discamb_result]
  actual = np.array([_actual[i] for i in isel_aniso]).flatten()
  assert(approx_equal(expected,actual, eps=1e-6))

  # Occupancy
  scatterers.flags_set_grads(state=False)
  scatterers.flags_set_grad_occupancy(iselection=iselection)
  expected = target.gradients_wrt_atomic_parameters().packed()
  actual = [res.occupancy_derivatives for res in discamb_result]
  assert(approx_equal(expected,actual, eps=1e-6))

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
