import iotbx.pdb
import mmtbx.f_model
from scitbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx

def run():
  xrs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    elements         = ["N"]*500,
    unit_cell        = (20, 30, 40, 70, 80, 120))
  xrs = xrs.set_b_iso(value=25)
  d_mins = [1.5]
  result = []
  for d_min in d_mins:
    f_obs = abs(xrs.structure_factors(d_min=d_min).f_calc())
    f_obs = f_obs.customized_copy(data = f_obs.data() * 135.)
    shifts = [0.0, 0.3]
    for xyz_shake_amount in shifts:
      xrs_shaken = xrs.deep_copy_scatterers()
      xrs_shaken.shake_sites_in_place(mean_distance = xyz_shake_amount)
      ml_err = flex.double()
      ml_err_new = flex.double()
      for trial in xrange(10):
        r_free_flags = f_obs.generate_r_free_flags()
        fmodel = mmtbx.f_model.manager(
          f_obs          = f_obs,
          r_free_flags   = r_free_flags,
          xray_structure = xrs_shaken)
        ml_err_ = fmodel.model_error_ml()
        ml_err.append(ml_err_)
      result.append(flex.mean(ml_err))
  assert result[0] > 0 and result[0] < 0.03
  assert result[1] > 0.2 and result[1] <= 0.32

if (__name__ == "__main__"):
  run()
