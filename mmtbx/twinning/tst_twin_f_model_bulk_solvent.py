from __future__ import absolute_import, division, print_function
"""
The twin model manager keeps a flat bulk-solvent and scale model (k_sol,
b_sol, b_cart, k_overall) for the twin target, while its R-factors come from
a regular fmodel with per-bin k_mask and per-reflection scales. The two must
describe the same model. With the second mask below the regular scaling of
this twinned 4 A data set reports k_sol=0, b_cart=0 and folds the overall
scale into k_anisotropic; taken literally that gave a twin target model with
R=28 and an ADP refinement that blew up (phenix_regression tst_091, ADP
weight 0.0008 instead of 2.3). With the first mask the dropped anisotropic
scale alone cost 0.035 in R.
"""
import os
import libtbx.load_env
import iotbx.pdb
import mmtbx.model
import mmtbx.masks
from mmtbx.twinning import twin_f_model
from cctbx import sgtbx
from iotbx import reflection_file_reader
from libtbx.utils import null_out

def run():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/refinement/data/tst_091.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/refinement/data/tst_091.mtz",
    test=os.path.isfile)
  if (None in [pdb_file, mtz_file]):
    print("phenix_regression not found, skipping")
    return
  f_obs = flags = None
  for a in reflection_file_reader.any_reflection_file(
      mtz_file).as_miller_arrays():
    if(a.info().label_string().startswith("FreeRflag")): flags = a
    if(a.info().label_string().startswith("I-obs")): f_obs = a
  f_obs = f_obs.f_sq_as_f().resolution_filter(d_min=4.0)
  flags = flags.customized_copy(data=(flags.data()==0))
  f_obs, flags = f_obs.common_sets(flags)
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(pdb_file),
    log=null_out())
  model.setup_scattering_dictionaries(scattering_table="n_gaussian")
  for r_solv, r_shrink in [(1.1, 0.9), (0.95, 0.7)]:
    mp = mmtbx.masks.mask_master_params.extract()
    mp.solvent_radius = r_solv
    mp.shrink_truncation_radius = r_shrink
    tm = twin_f_model.twin_model_manager(
      f_obs          = f_obs,
      r_free_flags   = flags,
      xray_structure = model.get_xray_structure(),
      twin_law       = sgtbx.rt_mx("h,-k,-l"),
      mask_params    = mp,
      out            = null_out())
    tm.update_all_scales(log=null_out())
    r_fmodel = tm.r_work()                      # regular fmodel, per-bin scales
    r_target, r_target_free = tm.r_values(table=False)  # twin target model
    print("mask %.2f/%.2f  r_work: fmodel %.4f twin target model %.4f  "
      "k_sol %.3f b_sol %.1f" % (r_solv, r_shrink, r_fmodel, r_target,
      tm.k_sol(), tm.b_sol()))
    assert abs(r_target - r_fmodel) < 0.01, (r_target, r_fmodel)
    assert abs(r_target_free - tm.r_free()) < 0.01, (r_target_free, tm.r_free())
    assert tm.k_sol() > 0.2, tm.k_sol()

if (__name__ == "__main__"):
  run()
  print("OK")
