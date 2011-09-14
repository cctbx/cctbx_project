from cctbx import xray
from cctbx.array_family import flex
import sys, random
import maxlik
from cctbx.development import random_structure
from cctbx.development import debug_utils


def exercise(space_group_info,
             n_sites=100,
             d_min=2.0,
             volume_per_atom=50,
             test_fraction = 0.3):

 #frm = (0.01,0.05,0.1,0.5,0.6) #,0.7,0.8,0.9)
 frm = (0.01,0.05,0.1,0.5,0.6)
 for fraction_missing in frm:
   print ">>> FRACTION MISSING = ", fraction_missing

   structure = random_structure.xray_structure(
      space_group_info=space_group_info,
      elements=["C"]*n_sites,
      volume_per_atom=volume_per_atom,
      random_u_iso=True)
   f_calc = structure.structure_factors(d_min=d_min,
                                anomalous_flag=False).f_calc()

   n_keep = int(round(structure.scatterers().size() * (1-fraction_missing)))
   partial_structure = xray.structure(special_position_settings=structure)
   partial_structure.add_scatterers(structure.scatterers()[:n_keep])
   f_calc_partial = partial_structure.structure_factors(d_min=d_min,
                                                anomalous_flag=False).f_calc()

   f_obs = abs(f_calc)
   f_calc= abs(f_calc_partial)

   flags=flex.bool(f_calc_partial.indices().size(), True)
   for i in xrange(f_calc_partial.indices().size()):
     if (random.random() < (1.0 - test_fraction)):
       flags[i]=False
     else:
       flags[i]=True

   if(1):
     alpha, beta = maxlik.alpha_beta_est(f_obs  = f_obs,
                                      f_calc = f_calc,
                                      free_reflections_per_bin = 50,
                                      flags = flags).alpha_beta(
                                               show_summary=True)

     f_star, wstar= maxlik.f_star_w_star(
                      f_obs = f_obs,
                      alpha = alpha,
                      beta  = beta).fs_ws(show_summary=True)

     assert flex.min(f_star.data()) >= 0.0
     assert f_obs.size() == f_star.size()
     assert f_obs.size() >= list(f_star.data()).count(0.0)
     assert len(alpha) == len(beta)
     assert len(alpha) == f_star.size()
     assert min(alpha) >= 0.0
     assert min(beta) >= 0.0
     assert min(wstar) > 0.0 and max(wstar) == 1.0

   if(1):
     alpha, beta = maxlik.alpha_beta_calc(
                    f   = f_obs,
                    n_atoms_absent  = structure.scatterers().size() - n_keep,
                    n_atoms_included= n_keep,
                    bf_atoms_absent = 25.0,
                    final_error     = 0.0).alpha_beta()

     f_star, wstar= maxlik.f_star_w_star(
                      f_obs = f_obs,
                      alpha = alpha,
                      beta  = beta).fs_ws(show_summary=True)

     assert flex.min(f_star.data()) >= 0.0
     assert f_obs.size() == f_star.size()
     assert f_obs.size() >= list(f_star.data()).count(0.0)
     assert len(alpha) == len(beta)
     assert len(alpha) == f_star.size()
     assert min(alpha) >= 0.0
     assert min(beta) >= 0.0
     assert min(wstar) > 0.0 and max(wstar) == 1.0

def run_call_back(flags, space_group_info):
  exercise(space_group_info=space_group_info)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
