from __future__ import absolute_import, division, print_function
import mmtbx.refinement.anomalous_scatterer_groups
import mmtbx.f_model
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import xray
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import random
import sys
from six.moves import zip

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def run_call_back(flags, space_group_info):
  d_min = 2.0
  structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=["N", "C", "O", "S"]*3 + ["Fe"]*2,
    volume_per_atom=100)
  if (not space_group_info.group().is_centric()):
    fp_fdp_targets = [(-1,2), (-2,6)]
  else:
    fp_fdp_targets = [(-1,0), (-2,0)]
  anomalous_scatterer_groups = [
    xray.anomalous_scatterer_group(
      iselection=flex.size_t(),
      f_prime=fp,
      f_double_prime=fdp,
      refine=["f_prime", "f_double_prime"]) for fp,fdp in fp_fdp_targets]
  for i_seq,scatterer in enumerate(structure.scatterers()):
    if (scatterer.scattering_type == "S"):
      anomalous_scatterer_groups[0].iselection.append(i_seq)
    if (scatterer.scattering_type == "Fe"):
      anomalous_scatterer_groups[1].iselection.append(i_seq)
  for group in anomalous_scatterer_groups:
    group.copy_to_scatterers_in_place(scatterers=structure.scatterers())
  if (flags.Verbose):
    structure.show_summary().show_scatterers()
  f_obs = abs(structure.structure_factors(
    d_min=2.0, anomalous_flag=True).f_calc())
  if (flags.Verbose):
    f_obs.show_comprehensive_summary()
  #
  for group in anomalous_scatterer_groups:
    group.f_prime = 0
    group.f_double_prime = 0
    group.copy_to_scatterers_in_place(scatterers=structure.scatterers())
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    xray_structure=structure,
    f_obs=f_obs,
    r_free_flags=f_obs.generate_r_free_flags(),
    sf_and_grads_accuracy_params = sfg_params,
    target_name="ls")
  #
  n_cycles = [0]
  def call_back(minimizer):
    n_cycles[0] += 1
    return True
  minimized = mmtbx.refinement.anomalous_scatterer_groups.minimizer(
    fmodel=fmodel,
    groups=anomalous_scatterer_groups,
    call_back_after_minimizer_cycle=call_back,
    number_of_finite_difference_tests=3)
  assert n_cycles == [3]
  #
  for group,(fp,fdp) in zip(anomalous_scatterer_groups, fp_fdp_targets):
    # Large eps because the minimization doesn't reliably converge.
    # We don't want to exercise the minimizer here, the important
    # test is the finite difference test embedded in the minimizer.
    assert approx_equal(group.f_prime, fp, eps=1)
    assert approx_equal(group.f_double_prime, fdp, eps=1)

def exercise():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  exercise()
