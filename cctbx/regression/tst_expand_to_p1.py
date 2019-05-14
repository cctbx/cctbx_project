from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx.development import debug_utils
import scitbx.math
import sys

def exercise(
      space_group_info,
      use_primitive_setting,
      anomalous_flag,
      use_u_aniso,
      n_elements=3,
      d_min=3.,
      verbose=0):
  if (use_primitive_setting):
    space_group_info = space_group_info.primitive_setting()
  structure = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=200,
    random_f_prime_d_min=1.0,
    random_f_double_prime=anomalous_flag,
    random_u_iso=True,
    use_u_aniso=use_u_aniso,
    random_occupancy=True)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
    print("n_special_positions:", \
          structure.special_position_indices().size())
  structure_p1 = structure.expand_to_p1()
  assert structure_p1.scatterers()[0].label == "Se1"
  assert structure_p1.scatterers()[-1].label == ("Se%d" % n_elements)
  structure_p1 = structure.expand_to_p1(append_number_to_labels=True)
  if (0 or verbose):
    structure_p1.show_summary().show_scatterers()
  l,i = structure_p1.scatterers()[0].label.split("_")
  assert l == "Se1"
  assert int(i) == 0
  l,i = structure_p1.scatterers()[-1].label.split("_")
  assert l == ("Se%d" % n_elements)
  assert int(i) \
      == structure.site_symmetry_table().get(n_elements-1).multiplicity()-1
  assert structure_p1.special_position_indices().size() == 0
  f_calc = structure.structure_factors(
    anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct").f_calc()
  miller_set_p1 = miller.set.expand_to_p1(f_calc)
  f_calc_p1 = f_calc.expand_to_p1()
  assert flex.order(miller_set_p1.indices(), f_calc_p1.indices()) == 0
  amplitudes = abs(f_calc)
  amplitudes_p1 = amplitudes.expand_to_p1()
  assert flex.order(miller_set_p1.indices(), amplitudes_p1.indices()) == 0
  c = flex.linear_correlation(abs(f_calc_p1).data(), amplitudes_p1.data())
  assert c.is_well_defined()
  assert c.n() > 20
  if (0 or verbose):
    print("correlation:", c.coefficient())
  assert c.coefficient() > 0.999
  for phase_deg in (False, True):
    phases = f_calc.arg(phase_deg)
    phases_p1 = phases.expand_to_p1(phase_deg)
    assert flex.order(miller_set_p1.indices(), phases_p1.indices()) == 0
    f_calc_p1_phases = f_calc_p1.arg(phase_deg)
    for i,phase in enumerate(f_calc_p1_phases.data()):
      e = scitbx.math.phase_error(phase, phases_p1.data()[i], deg=phase_deg)
      assert e < 1.e-6
  ctrl_amplitudes_p1 = abs(miller_set_p1.structure_factors_from_scatterers(
    xray_structure=structure_p1,
    algorithm="direct").f_calc())
  c = flex.linear_correlation(amplitudes_p1.data(), ctrl_amplitudes_p1.data())
  assert c.is_well_defined()
  if (0 or verbose):
    print("correlation:", c.coefficient())
  assert c.coefficient() > 0.999

def run_call_back(flags, space_group_info):
  use_primitive_setting_flags = [False]
  if (space_group_info.group().conventional_centring_type_symbol() != "P"):
    use_primitive_setting_flags.append(True)
  for use_primitive_setting in use_primitive_setting_flags:
    for anomalous_flag in (False, True)[:]: #SWITCH
      for use_u_aniso in (False, True)[:]: #SWITCH
        exercise(
          space_group_info=space_group_info,
          use_primitive_setting=use_primitive_setting,
          anomalous_flag=anomalous_flag,
          use_u_aniso=use_u_aniso,
          verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
