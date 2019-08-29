"""- Generates random structures
   - Computes structure factors
   - Randomizes phases given a fudge factor
   - Computes fft map
   - Determines skewness of map
"""
from __future__ import absolute_import, division, print_function

from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
from six.moves import range

def randomize_phases(f_calc, fudge_factor):
  assert 0 <= fudge_factor <= 1
  phases = flex.arg(f_calc.data(), True)
  centric_flags = f_calc.centric_flags().data()
  acentric_flags = ~centric_flags
  centric_phases = phases.select(centric_flags)
  acentric_phases = phases.select(acentric_flags)
  sel = flex.random_double(size=centric_phases.size()) < (0.5 * fudge_factor)
  centric_phases.set_selected(sel, centric_phases.select(sel) + 180)
  acentric_phases += (flex.random_double(size=acentric_phases.size())
                      * 360 - 180) * fudge_factor
  phases.set_selected(centric_flags, centric_phases)
  phases.set_selected(acentric_flags, acentric_phases)
  return f_calc.phase_transfer(phases, deg=True)

def skewness_calculation(space_group_info, n_test_points=10,
                         n_sites=20, d_min=3, volume_per_atom=200):
  structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=["Se"]*n_sites,
    volume_per_atom=volume_per_atom,
    random_u_iso=True)
  structure.show_summary()
  print()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  f_calc.show_summary()
  print()
  for i_fudge_factor in range(n_test_points+1):
    fudge_factor = i_fudge_factor/float(n_test_points)
    randomized_f_calc = randomize_phases(f_calc, fudge_factor)
    mwpe = f_calc.mean_weighted_phase_error(randomized_f_calc)
    rho = randomized_f_calc.fft_map().real_map_unpadded()
    # <(rho-rho_bar)**3>/<(rho-rho_bar)**2>**3/2
    rho_rho_bar = rho - flex.mean(rho)
    num = flex.mean(flex.pow(rho_rho_bar, 3))
    den = flex.mean(flex.pow(rho_rho_bar, 2))**(3/2.)
    assert den != 0
    skewness = num / den
    print("fudge factor, phase difference, map skewness:", end=' ')
    print("%4.2f, %5.2f, %.4g" % (fudge_factor, mwpe, skewness))
  print()

def run():
  for space_group_symbol in ("P 1", "C 2", "P 21 21 21", "R 32", "F 4 3 2"):
    skewness_calculation(sgtbx.space_group_info(space_group_symbol))

if (__name__ == "__main__"):
  run()
