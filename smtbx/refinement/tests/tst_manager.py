from cctbx.development import random_structure
from cctbx import sgtbx
from cctbx import xray
from cctbx import miller
from cctbx.array_family import flex
from smtbx.refinement.manager import manager
from cctbx.eltbx import wavelengths

import sys

def run(args):
  if args:
    # Load structure and reflections from file
    from iotbx.shelx.from_ins import from_ins
    from iotbx.shelx.hklf import fast_reader
    ideal_structure = xray.structure.from_shelx(
      filename=args[0],
      set_grad_flags=True)
    f_sq_obs = fast_reader(filename=args[1])\
             .as_miller_arrays(crystal_symmetry=ideal_structure,
                               merge_equivalents=False)[0]

  else:
    # Create random structure
    space_group_info = sgtbx.space_group_info(
      symbol="P21/c")
    d_min = 1.0
    ideal_structure = random_structure.xray_structure(
      space_group_info = space_group_info,
      elements = ("N", "C", "C", "O") * 5,
      volume_per_atom = 18,
      random_u_iso = True)
    f_calc = abs(ideal_structure.structure_factors(
      d_min=d_min,
      anomalous_flag=False).f_calc())*3
    mt = flex.mersenne_twister(seed=0)
    f = f_calc + (mt.random_double(size=f_calc.indices().size())*3+1)
    sigmas = mt.random_double(size=f_calc.indices().size())*9+1
    f_obs = miller.array(f_calc,
                         data=f.data(),
                         sigmas=sigmas)
    f_obs.set_observation_type_xray_amplitude()
    f_sq_obs = f_obs.f_as_f_sq()

  # Play with the structure
  model_structure = ideal_structure
  #model_structure = ideal_structure.random_shift_sites(max_shift_cart=0.1)
  #model_structure.shake_occupancies()

  lambda_ = wavelengths.characteristic('Mo').as_angstrom()
  refinement = manager(
    f_sq_obs=f_sq_obs,
    xray_structure=model_structure,
    lambda_=lambda_
  )

  refinement.start()
  refinement.show_final_summary()
  peaks = refinement.peak_search()

if (__name__ == "__main__"):
  run(sys.argv[1:])
