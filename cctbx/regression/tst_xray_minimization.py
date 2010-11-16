from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import sys, random

def shift_u_iso(structure, shift):
  for sc in structure.scatterers():
    if(sc.flags.use_u_iso() and sc.flags.grad_u_iso()):
      sc.u_iso += (shift * random.random())

def shift_u_aniso(structure, shift):
  for sc in structure.scatterers():
    if(sc.flags.use_u_aniso() and sc.flags.grad_u_aniso()):
      u_star = list(sc.u_star)
      for i in xrange(6):
        u_star[i] += (shift * random.random())
      sc.u_star = u_star

def exercise(target_functor, data_type, space_group_info, anomalous_flag,
             gradient_flags, occupancy_penalty,
             n_elements=9, d_min=None, shake_sigma=0.1,
             verbose=0,tan_u_iso=False, param = 0):
  assert data_type == 'F' or data_type == 'F^2'
  if (data_type == 'F^2'
      and not target_functor == xray.unified_least_squares_residual): return
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=(("O","N","C")*(n_elements))[:n_elements],#("Se",)*n_elements,
    volume_per_atom=200,
    random_u_iso=True,
    random_u_cart_scale=.3,
    random_occupancy=True,
    use_u_aniso = True)
  random_structure.random_modify_adp_and_adp_flags(
    scatterers         = structure_ideal.scatterers(),
    random_u_iso_scale = 0.3,
    random_u_iso_min   = 0.0)
  xray.set_scatterer_grad_flags(scatterers = structure_ideal.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp,
                                tan_u_iso  = tan_u_iso,
                                param      = param)
  if(0):
    print
    for sc in structure_ideal.scatterers():
      print sc.flags.use_u_iso(),sc.flags.grad_u_iso(),sc.flags.use_u_aniso(),\
            sc.flags.grad_u_aniso(),sc.u_iso, sc.u_star,sc.flags.tan_u_iso(),\
            sc.flags.param, sc.occupancy
  rnd_f_calc = structure_ideal.structure_factors(
    anomalous_flag=anomalous_flag,
    d_min=d_min,
    algorithm="direct",
    cos_sin_table=True).f_calc()
  if data_type == "F":
    y_obs = abs(rnd_f_calc)
  elif data_type == "F^2":
    y_obs = rnd_f_calc.norm()
    y_obs.set_observation_type_xray_intensity()
  else:
    raise "Error: invalid data type: %s" % data_type
  if (0 or verbose):
    print "structure_ideal:"
    structure_ideal.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_ideal.special_position_indices().size()
    print
  structure_ideal_cp = structure_ideal.deep_copy_scatterers()
  structure_shake = structure_ideal
  if (gradient_flags.site):
    structure_shake = structure_shake.random_modify_parameters(
      "site", shake_sigma)
  if (gradient_flags.occupancy):
    structure_shake = structure_shake.random_modify_parameters(
      "occupancy", shake_sigma)
    if (occupancy_penalty is not None):
      structure_shake.scatterers()[-1].occupancy = 0
  if (gradient_flags.u_aniso):
    shift_u_aniso(structure_shake, 0.001)
  if (gradient_flags.u_iso):
    shift_u_iso(structure_shake, 0.1)
  assert tuple(structure_ideal.special_position_indices()) \
         == tuple(structure_shake.special_position_indices())
  if (0 or verbose):
    print "structure_shake:"
    structure_shake.show_summary().show_scatterers()
    print
  for i_trial in xrange(10):
    try:
      minimizer = xray.minimization.lbfgs(
        target_functor=target_functor(y_obs),
        xray_structure=structure_shake,
        occupancy_penalty=occupancy_penalty,
        structure_factor_algorithm="direct")
    except RuntimeError, e:
      if (str(e).find("debye_waller_factor_exp: max_arg exceeded") < 0):
        raise
    else:
      break
  else:
    raise RuntimeError("Too many xray.minimization.lbfgs failures.")
  if (0 or verbose):
    print "first:", minimizer.first_target_value
    print "final:", minimizer.final_target_value
    print
  assert minimizer.final_target_value < minimizer.first_target_value
  if (0 or verbose):
    print "minimized structure_shake:"
    structure_shake.show_summary().show_scatterers()
    print
  f_final = y_obs.structure_factors_from_scatterers(
    xray_structure=structure_shake,
    algorithm="direct",
    cos_sin_table=True).f_calc()
  if data_type == 'F':
    f_final = abs(f_final)
  elif data_type == 'F^2':
    f_final = f_final.norm()
  c = flex.linear_correlation(y_obs.data(), f_final.data())
  assert c.is_well_defined()
  if (0 or verbose):
    label = gradient_flags.string_of_true()
    if (anomalous_flag):
      label += ",anomalous"
    print "correlation: %10.8f" % c.coefficient(), label
    print
  c_coefficient = c.coefficient()
  if(c_coefficient <= 0.999):
    print c_coefficient
  if data_type == 'F':
    assert c_coefficient > 0.999
  elif data_type == 'F^2':
    assert c_coefficient > 0.9

def run_call_back(flags, space_group_info):
  data_type = ('F', 'F^2')[hasattr(flags, 'F_sq')]
  options = (
    ( True,False,False,False),
    (False, True,False,False),
    (False,False, True,False),
    (False,False,False, True),
    (False, True, True, True))
  for target_functor in xray.target_functors.registry().values():
    for (fsite, fu_iso, foccupancy, fu_aniso) in options:
      gradient_flags = xray.structure_factors.gradient_flags(
        site      = fsite,
        u_iso     = fu_iso,
        u_aniso   = fu_aniso,
        occupancy = foccupancy)
      for anomalous_flag in (False, True)[:]: #SWITCH
        u_penalty_types = [None]
        tan_u_isos = [False]
        if (gradient_flags.u_iso):
          tan_u_isos.append(True)
        occupancy_penalty_types = [None]
        if (gradient_flags.occupancy):
          occupancy_penalty_types.append(
            xray.minimization.occupancy_penalty_exp())
        for tan_u_iso in tan_u_isos:
          if(tan_u_iso):
            param = 100
          else:
            param = 0
          for occupancy_penalty in occupancy_penalty_types:
            if(0):
              print fsite,fu_iso,foccupancy,fu_aniso,anomalous_flag,tan_u_iso
            do_exercise = lambda: exercise(
              target_functor,
              data_type,
              space_group_info,
              anomalous_flag,
              gradient_flags,
              occupancy_penalty=occupancy_penalty,
              verbose=flags.Verbose,
              tan_u_iso=tan_u_iso,
              param = param,
              d_min = 2.5)
            try:
              do_exercise()
            except AssertionError:
              print "Test did not pass: ruling out a random fluke..."
              do_exercise()
              do_exercise()
              print "Ruled out!"

def run():
  cmd_args = ['--F', '--F_sq', '--debugging']
  extra = []
  args = []
  for arg in sys.argv[1:]:
    if arg in cmd_args:
      extra.append(arg[2:])
    else:
      args.append(arg)
  extra = tuple(extra)
  debug = ()
  if 'debugging' in extra:
    extra = ()
    if 0:
      extra += ('F_sq',)
  assert not('F' in extra and 'F_sq' in extra)
  if 'F' in extra: print 'Refinement against F'
  if 'F_sq' in extra: print 'Refinement against F^2'
  debug_utils.parse_options_loop_space_groups(
    args, run_call_back, keywords=extra+debug)

if (__name__ == "__main__"):
  run()
