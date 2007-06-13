from smtbx import refinement
from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from cctbx import adptbx
import sys, random
from libtbx.test_utils import approx_equal
import scitbx.matrix
import cStringIO

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
             use_special_position_constraints=False,
             n_elements=9, d_min=None, shake_sigma=0.1,
             verbose=0,tan_u_iso=False, param = 0):
  assert data_type == 'F' or data_type == 'F^2'
  if (data_type == 'F^2'
      and not target_functor == xray.unified_least_squares_residual): return
  refined_structure = random_structure.xray_structure(
    space_group_info,
    elements=(("O","N","C")*(n_elements))[:n_elements],#("Se",)*n_elements,
    volume_per_atom=20,
    random_u_iso=True,
    random_u_cart_scale=.3,
    random_occupancy=True,
    use_u_aniso = True)
  sc = refined_structure.scatterers()
  sc[0].flags.set_use_u_iso(True ) ; sc[0].flags.set_use_u_aniso(True )
  sc[1].flags.set_use_u_iso(False) ; sc[1].flags.set_use_u_aniso(False)
  sc[2].flags.set_use_u_iso(True ) ; sc[2].flags.set_use_u_aniso(False)
  sc[3].flags.set_use_u_iso(False) ; sc[3].flags.set_use_u_aniso(True )
  sc[4].flags.set_use_u_iso(False) ; sc[4].flags.set_use_u_aniso(False)
  sc[5].flags.set_use_u_iso(False) ; sc[5].flags.set_use_u_aniso(True )
  sc[6].flags.set_use_u_iso(False) ; sc[6].flags.set_use_u_aniso(False)
  sc[7].flags.set_use_u_iso(True ) ; sc[7].flags.set_use_u_aniso(False)
  sc[8].flags.set_use_u_iso(False) ; sc[8].flags.set_use_u_aniso(False)
  xray.set_scatterer_grad_flags(scatterers = refined_structure.scatterers(),
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
    for sc in refined_structure.scatterers():
      print sc.flags.use_u_iso(),sc.flags.grad_u_iso(),sc.flags.use_u_aniso(),\
          sc.flags.grad_u_aniso(),sc.u_iso, sc.u_star,sc.flags.tan_u_iso(),\
          sc.flags.param, sc.occupancy
  rnd_f_calc = refined_structure.structure_factors(
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

  structure_ideal = refined_structure.deep_copy_scatterers()
  if (0 or verbose):
    print "structure_ideal:"
    structure_ideal.show_summary().show_scatterers()
    print
  if 1:
    print "special/general positions: %i/%i"\
          % (structure_ideal.special_position_indices().size(),
             structure_ideal.scatterers().size()
             - structure_ideal.special_position_indices().size(), )
  if (gradient_flags.site):
    refined_structure = refined_structure.random_modify_parameters(
    "site", shake_sigma)
  if (gradient_flags.occupancy):
    refined_structure = refined_structure.random_modify_parameters(
    "occupancy", shake_sigma)
    if (occupancy_penalty is not None):
      refined_structure.scatterers()[-1].occupancy = 0
  if (gradient_flags.u_aniso):
    shift_u_aniso(refined_structure, 0.001)
  if (gradient_flags.u_iso):
    shift_u_iso(refined_structure, 0.1)
  assert tuple(structure_ideal.special_position_indices()) \
    == tuple(refined_structure.special_position_indices())

  shaken_structure = refined_structure.deep_copy_scatterers()
  if (0 or verbose):
    print "shaken structure:"
    shaken_structure.show_summary().show_scatterers()
    print
  minimizer = refinement.minimization.lbfgs(
    target_functor=target_functor(y_obs),
    xray_structure=refined_structure,
    use_special_position_constraints=use_special_position_constraints,
    occupancy_penalty=occupancy_penalty,
    structure_factor_algorithm="direct",
    log=cStringIO.StringIO())
  if (0 or verbose):
    print "first:", minimizer.first_target_value
    print "final:", minimizer.final_target_value
    print
  assert minimizer.final_target_value < minimizer.first_target_value,\
         "trapped in local minimum"
  if (0 or verbose):
    print "minimized shaken structure:"
    refined_structure.show_summary().show_scatterers()
    print
  f_final = y_obs.structure_factors_from_scatterers(
    xray_structure=refined_structure,
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
    print 'Correlation initial-refined too small: %f' % c_coefficient
  if data_type == 'F':
    assert c_coefficient > 0.999, "bad fit"
  elif data_type == 'F^2':
    assert c_coefficient > 0.9, "bad fit"

def run_call_back(flags, space_group_info):
  print '*****'
  data_type = ('F', 'F^2')[hasattr(flags, 'F_sq')]
  use_spec_pos_constraints = hasattr(flags, 'special_position_constraints')
  options = ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(0,1,1,1))
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
          refinement.barriers.double_quadratic_barrier(lower_bound=0,
                                                       upper_bound=1))
          #occupancy_penalty_types.append(xray.occupancy_penalty_exp())
        for tan_u_iso in tan_u_isos:
          if(tan_u_iso):
            param = 100
          else:
            param = 0
          for occupancy_penalty in occupancy_penalty_types:
            if 0 and fu_aniso: continue
            if 0 and foccupancy: continue
            if 0 and not issubclass(target_functor,
                                    xray.intensity_correlation):
              continue
            if 1 and target_functor is not xray.unified_least_squares_residual:
              continue
            if(1):
              refined = []
              if fsite: refined.append("x's")
              if fu_iso: refined.append("Uiso")
              if foccupancy: refined.append("occ")
              if fu_aniso: refined.append("Uaniso")
              if anomalous_flag: refined.append("f' and f''")
              if tan_u_iso: refined.append("tan Uiso")
              if occupancy_penalty:
                refined.append("with occupancy penalty %s"
                               % occupancy_penalty.__class__.__name__)
              print "Refining %s (%s): %s" % (target_functor.__name__,
                                              data_type,
                                              ", ".join(refined))
            do_exercise = lambda: exercise(
              target_functor=target_functor,
              data_type=data_type,
              space_group_info=space_group_info,
              anomalous_flag=anomalous_flag,
              gradient_flags=gradient_flags,
              use_special_position_constraints=use_spec_pos_constraints,
              occupancy_penalty=occupancy_penalty,
              verbose=flags.Verbose,
              tan_u_iso=tan_u_iso,
              param = param,
              d_min = 0.7)
            try:
              do_exercise()
            except AssertionError, err:
              if 0:
                print "%s: %s" % (err.__doc__, err.args[0])
                continue
              if err.args[0] in ("bad fit", "trapped in local minimum"):
                print 'Ruling out a random fluke ...'
                do_exercise()
                do_exercise()
                print '... Ruled out!'
              else:
                raise

def run():
  cmd_args = ['--special_position_constraints', '--F_sq', '--debugging']
  messages = [ ('special position constraints', 'no','yes'),
               ('data type', 'F', 'F^2') ]
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
  for msg, arg in zip(messages, cmd_args):
    print '%s: %s' % (msg[0], msg[1:][arg[2:] in extra])
  debug_utils.parse_options_loop_space_groups(
    args, run_call_back, keywords=extra+debug)

if (__name__ == "__main__"):
  run()
