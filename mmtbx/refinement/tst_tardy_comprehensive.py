from __future__ import division
from mmtbx.refinement import tst_tardy_pdb
import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import show_times_at_exit
from libtbx.queuing_system_utils import chunk_manager
from libtbx import Auto
from cStringIO import StringIO
import pprint
import traceback
import sys, os
op = os.path

def report_exception(context_info):
  print ">Begin exception"
  print "Exception:", context_info
  sys.stdout.flush()
  sys.stderr.flush()
  traceback.print_exc()
  print ">End exception"
  sys.stdout.flush()
  sys.stderr.flush()
  print

class collector(object):

  def __init__(O):
    O.tardy_model = None
    O.rmsd_calculator = None
    O.rmsd = flex.double()

  def __call__(O, tardy_model=None, rmsd_calculator=None):
    if (tardy_model is not None):
      assert rmsd_calculator is not None
      assert O.tardy_model is None
      O.tardy_model = tardy_model
      O.rmsd_calculator = rmsd_calculator
    else:
      assert O.tardy_model is not None
    O.rmsd.append(O.rmsd_calculator(
      O.tardy_model.potential_obj.ideal_sites_cart,
      O.tardy_model.sites_moved()))

common_parameter_trial_table = [
  ("tardy_displacements_auto.rmsd_vs_high_resolution_factor", (1/3, 2/3, 1)),
  ("structure_factors_high_resolution", (1.25, 2.5, 3.75, 5)),
  ("real_space_target_weight", (10, 100, 500)),
  ("real_space_gradients_delta_resolution_factor", (1/3,)),
  ("emulate_cartesian", (False, True))
]
annealing_parameter_trial_table = common_parameter_trial_table + [
  ("start_temperature_kelvin", (2500, 5000)),
  ("number_of_cooling_steps", (250, 500))
]

def number_of_trials(table):
  result = 1
  for name,values in table:
    result *= len(values)
  return result

def set_parameters(params, trial_table, cp_i_trial):
  rest = cp_i_trial
  for i in reversed(range(len(trial_table))):
    name, values = trial_table[i]
    n = len(values)
    j = rest % n
    rest //= n
    phil_path = name.split(".")
    assert len(phil_path) > 0
    scope = params
    for scope_name in phil_path[:-1]:
      scope = getattr(scope, scope_name)
    setattr(scope, phil_path[-1], values[j])
  assert rest == 0

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""\
pdb_file = None
  .type = path
algorithm = *minimization annealing
  .type = choice
  .optional = False
orca_experiments = False
  .type = bool
random_displacements_parameterization = *constrained cartesian
  .type = choice
  .optional = False
%(dihedral_function_type_params_str)s
number_of_random_trials = 2
  .type = int
hot = False
  .type = bool
verbose = False
  .type = bool
keep_going = False
  .type = bool
chunk = 1 0
  .type = ints(size=2, value_min=0)
""" % tst_tardy_pdb.pdb_interpretation.__dict__)

def run(args):
  local_master_phil = get_master_phil()
  argument_interpreter = local_master_phil.command_line_argument_interpreter()
  phil_objects = []
  for arg in args:
    phil_objects.append(argument_interpreter.process(arg=arg))
  local_params = local_master_phil.fetch(sources=phil_objects).extract()
  chunk = chunk_manager(
    n=local_params.chunk[0],
    i=local_params.chunk[1]).easy_all()
  local_master_phil.format(local_params).show()
  print
  #
  assert local_params.pdb_file is not None
  assert op.isfile(local_params.pdb_file)
  #
  tst_tardy_pdb_master_phil = tst_tardy_pdb.get_master_phil()
  tst_tardy_pdb_params = tst_tardy_pdb_master_phil.extract()
  tst_tardy_pdb_params.tardy_displacements = Auto
  tst_tardy_pdb_params.tardy_displacements_auto.parameterization \
    = local_params.random_displacements_parameterization
  if (local_params.algorithm == "minimization"):
    parameter_trial_table = common_parameter_trial_table
  elif (local_params.algorithm == "annealing"):
    parameter_trial_table = annealing_parameter_trial_table
  else:
    raise AssertionError
  cp_n_trials = number_of_trials(table=parameter_trial_table)
  print "Number of parameter trials:", cp_n_trials
  print "parameter_trial_table:"
  pprint.pprint(parameter_trial_table)
  print
  #
  show_times_at_exit()
  #
  params_shown_once_already = False
  tst_tardy_pdb_log_shown_once_already = False
  for cp_i_trial in xrange(cp_n_trials):
    if (chunk.skip_iteration(i=cp_i_trial)): continue
    print "cp_i_trial: %d / %d = %.2f %%" % (
      cp_i_trial, cp_n_trials, 100 * (cp_i_trial+1) / cp_n_trials)
    if (local_params.verbose):
      print
    sys.stdout.flush()
    set_parameters(
      params=tst_tardy_pdb_params,
      trial_table=parameter_trial_table,
      cp_i_trial=cp_i_trial)
    if (local_params.algorithm == "minimization"):
      if (local_params.orca_experiments):
        tst_tardy_pdb_params.keep_all_restraints = True
        if (tst_tardy_pdb_params.emulate_cartesian):
          tst_tardy_pdb_params.orca_experiments = False
        else:
          tst_tardy_pdb_params.orca_experiments = True
      tst_tardy_pdb_params.number_of_cooling_steps = 0
      tst_tardy_pdb_params.minimization_max_iterations = None
    elif (local_params.algorithm == "annealing"):
      tst_tardy_pdb_params.number_of_time_steps = 1
      tst_tardy_pdb_params.time_step_pico_seconds = 0.001
      tst_tardy_pdb_params.minimization_max_iterations = 0
    else:
      raise AssertionError
    for random_seed in xrange(local_params.number_of_random_trials):
      tst_tardy_pdb_params.random_seed = random_seed
      tst_tardy_pdb_params.dihedral_function_type \
        = local_params.dihedral_function_type
      if (local_params.verbose or not params_shown_once_already):
        params_shown_once_already = True
        tst_tardy_pdb_master_phil.format(tst_tardy_pdb_params).show()
        print
        sys.stdout.flush()
      if (local_params.hot):
        if (local_params.verbose):
          tst_tardy_pdb_log = sys.stdout
        else:
          tst_tardy_pdb_log = StringIO()
        coll = collector()
        try:
          tst_tardy_pdb.run_test(
            params=tst_tardy_pdb_params,
            pdb_files=[local_params.pdb_file],
            other_files=[],
            callback=coll,
            log=tst_tardy_pdb_log)
        except KeyboardInterrupt: raise
        except:
          print
          print "tst_tardy_pdb_params leading to exception:"
          print
          tst_tardy_pdb_master_phil.format(tst_tardy_pdb_params).show()
          print
          if (not local_params.verbose):
            sys.stdout.write(tst_tardy_pdb_log.getvalue())
          sys.stdout.flush()
          if (not local_params.keep_going):
            raise
          report_exception(
            context_info="cp_i_trial=%d, random_seed=%d" % (
              cp_i_trial, random_seed))
        else:
          if (    not local_params.verbose
              and not tst_tardy_pdb_log_shown_once_already):
            tst_tardy_pdb_log_shown_once_already = True
            sys.stdout.write(tst_tardy_pdb_log.getvalue())
          print "RESULT_cp_i_trial_random_seed_rmsd:", \
            cp_i_trial, random_seed, list(coll.rmsd)
          sys.stdout.flush()
      first_pass = False
    if (local_params.hot):
      print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
