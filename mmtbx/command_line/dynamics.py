"""Shake up structure with simple molecular dynamics"""
# LIBTBX_SET_DISPATCHER_NAME phenix.dynamics

from __future__ import absolute_import, division, print_function
from mmtbx.command_line import geometry_minimization
import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import user_plus_sys_time, Usage
from libtbx import runtime_utils
from libtbx.str_utils import make_header
import os
import sys

master_params_str = """
%s
dynamics_type = *cartesian
  .type = choice
stop_at_diff = None
  .type = float
  .help = stop after reaching specified cutoff value
cartesian_dynamics
  .short_caption = Cartesian dynamics
{
  include scope mmtbx.dynamics.cartesian_dynamics.master_params
}""" % geometry_minimization.base_params_str

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def run_cartesian_dynamics(
    xray_structure,
    states_collector,
    restraints_manager,
    params,
    stop_at_diff,
    log):
  from mmtbx.dynamics import cartesian_dynamics
  make_header("Simple cartesian dynamics", out=log)
  sites_cart_start = xray_structure.sites_cart().deep_copy()
  gradients_calculator = \
    cartesian_dynamics.gradients_calculator_reciprocal_space(
      restraints_manager = restraints_manager,
      sites_cart         = xray_structure.sites_cart(),
      wc                 = 1)
  cartesian_dynamics.run(
    xray_structure=xray_structure,
    gradients_calculator=gradients_calculator,
    temperature=params.temperature,
    states_collector=states_collector,
    n_steps=params.number_of_steps,
    time_step=params.time_step,
    initial_velocities_zero_fraction=params.initial_velocities_zero_fraction,
    n_print=params.n_print,
    stop_cm_motion=params.stop_cm_motion,
    stop_at_diff=stop_at_diff,
    log=log,
    verbose=1)
  sites_cart_end = xray_structure.sites_cart()
  rmsd = sites_cart_end.rms_difference(sites_cart_start)
  print("", file=log)
  print("RMSD from starting structure: %.3f" % rmsd, file=log)
  return sites_cart_end

class run(geometry_minimization.run):
  _pdb_suffix = "shaken"
  def master_params(self):
    return master_params()

  def format_usage_message(self):
    raise Usage("""\
phenix.dynamics: perform simple dynamics to perturb a model
Usage examples:
  phenix.dynamics model.pdb
  phenix.dynamics model.pdb ligands.cif
""")

  def __execute(self):
    #
    self.caller(self.initialize,           "Initialization, inputs")
    self.caller(self.process_inputs,       "Processing inputs")
    self.caller(self.atom_selection,       "Atom selection")
    self.caller(self.get_restraints,       "Geometry Restraints")
    self.caller(self.dynamics,             "Dynamics")
    self.caller(self.write_pdb_file,       "Write PDB file")
    self.caller(self.write_geo_file,       "Write GEO file")
    #
    self.show_times()

  def atom_selection(self, prefix):
    self.selection = flex.bool(self.model.get_xray_structure().sites_cart().size(), True)
    geometry_minimization.broadcast(m=prefix, log = self.log)

  def dynamics(self, prefix):
    geometry_minimization.broadcast(m=prefix, log = self.log)
    self.sites_cart = self.model.get_sites_cart()
    if (self.params.dynamics_type == "cartesian"):
      sites_cart_result = run_cartesian_dynamics(
        xray_structure=self.model.get_xray_structure(),
        restraints_manager=self.model.get_restraints_manager(),
        states_collector=self.states_collector,
        params=self.params.cartesian_dynamics,
        stop_at_diff=self.params.stop_at_diff,
        log=self.log)
    else : # TODO
      raise NotImplementedError()
    self.model.set_sites_cart(sites_cart=sites_cart_result)

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    filename = run(args=self.args, log=sys.stdout,
                   use_directory_prefix=False).result_model_fname
    return os.path.join(self.output_dir, filename)

def validate_params(params):
  return geometry_minimization.validate_params(params)

def finish_job(result):
  return geometry_minimization.finish_job(result)

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  log = sys.stdout
  o = run(sys.argv[1:], log=log)
  tt = timer.elapsed()
  print("Overall runtime: %-8.3f" % tt, file=o.log)
  assert abs(tt-o.total_time) < 0.1 # guard against unaccounted times

