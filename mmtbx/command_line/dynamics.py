# LIBTBX_SET_DISPATCHER_NAME phenix.dynamics

from __future__ import division
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
cartesian_dynamics
  .short_caption = Cartesian dynamics
{
  include scope mmtbx.dynamics.cartesian_dynamics.master_params
}""" % geometry_minimization.base_params_str

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def run_cartesian_dynamics (
    xray_structure,
    restraints_manager,
    params,
    log) :
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
    n_steps=params.number_of_steps,
    time_step=params.time_step,
    initial_velocities_zero_fraction=params.initial_velocities_zero_fraction,
    n_print=params.n_print,
    log=log,
    verbose=1)
  sites_cart_end = xray_structure.sites_cart()
  rmsd = sites_cart_end.rms_difference(sites_cart_start)
  print >> log, ""
  print >> log, "RMSD from starting structure: %.3f" % rmsd

class run (geometry_minimization.run) :
  _pdb_suffix = "shaken"
  def master_params (self) :
    return master_params()

  def format_usage_message (self) :
    raise Usage("""\
phenix.dynamics: perform simple dynamics to perturb a model
Usage examples:
  phenix.dynamics model.pdb
  phenix.dynamics model.pdb ligands.cif
""")

  def __execute (self) :
    #
    self.caller(func = self.initialize,     prefix="Initialization, inputs")
    self.caller(func = self.process_inputs, prefix="Processing inputs")
    self.caller(func = self.atom_selection, prefix="Atom selection")
    self.caller(func = self.get_restraints, prefix="Geometry Restraints")
    self.caller(func = self.dynamics,       prefix="Dynamics")
    self.caller(func = self.write_pdb_file, prefix="Write PDB file")
    self.caller(func = self.write_geo_file, prefix="Write GEO file")
    #
    self.show_times()

  def atom_selection (self, prefix) :
    self.selection = flex.bool(self.xray_structure.sites_cart().size(), True)
    geometry_minimization.broadcast(m=prefix, log = self.log)
    self.generate_restrain_selection()

  def dynamics (self, prefix) :
    geometry_minimization.broadcast(m=prefix, log = self.log)
    self.sites_cart = self.xray_structure.sites_cart()
    if (self.params.dynamics_type == "cartesian") :
      run_cartesian_dynamics(
        xray_structure=self.xray_structure,
        restraints_manager=self.grm,
        params=self.params.cartesian_dynamics,
        log=self.log)
    else : # TODO
      raise NotImplementedError()
    self.sites_cart = self.xray_structure.sites_cart()

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=self.args, log=sys.stdout,
      use_directory_prefix=False).output_file_name

def validate_params (params) :
  return geometry_minimization.validate_params(params)

def finish_job (result) :
  return geometry_minimization.finish_job(result)

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  log = sys.stdout
  o = run(sys.argv[1:], log=log)
  tt = timer.elapsed()
  print >> o.log, "Overall runtime: %-8.3f" % tt
  assert abs(tt-o.total_time) < 0.1 # guard against unaccounted times
