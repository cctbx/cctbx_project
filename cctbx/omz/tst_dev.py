from __future__ import division
from cctbx.omz import dev
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import libtbx.phil.command_line
import random
import sys

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def run_call_back(flags, space_group_info, params):
  structure_shake = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "O"),
    volume_per_atom=200,
    min_distance=2.0,
    general_positions_only=params.general_positions_only,
    random_u_iso=True)
  structure_ideal = structure_shake.deep_copy_scatterers()
  structure_shake.shake_sites_in_place(rms_difference=params.shake_sites_rmsd)
  structure_shake.shake_adp(spread=params.shake_adp_spread)
  #
  r1 = dev.run_refinement(
    structure_ideal=structure_ideal,
    structure_shake=structure_shake,
    params=params).r1_factor()
  print "R1: %.4f" % r1
  print

def run(args):
  master_phil = dev.get_master_phil()
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil)
  phil_objects = []
  remaining_args = []
  for arg in args:
    if (arg.find("=") >= 0):
      phil_objects.append(argument_interpreter.process(arg=arg))
    else:
      remaining_args.append(arg)
  work_phil = master_phil.fetch(sources=phil_objects)
  work_phil.show()
  print
  params = work_phil.extract()
  debug_utils.parse_options_loop_space_groups(
    argv=remaining_args, call_back=run_call_back, params=params)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
