from iotbx.option_parser import option_parser
from iotbx import reflection_file_reader
from iotbx.shelx import crystal_symmetry_from_ins
from libtbx.utils import Sorry
from smtbx.ab_initio import charge_flipping

def run(args, command_name='smtbx.charge_flipping'):
  command_line = (option_parser(
    usage="%s [options] reflection_file" % command_name,
    description="Examples:\n"
                "  %s data.mtz\n"
                "  %s hklf4=data.hkl (implicit --symmetry=data.ins)\n"
                "  %s --symmetry=file.ins data.hkl=hklf4\n"
                "  %s --unit_cell='1 2 3 90. 105.1 90.' "
                "--space_group=P21/n hklf3=data.hkl\n"
                % ((command_name,)*4)
    )
    .enable_symmetry_comprehensive()
    .option(None, "--relative_delta",
      action="store",
      type="float",
      default=None,
      help="Parameter delta normalised to the highest reflection amplitude",
      metavar="FLOAT")
    .option(None, "--max_iterations",
      action="store",
      type="int",
      default=100,
      help="Number of iterations to perform",
      metavar="FLOAT")
  ).process(args)
  if not command_line.args:
    command_line.parser.show_help()
    return 1
  for file_name in command_line.args:
    reflection_file = reflection_file_reader.any_reflection_file(file_name)
    miller_arrays = None
    if reflection_file.file_type() is not None:
      try:
        miller_arrays = reflection_file.as_miller_arrays(
          crystal_symmetry=command_line.symmetry,
          force_symmetry=True,
          merge_equivalents=False)
      except Sorry, KeyboardInterrupt: raise
      except: pass
    if miller_arrays is None:
      print >> sys.stderr, "Error: unknown file format:", file_name
      print >> sys.stderr
      return 1
    for mi in miller_arrays:
      if mi.is_xray_intensity_array():
        print >> sys.stderr, "Warning: converting intensities to amplitudes"
        mi = mi.f_sq_as_f()
      elif not mi.is_xray_amplitude_array():
        continue
      solve(mi, command_line)

def solve(miller_array, command_line):
  flipped = charge_flipping.basic_charge_flipping_iterator(
    f_obs=miller_array,
    relative_delta=command_line.options.relative_delta)
  for i,state in enumerate(flipped):
    if i == command_line.options.max_iterations:
      print >> sys.stderr, "Maximum number of iterations reached: terminating"
      break
    r1 = state.r1_factor()
  print "R_1 = %.3f" % r1


if __name__ == '__main__':
  import sys
  import os.path
  sys.exit(run(args=sys.argv[1:]))
