# LIBTBX_SET_DISPATCHER_NAME smtbx.refine

import os
from scitbx import lstbx
import scitbx.lstbx.normal_eqns_solving
from smtbx import refinement

allowed_input_file_extensions = ('.ins', '.res', '.cif')

class command_line_error(RuntimeError): pass

class number_of_arguments_error(command_line_error): pass

def run(filenames, options):
  #interpret file names
  if not 1 <= len(filenames) <= 3:
    raise number_of_arguments_error()
  input_filename = filenames[0]
  output_filename = None
  reflections_filename = None
  if len(filenames) == 3:
    reflections_filename, output_filename = filenames[1:]
  elif len(filenames) == 2:
    _, ext = os.path.splitext(filenames[1])
    if ext in allowed_input_file_extensions:
      output_filename = filenames[1]
    else:
      reflections_filename = filenames[1]

  # adjust file names
  in_root, in_ext = os.path.splitext(input_filename)
  if not in_ext: in_ext = '.ins'
  if reflections_filename is None:
    reflections_filename = in_root + '.hkl'
  if output_filename is None:
    if in_ext == '.ins':
      out_root, out_ext = in_root, '.res'
    else:
      out_root, out_ext = in_root + '-out', in_ext
    output_filename = out_root + out_ext
  else:
    out_root, out_ext = os.path.splitext(output_filename)

  # check extensions are supported
  for ext in (in_ext, out_ext):
    if ext not in allowed_input_file_extensions:
      raise command_line_error("unsupported extension: %s" % ext)

  # Investigate whether input and ouput files do exist, are the same, etc
  for filename in (input_filename, reflections_filename):
    if not os.path.isfile(filename):
      raise command_line_error("No such file %s" % filename)
  if os.path.isfile(output_filename):
    if options.overwrite:
      print >>sys.stderr, ("\nWARNING: file %s will be overwritten.\n"
                           % output_filename)
    else:
      raise command_line_error(
        "refuse to overwrite file %s (use option 'overwrite' to force it)"
        % output_filename)

  # Load input model and reflections
  if in_ext == '.cif':
    xm = refinement.model.from_cif(model=input_filename,
                                   reflections=reflections_filename)
  else:
    xm = refinement.model.from_shelx(ins_or_res=input_filename,
                                     hkl=reflections_filename)

  # At last...
  for sc in xm.xray_structure.scatterers():
    sc.flags.set_grad_site(True)
    if sc.flags.use_u_iso(): sc.flags.set_grad_u_iso(True)
    if sc.flags.use_u_aniso(): sc.flags.set_grad_u_aniso(True)
  steps = lstbx.normal_eqns_solving.naive_iterations(
    non_linear_ls=xm.least_squares(),
    n_max_iterations=options.max_cycles,
    gradient_threshold=options.stop_if_max_derivative_below,
    step_threshold=options.stop_if_shift_norm_below)
  steps.do()

  # Write result to disk
  if out_ext != '.cif':
    raise NotImplementedError("Write refined structure to %s file" % out_ext)
  with open(output_filename, 'w') as out:
    xm.xray_structure.as_cif_simple(out)

here_usage="""\
refine [options] INPUT
refine [options] INPUT REFLECTIONS
refine [options] INPUT OUTPUT
refine [options] INPUT REFLECTIONS OUTPUT
"""

here_description = """\
Supported formats for the INPUT and OUTPUT files are CIF and ShelX
whereas REFLECTIONS file format may be any of those supported
by iotbx.reflection_file_reader.

A certain amount of guessing is done by the program:
- if INPUT lacks any extension, a .ins extension is appended;
- if OUTPUT is not specified,
  o if INPUT is a .ins file, then OUTPUT is a .res file with the same root;
  o otherwise OUTPUT has the same extension as INPUT but with "-out"
    appended to INPUT root.
- if REFLECTIONS is not specified, it is INPUT root with a .hkl extension
"""

if __name__ == '__main__':
  import sys, optparse
  parser = optparse.OptionParser(
    usage=here_usage,
    description=here_description)
  parser.add_option(
    '--overwrite',
    action='store_true',
    help='allows OUTPUT to be overwritten')
  parser.add_option(
    '--stop-if-max-derivative-below',
    type='float',
    default=1e-7,
    help='Stop refinement as soon as the largest absolute value of objective '
         'derivatives is below the given threshold.')
  parser.add_option(
    '--stop-if-shift-norm-below',
    type='float',
    default=1e-7,
    help='Stop refinement as soon as the Euclidean norm of the vector '
         'of parameter shifts is below the given threshold.')
  parser.add_option(
    '--max-cycles',
    type='int',
    default=8,
    help='Stop refinement as soon as the given number of cycles have been '
         'performed')
  options, args = parser.parse_args()
  try:
    run(args, options)
  except number_of_arguments_error:
    parser.print_usage()
    sys.exit(1)
  except command_line_error, err:
    print >>sys.stderr, "\nERROR: %s\n" % err
    parser.print_help()
    sys.exit(1)
