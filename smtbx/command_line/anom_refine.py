from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME smtbx.anom_refine

import os, sys, argparse
from scitbx import lstbx
import scitbx.lstbx.normal_eqns_solving
from smtbx import refinement
from timeit import default_timer as current_time
from math import sqrt


allowed_input_file_extensions = ('.ins', '.res', '.cif')

class command_line_error(RuntimeError): pass

class number_of_arguments_error(command_line_error): pass

class energy_missing_error(RuntimeError): pass

def make_parser():
  parser = argparse.ArgumentParser(
      description='''Refinement of inelastic scattering factors. Contact: Daniel
Paley (dwp2111@columbia.edu)''')
  parser.add_argument(
      'ref_structure',
      type=str,
      help='''A structure in .cif or .ins format determined far from any
absorption edge. Coordinates, displacements, and occupancies will be fixed
at their values in this file.'''
      )
  parser.add_argument(
      'reflections',
      type=str,
      help='''A .hkl file containing intensities (ShelX HKLF 4 format).'''
      )
  parser.add_argument(
      'anom_atom',
      type=str,
      help='''The atom type for which f' and f" will be refined. All atoms of
this type will be refined independently.'''
      )
  group = parser.add_mutually_exclusive_group()
  group.add_argument(
    '-e', '--energy',
    type=float,
    default=None,
    help='''Beam energy in eV.'''
    )
  group.add_argument(
    '-E', '--energy-in-fname',
    type=str,
    default=None,
    help='''First digit and length of energy if given in hkl filename. Format
start,length.'''
    )
  parser.add_argument(
    '-t', '--table',
    action='store_true',
    help='Output in condensed format suitable for further processing')
  parser.add_argument(
    '-T', '--table-with-su',
    action='store_true',
    help='Output in condensed format with SUs under the values')
  parser.add_argument(
    '-c', '--max-cycles',
    type=int,
    default=100,
    help='Stop refinement as soon as the given number of cycles have been '
         'performed')
  parser.add_argument(
    '-o', '--outfile',
    type=str,
    default=None,
    help='Write output to filename OUTFILE (default: print to stdout)')
  parser.add_argument(
    '-O', '--overwrite',
    action='store_true',
    help='No error if OUTFILE exists')
  parser.add_argument(
    '-d', '--stop-deriv',
    type=float,
    default=1e-7,
    help='Stop refinement as soon as the largest absolute value of objective '
         'derivatives is below the given threshold.')
  parser.add_argument(
    '-s', '--stop-shift',
    type=float,
    default=1e-7,
    help='Stop refinement as soon as the Euclidean norm of the vector '
         'of parameter shifts is below the given threshold.')
  parser.add_argument(
    '-A', '--adp_global',
    action='store_true',
    help='Refine a global ADP scaling parameter (an extremely simple model for'
         'beam damage)')
  return parser

def run(args):

  # adjust file names
  in_root, in_ext = os.path.splitext(args.ref_structure)

  # Check that input files exist
  for filename in (args.ref_structure, args.reflections):
    if not os.path.isfile(filename):
      raise command_line_error("No such file %s" % filename)

  # Check output file
  if args.outfile and os.path.isfile(args.outfile) and not args.overwrite:
    raise command_line_error("Output file {} exists.".format(args.outfile))

  # Load input model and reflections
  if in_ext == '.cif':
    xm = refinement.model.from_cif(model=args.ref_structure,
                                   reflections=args.reflections + '=hklf4')
  else:
    xm = refinement.model.from_shelx(ins_or_res=args.ref_structure,
                                     hkl=args.reflections,
                                     strictly_shelxl=False)

  # Look for beam energy
  if args.energy:
    energy = args.energy
    wvl = 12398 / energy
  elif args.energy_in_fname:
    estart,elength = args.energy_in_fname.split(',')
    estart = int(estart)
    elength = int(elength)
    energy = float(args.reflections[estart:estart+elength])
    wvl = 12398 / energy
  else:
    energy = None
    sys.stderr.write('''WARNING: Using beam energy from reference model. \
Inelastic form factors for \n non-refined atoms may be inaccurate.\n''')
    wvl = xm.wavelength




  # Load default anomalous scattering factors if wavelength is available
  if wvl:
    xm.xray_structure.set_inelastic_form_factors(wvl, 'sasaki')
  else:
    raise energy_missing_error()

  # At last...
  anom_sc_list=[]
  for i, sc in enumerate(xm.xray_structure.scatterers()):
    sc.flags.set_grad_site(False)
    sc.flags.set_grad_u_iso(False)
    sc.flags.set_grad_u_aniso(False)
    sc.flags.set_grad_occupancy(False)
    if sc.element_symbol().lower() == args.anom_atom.lower():
      sc.flags.set_use_fp_fdp(True)
      sc.flags.set_grad_fp(True)
      sc.flags.set_grad_fdp(True)
      anom_sc_list.append((i, sc))

  if args.adp_global:
    from smtbx.refinement.constraints.adp import scalar_scaled_u
    for sc in xm.xray_structure.scatterers():
      if sc.flags.use_u_aniso():
        sc.flags.set_grad_u_aniso(True)
      else:
        sc.flags.set_grad_u_iso(True)
    adp_scale = scalar_scaled_u(range(len(xm.xray_structure.scatterers())))
    xm.constraints.append(adp_scale)


  ls = xm.least_squares()
  steps = lstbx.normal_eqns_solving.levenberg_marquardt_iterations(
    non_linear_ls=ls,
    n_max_iterations=args.max_cycles,
    gradient_threshold=args.stop_deriv,
    step_threshold=args.stop_shift)
  cov_an = ls.covariance_matrix_and_annotations()

  # Prepare output
  result = ''

  if args.table or args.table_with_su:
    if energy: label = "{:.1f} ".format(energy)
    else: label = "{} ".format(args.reflections)
    result += label
    for i, sc in anom_sc_list:
      result += "{:6.3f} ".format(sc.fp)
    for i, sc in anom_sc_list:
      result += "{:6.3f} ".format(sc.fdp)
    result += '\n'

  if args.table_with_su:
    result += ' ' * len(label)
    for i, sc in anom_sc_list:
      label_fp = '{}.fp'.format(sc.label)
      sigma_fp = sqrt(cov_an.variance_of(label_fp))
      result += "{:6.3f} ".format(sigma_fp)
    for i, sc in anom_sc_list:
      label_fdp = '{}.fdp'.format(sc.label)
      sigma_fdp = sqrt(cov_an.variance_of(label_fdp))
      result += "{:6.3f} ".format(sigma_fdp)
    result += '\n'

  if not (args.table or args.table_with_su):
    from libtbx.utils import format_float_with_standard_uncertainty \
        as format_float_with_su
    result += "\n### REFINE ANOMALOUS SCATTERING FACTORS ###\n"
    result += "Reflections: {}\n\n".format(args.reflections)

    if args.adp_global:
      sigma_adp_global = adp_scale.esd(ls)
      result += "ADP scale: {}\n\n".format(
          format_float_with_su(adp_scale.scalar.value, sigma_adp_global))

    for i, sc in anom_sc_list:
      label_fp = '{}.fp'.format(sc.label)
      label_fdp = '{}.fdp'.format(sc.label)
      sigma_fp = sqrt(cov_an.variance_of(label_fp))
      sigma_fdp = sqrt(cov_an.variance_of(label_fdp))
      result += "{}:\n\tfp: {}\n\tfdp: {}\n".format(
          sc.label,
          format_float_with_su(sc.fp, sigma_fp),
          format_float_with_su(sc.fdp, sigma_fdp))

  # Write to file or stdout
  if args.outfile:
    with open(args.outfile, 'w') as f:
      f.write(result)
  else:
    print(result)



if __name__ == '__main__':
  from timeit import default_timer as current_time
  t0 = current_time()
  parser = make_parser()
  args = parser.parse_args()
  try:
    run(args)
  except number_of_arguments_error:
    parser.print_usage()
    sys.exit(1)
  except command_line_error as err:
    print("\nERROR: %s\n" % err, file=sys.stderr)
    parser.print_help()
    sys.exit(1)
  except energy_missing_error as err:
    print('Must provide beam energy on the command line or in the reference '
        'structure.')
    sys.exit(1)
  t1 = current_time()
  if not (args.table or args.table_with_su):
    print("Total time: %.3f s" % (t1 - t0))
