# LIBTBX_SET_DISPATCHER_NAME smtbx.absolute_structure

from cctbx.array_family import flex

import iotbx.builders
import iotbx.cif
import iotbx.cif.builders
from iotbx import shelx
from iotbx.shelx import hklf
from iotbx.option_parser import option_parser

from smtbx import absolute_structure

import glob, os, sys

def crawl(directory, ext='cif', log=None):
  assert ext in ('res', 'ins', 'fcf', 'cif')
  for root, dirs, files in os.walk(directory):
    if '.olex' in root: continue # ignore Olex2 strdir subdirectories
    g = glob.glob(os.path.join(root, "*.%s" %ext))
    for path in g:
      try:
        run_once(path, log=log)
      except Exception, e:
        continue

def run_once(file_path, nu=None, log=None):
  if log is None:
    log = sys.stdout
  file_root, file_ext = os.path.splitext(file_path)
  hkl_path = file_root + '.hkl'
  fcf_path = file_root + '.fcf'
  if file_ext == '.fcf':
    xs, fo2, fc, scale = structure_factors_from_fcf(file_path)
  elif file_ext == '.cif':
    cif = iotbx.cif.reader(file_path=file_path).model()
    cif_block = cif.values()[0]
    wavelength = float(cif_block['_diffrn_radiation_wavelength'])
    xs = iotbx.cif.builders.crystal_structure_builder(cif_block).structure
    xs.set_inelastic_form_factors(photon=wavelength, table="sasaki")
    if os.path.exists(fcf_path):
      xs, fo2, fc, scale = structure_factors_from_fcf(fcf_path, xs)
    elif os.path.exists(hkl_path):
      xs, fo2, fc, scale = structure_factors_from_hkl(hkl_path, xs)
    else: return
  else:
    if not os.path.exists(hkl_path): return
    xs, fo2, fc, scale = structure_factors_from_ins_res(file_path)
  if not fc.space_group().is_centric():
    print >> log, file_path
    absolute_structure_analysis(xs, fo2, fc, scale, nu=nu, log=log)
  log.flush()

def structure_factors_from_fcf(file_path, xs=None):
  cif = iotbx.cif.reader(file_path=file_path).model()
  cif_block = cif.values()[0]
  if '_shelx_refln_list_code' in cif_block:
    assert cif_block['_shelx_refln_list_code'] == '4'
  arrays = iotbx.cif.builders.miller_array_builder(cif_block).arrays()
  fo2 = arrays['_refln_F_squared_meas']
  fc2 = arrays['_refln_F_squared_calc']
  if xs is None:
    fc = fc2.f_sq_as_f().phase_transfer(flex.double(fc2.size(), 0))
    scale = 1
  else:
    fc = fo2.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
    scale = fo2.scale_factor(fc)
  return xs, fo2, fc, scale

def structure_factors_from_hkl(xs, hkl_path, weighting_scheme=None):
  fo2 = hklf.reader(filename=hkl_path).as_miller_arrays(
    crystal_symmetry=xs)[0]
  fo2.set_observation_type_xray_intensity()
  merging = fo2.merge_equivalents()
  fo2 = merging.array()
  xs.scattering_type_registry(table="it1992", d_min=fo2.d_min())
  fc = fo2.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
  scale = fo2.scale_factor(fc)
  if weighting_scheme is not None:
    weights = weighting_scheme(
      fo2.data(), fo2.sigmas(), fc.as_intensity_array().data(), scale)
    scale = fo2.scale_factor(fc, weights=weights)
  return xs, fo2, fc, scale

def structure_factors_from_ins_res(file_path):
  from iotbx.builders \
       import weighted_constrained_restrained_crystal_structure_builder
  hkl_path = os.path.splitext(file_path)[0] + ".hkl"
  if not os.path.exists(hkl_path): return [None]*3
  builder = weighted_constrained_restrained_crystal_structure_builder()
  stream = iotbx.shelx.command_stream(filename=file_path)
  l_ins = iotbx.shelx.instruction_parser(stream, builder)
  stream = iotbx.shelx.crystal_symmetry_parser(l_ins.filtered_commands(),
                                               builder)
  stream = iotbx.shelx.wavelength_parser(stream.filtered_commands(), builder)
  stream = iotbx.shelx.afix_parser(stream.filtered_commands(), builder)
  stream = iotbx.shelx.atom_parser(stream.filtered_commands(), builder)
  stream = iotbx.shelx.restraint_parser(stream.filtered_commands(), builder)
  stream.parse()
  xs = builder.structure
  twin = l_ins.instructions.get('twin')
  if twin is not None and not xs.space_group().is_centric():
    print file_path
    print 'twin: ', twin['matrix']
  xs.set_inelastic_form_factors(
    photon=builder.wavelength_in_angstrom, table="sasaki")
  return structure_factors_from_hkl(
    xs, hkl_path, weighting_scheme=builder.weighting_scheme)

def absolute_structure_analysis(xs, fo2, fc, scale, nu=None, log=None):
  if log is None:
    log = sys.stdout
  if not fc.space_group().is_centric():
    hooft_analysis = absolute_structure.hooft_analysis(
      fo2, fc, scale_factor=scale)
    print >> log, "Gaussian analysis:"
    hooft_analysis.show(out=log)
    NPP = absolute_structure.bijvoet_differences_probability_plot(
      hooft_analysis)
    print >> log, "Probability plot:"
    NPP.show(out=log)
    print >> log
    if nu is None:
      nu = absolute_structure.maximise_students_t_correlation_coefficient(
        NPP.y, min_nu=1, max_nu=200)
    t_analysis = absolute_structure.students_t_hooft_analysis(
      fo2, fc, nu, scale_factor=scale, probability_plot_slope=NPP.fit.slope())
    tPP = absolute_structure.bijvoet_differences_probability_plot(
      t_analysis, use_students_t_distribution=True, students_t_nu=nu)
    print >> log, "Student's t analysis:"
    print >> log, "nu: %.2f" %nu
    t_analysis.show(out=log)
    print >> log, "Probability plot:"
    tPP.show(out=log)
    print >> log
    if xs is not None:
      flack = absolute_structure.flack_analysis(xs, fo2)
      flack.show(out=log)

def run(args):
  command_line = (option_parser(
    usage="smtbx.absolute_structure directory|cif|fcf|ins/res [options]")
                  .enable_symmetry_comprehensive()
                  .option(None, "--ext",
                          action="store")
                  .option(None, "--nu",
                          action="store",
                          type="float")
                  .option(None, "--debug",
                          action="store_true")
                  .option(None, "--verbose",
                          action="store_true")
                  .option(None, "--log",
                          action="store")
                  ).process(args=args)
  if len(command_line.args) != 1:
    command_line.parser.show_help()
    return
  if command_line.options.log is not None:
    log = open(command_line.options.log, 'wb')
  else:
    log = None
  if os.path.isdir(command_line.args[0]):
    crawl(command_line.args[0], ext=command_line.options.ext, log=log)
  elif os.path.isfile(command_line.args[0]):
    run_once(command_line.args[0], nu=command_line.options.nu, log=log)
  else:
    print "Please provide a valid file or directory"


if __name__ == '__main__':
  run(sys.argv[1:])
