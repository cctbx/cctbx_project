# LIBTBX_SET_DISPATCHER_NAME smtbx.absolute_structure

from cctbx.array_family import flex
from cctbx import miller
from cctbx import xray

import iotbx.builders
import iotbx.cif
import iotbx.cif.builders
from iotbx import shelx
from iotbx.shelx import hklf
from iotbx.option_parser import option_parser

from smtbx import absolute_structure

import glob, os, sys

def crawl(directory, ext='cif'):
  assert ext in ('res', 'ins', 'fcf', 'cif')
  for root, dirs, files in os.walk(directory):
    if '.olex' in root: continue # ignore Olex2 strdir subdirectories
    g = glob.glob(os.path.join(root, "*.%s" %ext))
    for path in g:
      try:
        run_once(path)
      except Exception, e:
        continue

def run_once(file_path):
  file_root, file_ext = os.path.splitext(file_path)
  hkl_path = file_root + '.hkl'
  fcf_path = file_root + '.fcf'
  if file_ext == '.fcf':
    fo2, fc, scale = structure_factors_from_fcf(file_path)
  elif file_ext == '.cif':
    cif = iotbx.cif.reader(file_path=file_path).model()
    cif_block = cif.values()[0]
    wavelength = float(cif_block['_diffrn_radiation_wavelength'])
    xs = iotbx.cif.builders.crystal_structure_builder(cif_block).structure
    xs.set_inelastic_form_factors(photon=wavelength, table="sasaki")
    if os.path.exists(fcf_path):
      fo2, fc, scale = structure_factors_from_fcf(fcf_path, xs)
    elif os.path.exists(hkl_path):
      fo2, fc, scale = structure_factors_from_hkl(hkl_path, xs)
  else:
    if not os.path.exists(hkl_path): return
    fo2, fc, scale = structure_factors_from_ins_res(file_path)
  if not fc.space_group().is_centric():
    print file_path
    absolute_structure_analysis(fo2, fc, scale)

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
  return fo2, fc, scale

def structure_factors_from_hkl(xs, hkl_path):
  fo2 = hklf.reader(filename=hkl_path).as_miller_arrays(
    crystal_symmetry=xs)[0]
  fo2.set_observation_type_xray_intensity()
  merging = fo2.merge_equivalents()
  fo2 = merging.array()
  fc = fo2.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
  scale = fo2.scale_factor(fc)
  return fo2, fc, scale

def structure_factors_from_ins_res(file_path):
  hkl_path = os.path.splitext(file_path)[0] + ".hkl"
  if not os.path.exists(hkl_path): return [None]*3
  builder = iotbx.builders.crystal_structure_builder()
  stream = shelx.command_stream(filename=file_path)
  l_ins = shelx.instruction_parser(stream, builder)
  l_cs = shelx.crystal_symmetry_parser(l_ins.filtered_commands(), builder)
  l_xs = shelx.atom_parser(l_cs.filtered_commands(), builder)
  l_xs.parse()
  xs = builder.structure
  xs.set_inelastic_form_factors(
    photon=l_ins.instructions['cell'][0], table="sasaki")
  return structure_factors_from_hkl(xs, hkl_path)

def absolute_structure_analysis(fo2, fc, scale):
  if not fc.space_group().is_centric():
    hooft_analysis = absolute_structure.hooft_analysis(
      fo2, fc, scale_factor=scale)
    print "Gaussian analysis:"
    hooft_analysis.show()
    NPP = absolute_structure.bijvoet_differences_probability_plot(
      hooft_analysis)
    print "Probability plot:"
    NPP.show()
    print

def run(args):
  command_line = (option_parser(
    usage="smtbx.absolute_structure directory|cif|fcf|ins/res [options]")
                  .enable_symmetry_comprehensive()
                  .option(None, "--ext",
                          action="store")
                  .option(None, "--debug",
                          action="store_true")
                  .option(None, "--verbose",
                          action="store_true")
                  ).process(args=args)
  if len(command_line.args) != 1:
    command_line.parser.show_help()
    return
  if os.path.isdir(command_line.args[0]):
    crawl(command_line.args[0], ext=command_line.options.ext)
  else:
    run_once(command_line.args[0])


if __name__ == '__main__':
  run(sys.argv[1:])
