# LIBTBX_SET_DISPATCHER_NAME phenix.map_value_at_point

import sys
import mmtbx.utils
from libtbx.utils import Sorry
import iotbx.phil
import iotbx.pdb
from iotbx import reflection_file_reader
from iotbx.pdb import combine_unique_pdb_files
from scitbx.array_family import flex

msg="""\

phenix.map_value_at_point: tool to compute map value at a given point using
eight-point interpolation.

Inputs:
  - MTZ file containing Fourier map coefficients;
  - label selecting which map coefficients should be used (in case there are
    multiple choices in input file, there is no need to provide label otherwise);
  - the coordinates of points at which the map value will be computed can be
    defined in two ways: either providing a triplet of xyz coordinates using
    'point' keyword, or giving a PDB file in the command line in which case the
    atomic coordinates will serve as points;
  - optionally the 'resolution_factor' can be specified. It defines the map
    calculation grid step as grid_step=resolution_factor*d_min
  - optionally, type of map scaling can be defined: sigma or volume scaled
    (sigma is the default).

Usage examples:
  1. phenix.map_value_at_point map_coeffs.mtz point="0.5 1.2 3" point="0.53 3.3 2.5" label="2mFo-DFc"
  2. phenix.map_value_at_point map_coeffs.mtz model.pdb
  3. phenix.map_value_at_point map_coeffs.mtz model.pdb point="1 2 3" resolution_factor=0.3

  IMPORTANT: point must be a triplet of Cartesian coordinates (not fractional)!

Output:
  list of map values at specified points.
"""

master_params_str="""\
label = None
  .type = str
point = None
  .type = floats(3)
  .multiple = True
resolution_factor = 0.25
  .type = float
scale = *sigma volume
  .type = choice(multi=False)
"""

def defaults(log):
  print >> log, "Default params::\n"
  parsed = iotbx.phil.parse(master_params_str)
  parsed.show(prefix="  ", out=log)
  print >> log
  return parsed

def run(args, log = sys.stdout):
  if(len(args)==0):
    print >> log, msg
    defaults(log=log)
    return
  #
  parsed = defaults(log=log)
  processed_args = mmtbx.utils.process_command_line_args(args = args,
    log = sys.stdout, master_params = parsed)
  params = processed_args.params.extract()
  reflection_files = processed_args.reflection_files
  #
  atoms_with_labels = None
  if(len(processed_args.pdb_file_names)==1):
    pdb_combined = combine_unique_pdb_files(
      file_names=processed_args.pdb_file_names)
    pdb_combined.report_non_unique()
    pdb_inp = iotbx.pdb.input(source_info = None,
      lines = flex.std_string(pdb_combined.raw_records))
    atoms_with_labels = pdb_inp.atoms_with_labels()
  #
  if(len(reflection_files) == 0):
    raise Sorry("No reflection file found.")
  if(len(reflection_files) > 1):
    raise Sorry("More than one reflection file found.")
  crystal_symmetry = processed_args.crystal_symmetry
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  reflection_file_name = reflection_files[0].file_name()
  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    reflection_file_name).as_miller_arrays()
  #
  if(len(miller_arrays)==1 and params.label is None):
    ma = miller_arrays[0]
  elif(len(miller_arrays)>1 and params.label is None):
    raise Sorry("Multiple data columns found in input file. Use label keyword to select the one.")
  elif(len(miller_arrays)==1 and params.label is not None):
    ma = miller_arrays[0]
    if(ma.info().labels[0].lower() != params.label.lower()):
      raise Sorry("Specified label 'label=%s' does not match any label in input file."%params.label)
  elif(len(miller_arrays)>1 and params.label is not None):
    for ma in miller_arrays:
      if(ma.info().labels[0].lower() == params.label.lower()):
        break
  if(not ma.is_complex_array()):
    raise Sorry("Data must be complex type (real provided).")
  print "Input reflection data (Fourier map coefficients):"
  ma.show_comprehensive_summary(prefix="  ")
  print
  #
  if(len(params.point)==0 and atoms_with_labels is None):
    raise Sorry("No points given to compute map value at.")
  else:
    fft_map = ma.fft_map(resolution_factor=params.resolution_factor)
    if(params.scale == "sigma"):
      fft_map.apply_sigma_scaling()
      print "Using sigma scaled map.\n"
    else:
      fft_map.apply_volume_scaling()
      print "Using volume scale map.\n"
    map_3d = fft_map.real_map_unpadded()
    print "Map values at specified points:"
    for point in params.point:
      point_frac = ma.unit_cell().fractionalize(point)
      point_formatted = ",".join([str("%10.3f"%p).strip() for p in point])
      point_frac_formatted = \
        ",".join([str("%10.3f"%p).strip() for p in point_frac])
      map_value = str(
        "%10.3f"%map_3d.eight_point_interpolation(point_frac)).strip()
      print "  Input point: (%s) Fractionalized: (%s) Map value: %s"%(
        point_formatted, point_frac_formatted, map_value)
    #
    if(atoms_with_labels is not None):
      for point in atoms_with_labels:
        point_frac = ma.unit_cell().fractionalize(point.xyz)
        point_formatted = ",".join([str("%8.3f"%p) for p in point.xyz])
        map_value = str(
          "%10.3f"%map_3d.eight_point_interpolation(point_frac)).strip()
        print point.quote(), "Point: %s Map value: %s"%(point_formatted,map_value)
  #
  print
  print "All done."

if(__name__ == "__main__"):
  run(sys.argv[1:])
