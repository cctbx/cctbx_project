from iotbx import pdb
from iotbx.option_parser import option_parser
from cctbx.array_family import flex
from libtbx.str_utils import show_string
from libtbx import easy_pickle
import sys, os

def run(args, command_name="iotbx.pdb.as_xray_structure"):
  command_line = (option_parser(
    usage=command_name+" [options] pdb_file ...",
    description="Example: %s pdb1ab1.ent" % command_name)
    .enable_symmetry_comprehensive()
    .option(None, "--weak_symmetry",
      action="store_true",
      default=False,
      help="symmetry on command line is weaker than symmetry found in files")
    .option(None, "--ignore_occ_for_site_symmetry",
      action="store_true",
      default=False,
      help="disables non_unit_occupancy_implies_min_distance_sym_equiv_zero")
    .option("-v", "--verbose",
      action="store_true",
      default=False,
      help="show scatterers")
    .option(None, "--pickle",
      action="store",
      type="string",
      help="write all data to FILE ('--pickle .' copies name of input file)",
      metavar="FILE")
    .option(None, "--fake_f_obs_and_r_free_flags_d_min",
      action="store",
      type="float",
      help="write F-calc as F-obs, add random R-free flags (MTZ format)",
      metavar="FLOAT")
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
  co = command_line.options
  d_min = co.fake_f_obs_and_r_free_flags_d_min
  all_structures = []
  for file_name in command_line.args:
    print "file_name:", file_name
    sys.stdout.flush()
    pdb_inp = pdb.input(file_name=file_name)
    structure = pdb_inp.xray_structure_simple(
      crystal_symmetry=command_line.symmetry,
      weak_symmetry=co.weak_symmetry,
      non_unit_occupancy_implies_min_distance_sym_equiv_zero=
        not co.ignore_occ_for_site_symmetry)
    structure.show_summary()
    if (structure.special_position_indices().size() != 0):
      structure.show_special_position_shifts(
        sites_cart_original=pdb_inp.atoms().extract_xyz())
    structure.scattering_type_registry().show(show_gaussians=False)
    if (co.verbose):
      structure.show_scatterers()
    if (d_min is not None and d_min > 0):
      f_obs = abs(structure.structure_factors(
        d_min=d_min, anomalous_flag=False).f_calc())
      f_obs = f_obs.customized_copy(sigmas=flex.sqrt(f_obs.data()))
      r_free_flags = f_obs.generate_r_free_flags(fraction=0.05, max_free=None)
      mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
      mtz_dataset.add_miller_array(
        miller_array=r_free_flags,
        column_root_label="R-free-flags")
      mtz_object = mtz_dataset.mtz_object()
      history = "%s %s" % (command_name, show_string(file_name))
      lines = flex.std_string(["Fake F-obs, R-free-flags"])
      while (len(history) != 0):
        lines.append(history[:77])
        history = history[77:]
      mtz_object.add_history(lines=lines)
      mtz_object.show_summary()
      mtz_file_name = os.path.basename(file_name).replace(".","_") \
                    + "_fake.mtz"
      print "Writing file:", mtz_file_name
      mtz_object.write(file_name=mtz_file_name)
    all_structures.append(structure)
    print
  pickle_file_name = co.pickle
  if (pickle_file_name is not None and len(all_structures) > 0):
    if (pickle_file_name == "."):
      if (len(command_line.args) > 1):
        raise Sorry(
          "Ambiguous name for pickle file (more than one input file).")
      pickle_file_name = os.path.basename(command_line.args[0])
    if (not pickle_file_name.lower().endswith(".pickle")):
      pickle_file_name += ".pickle"
    if (len(all_structures) == 1):
      all_structures = all_structures[0]
    else:
      print
    print "Writing all xray structures to file:", pickle_file_name
    easy_pickle.dump(pickle_file_name, all_structures)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
