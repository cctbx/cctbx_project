from iotbx.kriber import strudat
import iotbx.pdb
import iotbx.cif
from iotbx.option_parser import option_parser
from cctbx.crystal import coordination_sequences
from cctbx import crystal, xray
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
import sys

def display(
      distance_cutoff,
      show_cartesian,
      max_shell,
      coseq_dict,
      xray_structure):
  xray_structure.show_summary().show_scatterers()
  print
  pairs = xray_structure.show_distances(
    distance_cutoff=distance_cutoff,
    show_cartesian=show_cartesian,
    keep_pair_asu_table=True)
  print
  if (pairs.pair_counts.size() <= 15):
    print "Pair counts:", list(pairs.pair_counts)
    print
  if (max_shell is None):
    term_table = None
  else:
    term_table = crystal.coordination_sequences.simple(
      pair_asu_table=pairs.pair_asu_table,
      max_shell=max_shell)
    coordination_sequences.show_terms(
      structure=xray_structure,
      term_table=term_table,
      coseq_dict=coseq_dict)
    print
  return pairs, term_table

def run(args):
  command_line = (option_parser(
    usage="iotbx.show_distances [options] studat_file [...]",
    description="Example: iotbx.show_distances strudat --tag=SOD")
    .option(None, "--cif_data_block_name",
      action="store",
      type="string",
      default=None,
      help="data block name as it appears in the CIF file")
    .option(None, "--tag",
      action="store",
      type="string",
      help="tag as it appears in the strudat file")
    .option(None, "--distance_cutoff",
      action="store",
      type="float",
      default=5,
      help="Maximum distance to be considered",
      metavar="FLOAT")
    .option(None, "--min_distance_sym_equiv",
      action="store",
      type="float",
      default=0.5,
      help="Minimum distance between symmetry mates"
           " (for special position analysis)",
      metavar="FLOAT")
    .option(None, "--show_cartesian",
      action="store_true",
      help="Show Cartesian coordinates (instead of fractional)")
    .enable_symmetry_comprehensive()
    .option(None, "--cs",
      action="store",
      type="int",
      help="Compute N terms of the coordination sequences",
      metavar="N")
    .option(None, "--coseq",
      action="store",
      type="string",
      help="name of file with known coordination sequences",
      metavar="FILE")
  ).process(args=args)
  if (len(command_line.args) == 0):
    command_line.parser.show_help()
    return
  co = command_line.options
  max_shell = co.cs
  if (co.coseq is not None):
    coseq_dict = coordination_sequences.get_kriber_coseq_file(
      file_name=co.coseq)
    if (max_shell is None): max_shell = 10
  else:
    coseq_dict = None
  def call_display(xray_structure):
    display(
      distance_cutoff=co.distance_cutoff,
      show_cartesian=co.show_cartesian,
      max_shell=max_shell,
      coseq_dict=coseq_dict,
      xray_structure=xray_structure)
  cif_data_block_name_use_counter = 0
  for file_name in command_line.args:
    xray_structure = None
    if (iotbx.pdb.is_pdb_file(file_name=file_name)):
      xray_structure = iotbx.pdb.input(
        file_name=file_name).xray_structure_simple(
          crystal_symmetry=command_line.symmetry,
          cryst1_substitution_buffer_layer=max(5,
            co.distance_cutoff+1),
          min_distance_sym_equiv=co.min_distance_sym_equiv,
          enable_scattering_type_unknown=True)
      call_display(xray_structure)
      continue
    if (file_name.lower().endswith(".cif")):
      xray_structures = iotbx.cif.reader(
        file_path=file_name).build_crystal_structures()
      if (co.cif_data_block_name is not None):
        xray_structure = xray_structures.get(co.cif_data_block_name)
        if (xray_structure is None):
          continue
        cif_data_block_name_use_counter += 1
        xray_structures = [xray_structure]
      else:
        xray_structures = xray_structures.values()
      for xray_structure in xray_structures:
        call_display(xray_structure)
      continue
    if (   file_name.endswith(".ins")
        or file_name.endswith(".res")):
      xray_structure = xray.structure.from_shelx(
        filename=file_name, strictly_shelxl=False)
      call_display(xray_structure)
      continue
    if (command_line.symmetry is not None
        and (command_line.symmetry.unit_cell() is not None
          or command_line.symmetry.space_group_info() is not None)):
      raise Sorry(
        "Command-line symmetry options not supported for strudat files.")
    strudat_entries = strudat.read_all_entries(open(file_name))
    for entry in strudat_entries.entries:
      if (    co.tag is not None
          and co.tag != entry.tag):
        continue
      print "strudat tag:", entry.tag
      print
      xray_structure = entry.as_xray_structure(
        min_distance_sym_equiv=co.min_distance_sym_equiv)
      call_display(xray_structure)
  if (    co.cif_data_block_name is not None
      and cif_data_block_name_use_counter == 0):
    raise Sorry(
      "cif_data_block_name %s not found in any input files"
        % show_string(co.cif_data_block_name))

if (__name__ == "__main__"):
  run(sys.argv[1:])
