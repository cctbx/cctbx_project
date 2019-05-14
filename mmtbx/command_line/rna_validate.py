# LIBTBX_SET_DISPATCHER_NAME phenix.rna_validate
# LIBTBX_SET_DISPATCHER_NAME molprobity.rna_validate

from __future__ import absolute_import, division, print_function
from mmtbx.validation.rna_validate import rna_validation
import iotbx.phil
import os, sys

def get_master_phil():
  return iotbx.phil.parse(input_string="""
  model = None
      .type = path
      .help = "Model file (PDB or mmCIF)"

  outliers_only = True
    .type = bool
    .help = "Only display outliers"

  verbose = True
    .type = bool
    .help = '''Verbose'''
  rna_sugar_pucker_analysis
    .short_caption = RNA sugar pucker analysis
    .style = box noauto auto_align menu_item parent_submenu:advanced
  {
    include scope mmtbx.monomer_library.rna_sugar_pucker_analysis.master_phil
  }
""", process_includes=True)
prog = os.getenv('LIBTBX_DISPATCHER_NAME')
usage_string = """
%(prog)s file.pdb [params.eff] [options ...]

Options:

  pdb=input_file        input PDB file
  outliers_only=False   print all suites, not just outliers

Example:

  %(prog)s pdb=1u8d.pdb outliers_only=False

""" % locals()

def run(args, out=sys.stdout, quiet=False):
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx.monomer_library import server
  import iotbx.pdb
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()
  if (params.model is None):
    raise Usage(usage_string)
  pdb_in = iotbx.pdb.input(source_info=params.model, file_name=params.model)
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    pdb_inp=pdb_in,
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  geometry = processed_pdb_file.geometry_restraints_manager()
  result = rna_validation(
    pdb_hierarchy=pdb_hierarchy,
    geometry_restraints_manager=geometry,
    params=params,
    outliers_only=params.outliers_only)
  result.show(out=out)
  return result

if (__name__ == "__main__"):
  run(sys.argv[1:])
