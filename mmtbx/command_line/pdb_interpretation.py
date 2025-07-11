"""Read PDB file and build restraints for refinement (useful for trouble-shooting)"""
# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_interpretation
# LIBTBX_SET_DISPATCHER_NAME mmtbx.pdb_interpretation
from __future__ import absolute_import, division, print_function
from iotbx.file_reader import any_file
import sys

def get_master_phil():
  import iotbx.phil
  return iotbx.phil.parse("""\
strict_processing = False
  .type = bool
build_geometry_restraints_manager = True
  .type = bool
build_xray_structure = True
  .type = bool
max_atoms = None
  .type = int
write_geo_files = False
  .type = bool
write_tardy_geo_files = False
  .type = bool
model_file_name = None
  .type = path
  .multiple = True
restraints_cif_file_name = None
  .type = path
  .multiple = True
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
""", process_includes=True)

def run(args):
  from libtbx.str_utils import show_string
  #
  master_phil = get_master_phil()
  import iotbx.phil
  input_objects = iotbx.phil.process_command_line_with_files(
    args=args,
    pdb_file_def="model_file_name",
    cif_file_def="restraints_cif_file_name",
    master_phil=master_phil)
  input_objects.work.show()
  work_params = input_objects.work.extract()
  #
  import mmtbx.monomer_library.server
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  for file_name in work_params.restraints_cif_file_name:
    print("Processing CIF file: %s" % show_string(file_name))
    af = any_file(file_name = file_name)
    for i, srv in enumerate([mon_lib_srv, ener_lib]):
      srv.process_cif_object(
        cif_object=af.file_object.model(),
        file_name=af.file_name,
        process_tor=not i)
  #
  import mmtbx.monomer_library.pdb_interpretation
  from libtbx.utils import Sorry
  processed_pdb_files = []
  for file_name in work_params.model_file_name:
    processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=work_params.pdb_interpretation,
      file_name=file_name,
      strict_conflict_handling=work_params.strict_processing,
      substitute_non_crystallographic_unit_cell_if_necessary=True,
      max_atoms=work_params.max_atoms,
      log=sys.stdout)
    if (work_params.strict_processing):
      msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
      if (msg is not None):
        raise Sorry(msg)
    if (work_params.build_geometry_restraints_manager):
      processed_pdb_file.geometry_restraints_manager(
        params_edits=work_params.geometry_restraints.edits,
        params_remove=work_params.geometry_restraints.remove,
        )
    if (work_params.build_xray_structure):
      processed_pdb_file.xray_structure()
    processed_pdb_files.append(processed_pdb_file)
  print()
  if (   work_params.write_geo_files
      or work_params.write_tardy_geo_files):
    import os.path as op
    for processed_pdb_file in processed_pdb_files:
      acp = processed_pdb_file.all_chain_proxies
      source = acp.pdb_inp.source_info()
      assert source.startswith("file ")
      pdb_file_name = source[5:]
      sites_cart = acp.sites_cart_exact()
      site_labels = [atom.id_str() for atom in acp.pdb_atoms]
      def write_geo(label, geo, geo_file_name):
        from libtbx.utils import date_and_time
        header = "# %sgeometry restraints for file:\n" % label
        header += "#   %s\n# %s\n" % (show_string(pdb_file_name),
            date_and_time())
        geo.write_geo_file(
            sites_cart=sites_cart,
            site_labels=site_labels,
            file_name=geo_file_name,
            header=header)
      geo = processed_pdb_file.geometry_restraints_manager()
      if (work_params.write_geo_files):
        geo_file_name = op.basename(pdb_file_name) + ".geo"
        print("Writing file: %s" % show_string(geo_file_name))
        write_geo(label="", geo=geo, geo_file_name=geo_file_name)
        print()
      if (work_params.write_tardy_geo_files):
        geo_file_name = op.basename(pdb_file_name) + ".tardy.geo"
        print("Writing file: %s" % show_string(geo_file_name))
        tardy_tree = geo.construct_tardy_tree(sites_cart=sites_cart)
        reduced_geo = geo.reduce_for_tardy(tardy_tree=tardy_tree)
        write_geo(label="tardy ", geo=reduced_geo, geo_file_name=geo_file_name)
        print()
  return processed_pdb_files

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

