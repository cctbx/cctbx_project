# LIBTBX_SET_DISPATCHER_NAME phenix.pdb_interpretation

def get_master_phil():
  import iotbx.phil
  return iotbx.phil.parse("""\
atom_selection = None
  .type = str
  .help = "Limit all analysis of restraints to this selection only."
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

include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
""", process_includes=True)

def run(args):
  from libtbx.str_utils import show_string
  #
  master_phil = get_master_phil()
  import iotbx.utils
  input_objects = iotbx.utils.process_command_line_inputs(
    args=args,
    master_phil=master_phil,
    input_types=("pdb", "cif"))
  work_phil = master_phil.fetch(sources=input_objects["phil"])
  work_phil.show()
  print
  work_params = work_phil.extract()
  #
  import mmtbx.monomer_library.server
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  for file_obj in input_objects["cif"]:
    print "Processing CIF file: %s" % show_string(file_obj.file_name)
    for srv in [mon_lib_srv, ener_lib]:
      srv.process_cif_object(
        cif_object=file_obj.file_content,
        file_name=file_obj.file_name)
  #
  import mmtbx.monomer_library.pdb_interpretation
  from libtbx.utils import Sorry
  processed_pdb_files = []
  for input_obj in input_objects["pdb"]:
    processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=work_params.pdb_interpretation,
      file_name=input_obj.file_name,
      atom_selection_string=work_params.atom_selection,
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
        params_remove=work_params.geometry_restraints.remove)
    if (work_params.build_xray_structure):
      processed_pdb_file.xray_structure()
    processed_pdb_files.append(processed_pdb_file)
  print
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
        f = open(geo_file_name, "w")
        print >> f, "# %sgeometry restraints for file:" % label
        print >> f, "#   %s" % show_string(pdb_file_name)
        from libtbx.utils import date_and_time
        print >> f, "#", date_and_time()
        print >> f
        geo.show_sorted(
          sites_cart=sites_cart, site_labels=site_labels, f=f)
        f.close()
      geo = processed_pdb_file.geometry_restraints_manager()
      if (work_params.write_geo_files):
        geo_file_name = op.basename(pdb_file_name) + ".geo"
        print "Writing file: %s" % show_string(geo_file_name)
        write_geo(label="", geo=geo, geo_file_name=geo_file_name)
        print
      if (work_params.write_tardy_geo_files):
        geo_file_name = op.basename(pdb_file_name) + ".tardy.geo"
        print "Writing file: %s" % show_string(geo_file_name)
        tardy_tree = geo.construct_tardy_tree(sites_cart=sites_cart)
        reduced_geo = geo.reduce_for_tardy(tardy_tree=tardy_tree)
        write_geo(label="tardy ", geo=reduced_geo, geo_file_name=geo_file_name)
        print

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
