
import sys, os

#-----------------------------------------------------------------------
# XXX this isn't really very useful as a user-space program - it is
# mostly just here for testing purposes
def run (args, out=sys.stdout) :
  from libtbx.utils import Sorry, Usage
  if (len(args) == 0) :
    print "Warning: this is not intended for general use."
    raise Usage("mmtbx.flip_peptides [model.pdb] [data.mtz] [params ...]")
  from mmtbx.refinement import flip_peptides, print_statistics
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx import utils
  from iotbx import reflection_file_utils
  import iotbx.phil
  import libtbx.phil.command_line
  from cStringIO import StringIO
  master_phil = iotbx.phil.parse("""
    input {
      xray_data {
        include scope mmtbx.utils.data_and_flags_str
      }
      pdb_file = None
        .type = path
    }
    flip_peptides {
      include scope mmtbx.refinement.flip_peptides.master_params_str
    }
  """, process_includes=True)
  cmdline = utils.process_command_line_args(
    args=args,
    master_params=master_phil)
  params = cmdline.params.extract()
  if params.input.pdb_file is None :
    if (len(cmdline.reflection_file_names) != 1) :
      raise Sorry("A single PDB file is required as input.")
    params.input.pdb_file = cmdline.pdb_file_names[0]
  reflection_file_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry=cmdline.crystal_symmetry,
    force_symmetry=True,
    reflection_files=cmdline.reflection_files,
    err=sys.stderr)
  data_and_flags = utils.determine_data_and_flags(
    reflection_file_server=reflection_file_server,
    parameters=params.input.xray_data,
    data_parameter_scope="input.xray_data",
    flags_parameter_scope="input.xray_data.r_free_flags")
  pdb_file = params.input.pdb_file
  cif_files = [ file_name for (file_name, cif_object) in cmdline.cif_objects ]
  processed_pdb_file = pdb_interpretation.run([pdb_file]+cif_files, log=out)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies=False)
  xray_structure = processed_pdb_file.xray_structure()
  chain_proxies = processed_pdb_file.all_chain_proxies
  pdb_hierarchy = chain_proxies.pdb_hierarchy
  fmodel = utils.fmodel_simple(
    xray_structures=[xray_structure],
    f_obs=data_and_flags.f_obs,
    r_free_flags=data_and_flags.r_free_flags)
  print >> out, ""
  print_statistics.make_header("Analyzing peptide bonds", out=out)
  flip_peptides.run(
    fmodel=fmodel,
    geometry_restraints_manager=geometry,
    pdb_hierarchy=pdb_hierarchy,
    solvent_selection=None)
  pdb_out = os.path.splitext(os.path.basename(pdb_file))[0] + "_new.pdb"
  open(pdb_out, "w").write(pdb_hierarchy.as_pdb_string())
  print >> out, ""
  print "wrote new PDB file to %s" % pdb_out

if __name__ == "__main__" :
  run(sys.argv[1:])
