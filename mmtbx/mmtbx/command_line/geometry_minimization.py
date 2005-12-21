from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from iotbx import pdb
import iotbx.pdb.interpretation
from cctbx import geometry_restraints
import cctbx.geometry_restraints.lbfgs
import scitbx.lbfgs
from libtbx.str_utils import show_string
from libtbx.utils import Sorry
from iotbx.option_parser import iotbx_option_parser
import sys, os

class lbfgs(geometry_restraints.lbfgs.lbfgs):

  def __init__(self,
        stage_1,
        sites_cart,
        geometry_restraints_manager,
        lbfgs_termination_params,
        intermediate_pdb_file_name_format,
        i_macro_cycle,
        write_pdb_rms_threshold,
        pymol_commands):
    self.stage_1 = stage_1
    self.intermediate_pdb_file_name_format = intermediate_pdb_file_name_format
    self.i_macro_cycle = i_macro_cycle
    self.write_pdb_rms_threshold = write_pdb_rms_threshold
    self.pymol_commands = pymol_commands
    self.sites_cart_written = sites_cart.deep_copy()
    geometry_restraints.lbfgs.lbfgs.__init__(self,
      sites_cart=sites_cart,
      geometry_restraints_manager=geometry_restraints_manager,
      lbfgs_termination_params=lbfgs_termination_params)

  def callback_after_step(self, minimizer):
    geometry_restraints.lbfgs.lbfgs.callback_after_step(self, minimizer)
    self.apply_shifts()
    incremental_rms_difference = self.tmp.sites_shifted.rms_difference(
      self.sites_cart_written)
    if (self.write_pdb_rms_threshold is not None
        and incremental_rms_difference >= self.write_pdb_rms_threshold):
      self.sites_cart_written = self.tmp.sites_shifted.deep_copy()
      total_rms_difference = self.tmp.sites_shifted.rms_difference(
        self.tmp.sites_cart)
      output_pdb_file_name = self.intermediate_pdb_file_name_format % (
        self.i_macro_cycle, minimizer.iter())
      output_pdb_file_object = open(output_pdb_file_name, "w")
      print "Writing:", output_pdb_file_name, "rms difference: %.6g" % (
        total_rms_difference)
      print >> output_pdb_file_object, "REMARK geometry minimization" \
        " macro cycle %d, step %d, rms difference: %.6g" % (
          self.i_macro_cycle, minimizer.iter(), total_rms_difference)
      self.stage_1.write_modified(
        out=output_pdb_file_object,
        new_sites_cart=self.tmp.sites_shifted,
        crystal_symmetry=self.tmp.geometry_restraints_manager.crystal_symmetry)
      output_pdb_file_object.close()
      self.pymol_commands.append("load %s, mov" % output_pdb_file_name)

def write_pdb_file(
      output_pdb_file_name,
      remark,
      stage_1,
      crystal_symmetry,
      sites_cart,
      pymol_commands):
  output_pdb_file_object = open(output_pdb_file_name, "w")
  print "Writing:", output_pdb_file_name
  print >> output_pdb_file_object, "REMARK geometry minimization:", remark
  stage_1.write_modified(
    out=output_pdb_file_object,
    new_sites_cart=sites_cart,
    crystal_symmetry=crystal_symmetry)
  output_pdb_file_object.close()
  pymol_commands.append("load %s, mov" % output_pdb_file_name)
  print

def write_pymol_commands(file_name, commands):
  print "Writing:", file_name
  f = open(file_name, "w")
  print >> f, commands[-1]
  print >> f, commands[0]+", 1"
  for cmd in commands[1:]:
    print >> f, cmd
  print >> f, "mplay"
  f.close()
  print

def run(args, this_command="mmtbx.geometry_minimization"):
  command_line = (iotbx_option_parser(
    usage=this_command+" [options] pdb_file [output_pdb_file]")
    .enable_symmetry_comprehensive()
    .option(None, "--max_iterations",
      action="store",
      type="int",
      dest="max_iterations",
      default=500,
      metavar="INT")
    .option(None, "--macro_cycles",
      action="store",
      type="int",
      dest="macro_cycles",
      default=1,
      metavar="INT")
    .option(None, "--write_pdb_rms_threshold",
      action="store",
      type="float",
      dest="write_pdb_rms_threshold",
      default=None,
      metavar="FLOAT")
    .option(None, "--show_geometry_restraints",
      action="store_true",
      dest="show_geometry_restraints")
  ).process(args=args, min_nargs=1, max_nargs=2)
  input_pdb_file_name = None
  output_pdb_file_name = None
  cif_objects = []
  arg_is_processed = False
  for arg in command_line.args:
    if (not os.path.isfile(arg)):
      if (output_pdb_file_name is not None):
        raise Sorry(
          "More than one output PDB file name:\n"
          "  %s\n"
          "  %s\n" % (
            show_string(output_pdb_file_name),
            show_string(arg)))
      output_pdb_file_name = arg
      arg_is_processed = True
    else:
      if (pdb.interpretation.is_pdb_file(file_name=arg)):
        if (input_pdb_file_name is not None):
          raise Sorry(
            "More than one input PDB file name:\n"
            "  %s\n"
            "  %s\n" % (
              show_string(input_pdb_file_name),
              show_string(arg)))
        input_pdb_file_name = arg
        arg_is_processed = True
      else:
        try: cif_object = mmtbx.monomer_library.server.read_cif(file_name=arg)
        except KeyboardInterrupt: raise
        except: pass
        else:
          if (len(cif_object) > 0):
            cif_objects.append((arg,cif_object))
            arg_is_processed = True
    if (not arg_is_processed):
      raise Sorry(
        "Command line argument not recognized: %s" %
          show_string(arg))
  log = sys.stdout
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  for file_name,cif_object in cif_objects:
    mon_lib_srv.process_cif_object(cif_object=cif_object, file_name=file_name)
    ener_lib.process_cif_object(cif_object=cif_object, file_name=file_name)
  del cif_objects
  if (output_pdb_file_name is None):
    output_pdb_file_name = list(os.path.splitext(os.path.basename(
      input_pdb_file_name)))
    output_file_name_prefix = output_pdb_file_name[0]
    intermediate_pdb_file_name_format = "".join([
      output_pdb_file_name[0] + "_macro_%03d_step_%05d",
      output_pdb_file_name[1]])
    output_pdb_file_name[0] += "_geometry_minimized"
    output_pdb_file_name = "".join(output_pdb_file_name)
  output_pdb_file_object = open(output_pdb_file_name, "w")
  all_chain_proxies=monomer_library.pdb_interpretation.build_all_chain_proxies(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=input_pdb_file_name,
    crystal_symmetry=command_line.symmetry,
    force_symmetry=True,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=log)
  geometry_restraints_manager = \
    all_chain_proxies.construct_geometry_restraints_manager(
      ener_lib=ener_lib,
      disulfide_link=mon_lib_srv.link_link_id_dict["SS"],
      log=log)
  special_position_settings = all_chain_proxies.special_position_settings
  sites_cart = all_chain_proxies.sites_cart_exact().deep_copy()
  atom_labels = [atom.pdb_format()
    for atom in all_chain_proxies.stage_1.atom_attributes_list]
  geometry_restraints_manager.site_symmetry_table \
    .show_special_position_shifts(
      special_position_settings=special_position_settings,
      site_labels=atom_labels,
      sites_cart_original=all_chain_proxies.sites_cart,
      sites_cart_exact=sites_cart,
      out=log,
      prefix="  ")
  if (command_line.options.show_geometry_restraints):
    geometry_restraints_manager.show_interactions(
      sites_cart=sites_cart,
      site_labels=atom_labels)
  geometry_restraints_manager.pair_proxies(
    sites_cart=all_chain_proxies.sites_cart)\
      .bond_proxies.show_sorted_by_residual(
        sites_cart=sites_cart,
        labels=atom_labels,
        f=log,
        max_lines=10)
  print
  log.flush()
  pymol_commands = []
  if (command_line.options.write_pdb_rms_threshold is not None):
    write_pdb_file(
      output_pdb_file_name=output_file_name_prefix+"_initial.pdb",
      remark="initial coordinates",
      stage_1=all_chain_proxies.stage_1,
      crystal_symmetry=geometry_restraints_manager.crystal_symmetry,
      sites_cart=sites_cart,
      pymol_commands=pymol_commands)
  for i_macro_cycle in xrange(command_line.options.macro_cycles):
    minimized = lbfgs(
      stage_1=all_chain_proxies.stage_1,
      sites_cart=sites_cart,
      geometry_restraints_manager=geometry_restraints_manager,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=command_line.options.max_iterations),
      intermediate_pdb_file_name_format=intermediate_pdb_file_name_format,
      i_macro_cycle=i_macro_cycle,
      write_pdb_rms_threshold=command_line.options.write_pdb_rms_threshold,
      pymol_commands=pymol_commands)
    print "Energies at start of minimization:"
    minimized.first_target_result.show()
    print
    print "Number of minimization iterations:", minimized.minimizer.iter()
    print "Root-mean-square coordinate difference: %.3f" % (
      all_chain_proxies.sites_cart.rms_difference(sites_cart))
    print
    print "Energies at end of minimization:"
    minimized.final_target_result.show()
    print
    geometry_restraints_manager.pair_proxies(
      sites_cart=all_chain_proxies.sites_cart)\
        .bond_proxies.show_sorted_by_residual(
          sites_cart=sites_cart,
          labels=atom_labels,
          f=log,
          max_lines=10)
    print
  print "Writing:", output_pdb_file_name
  print >> output_pdb_file_object, "REMARK", this_command, " ".join(args)
  all_chain_proxies.stage_1.write_modified(
    out=output_pdb_file_object,
    new_sites_cart=sites_cart,
    crystal_symmetry=special_position_settings)
  output_pdb_file_object.close()
  print
  if (len(pymol_commands) > 0):
    pymol_commands.append("load %s, mov" % output_pdb_file_name)
    write_pymol_commands(
      os.path.splitext(output_pdb_file_name)[0]+".pml",
      pymol_commands)

if (__name__ == "__main__"):
  run(sys.argv[1:])
