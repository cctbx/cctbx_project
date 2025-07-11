"""Regularize model geometry against existing restraints"""
# LIBTBX_SET_DISPATCHER_NAME phenix.geometry_minimization

from __future__ import absolute_import, division, print_function
import mmtbx.refinement.geometry_minimization
import mmtbx.utils
from iotbx.pdb import combine_unique_pdb_files
import iotbx.phil
from cctbx.array_family import flex
from libtbx.utils import user_plus_sys_time, Sorry
from libtbx import runtime_utils
import os
import sys
from six.moves import cStringIO as StringIO
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.hydrogens import riding
import mmtbx.model
from cctbx import uctbx

base_params_str = """\
silent = False
  .type = bool
write_geo_file = True
  .type = bool
file_name = None
  .type = path
  .short_caption = Model file
  .style = file_type:pdb bold input_file
show_states = False
  .type = bool
restraints = None
  .type = path
  .multiple = True
  .short_caption = Restraints
  .style = file_type:cif bold input_file
restraints_directory = None
  .type = path
  .style = directory
output_file_name_prefix = None
  .type = str
  .input_size = 400
  .style = bold
directory = None
  .type = path
  .short_caption = Output directory
  .style = output_dir
include scope libtbx.phil.interface.tracking_params
fix_rotamer_outliers = True
  .type = bool
  .help = Remove outliers
allow_allowed_rotamers = True
  .type = bool
  .help = More strict fixing outliers
stop_for_unknowns = True
  .type = bool
  .short_caption = Stop for unknown residues
  .style = noauto
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
include scope \
    mmtbx.geometry_restraints.torsion_restraints.reference_model.reference_model_params
"""

master_params_str = """
%s
selection = all
  .type = str
  .help = Atom selection string: selected atoms are subject to move
  .short_caption = Atom selection
  .input_size = 400
minimization
  .help = Geometry minimization parameters
  .short_caption = Minimization parameters
  .expert_level=1
{
  max_iterations = 500
    .type = int
    .help = Maximun number of minimization iterations
    .short_caption = Max. iterations
    .style = noauto
  macro_cycles = 5
    .type = int
    .help = Number of minimization macro-cycles
  alternate_nonbonded_off_on = False
    .type = bool
    .short_caption = Macro cycles
    .style = noauto
  rmsd_bonds_termination_cutoff = 0
    .type = float
    .help = stop after reaching specified cutoff value
  rmsd_angles_termination_cutoff = 0
    .type = float
    .help = stop after reaching specified cutoff value
  grms_termination_cutoff = 0
    .type = float
    .help = stop after reaching specified cutoff value
  correct_special_position_tolerance = 1.0
    .type = float
  riding_h = True
    .type = bool
    .help = Use riding model for H
  move
    .help = Define what to include into refinement target
    .short_caption = Geometry terms
    .style = box auto_align columns:4 noauto
  {
  bond = True
    .type = bool
    .short_caption = Bond lengths
  nonbonded = True
    .type = bool
    .short_caption = Nonbonded distances
  angle = True
    .type = bool
    .short_caption = Bond angle
  dihedral = True
    .type = bool
    .short_caption = Dihedral angle
  chirality = True
    .type = bool
    .short_caption = Chirality
  planarity = True
    .type = bool
    .short_caption = Planarity
  parallelity = True
    .type = bool
    .short_caption = Parallelity
  }
}
  include scope mmtbx.geometry_restraints.external.external_energy_params_str
""" % base_params_str

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def broadcast(m, log):
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

def format_usage_message(log):
  print("-"*79, file=log)
  msg = """\
phenix.geometry_minimization: regularize model geometry

Usage examples:
  phenix.geometry_minimization model.pdb
  phenix.geometry_minimization model.pdb ligands.cif
"""
  print(msg, file=log)
  print("-"*79, file=log)

def run_minimization(
      selection,
      restraints_manager,
      riding_h_manager,
      pdb_hierarchy,
      params,
      cdl,
      rdl,
      correct_hydrogens,
      states_collector,
      fix_rotamer_outliers,
      allow_allowed_rotamers,
      log,
      ncs_restraints_group_list = [],
      mon_lib_srv = None):
  o = mmtbx.refinement.geometry_minimization.run2(
    restraints_manager             = restraints_manager,
    riding_h_manager               = riding_h_manager,
    pdb_hierarchy                  = pdb_hierarchy,
    ncs_restraints_group_list      = ncs_restraints_group_list,
    max_number_of_iterations       = params.max_iterations,
    number_of_macro_cycles         = params.macro_cycles,
    selection                      = selection,
    correct_special_position_tolerance = params.correct_special_position_tolerance,
    bond                           = params.move.bond,
    nonbonded                      = params.move.nonbonded,
    angle                          = params.move.angle,
    dihedral                       = params.move.dihedral,
    chirality                      = params.move.chirality,
    planarity                      = params.move.planarity,
    parallelity                    = params.move.parallelity,
    rmsd_bonds_termination_cutoff  = params.rmsd_bonds_termination_cutoff,
    rmsd_angles_termination_cutoff = params.rmsd_angles_termination_cutoff,
    alternate_nonbonded_off_on     = params.alternate_nonbonded_off_on,
    cdl                            = cdl,
    rdl                            = rdl,
    states_collector               = states_collector,
    correct_hydrogens              = correct_hydrogens,
    fix_rotamer_outliers           = fix_rotamer_outliers,
    allow_allowed_rotamers         = allow_allowed_rotamers,
    log                            = log,
    mon_lib_srv                    = mon_lib_srv)

def run_minimization_amber(
      selection,
      restraints_manager,
      pdb_hierarchy,
      params,
      log,
      prmtop,
      ambcrd,
      ):
  import amber_adaptbx.amber_geometry_minimization
  o = amber_adaptbx.amber_geometry_minimization.run(
    restraints_manager             = restraints_manager,
    pdb_hierarchy = pdb_hierarchy,
    max_number_of_iterations       = params.max_iterations,
    number_of_macro_cycles         = params.macro_cycles,
    selection                      = selection,
    bond                           = params.move.bond,
    nonbonded                      = params.move.nonbonded,
    angle                          = params.move.angle,
    dihedral                       = params.move.dihedral,
    chirality                      = params.move.chirality,
    planarity                      = params.move.planarity,
    parallelity                    = params.move.parallelity,
    grms_termination_cutoff       = params.grms_termination_cutoff,
    alternate_nonbonded_off_on     = params.alternate_nonbonded_off_on,
    log                            = log,
    prmtop                         = prmtop,
    ambcrd                         = ambcrd,
    )

class run(object):
  _pdb_suffix = "minimized"
  def __init__(self, args, log, use_directory_prefix=True):
    # You are not supposed to put here (in __init__) any time-consuming stuff,
    # otherwise self.total_time would be unaccurate. It's not clear
    # why it is important.
    self.model                = None
    self.log                  = log
    self.params               = None
    self.inputs               = None
    self.args                 = args
    self.selection            = None
    self.restrain_selection   = None
    self.time_strings         = []
    self.total_time           = 0
    self.pdb_file_names       = []
    self.use_directory_prefix = use_directory_prefix
    self.sites_cart_start     = None
    self.states_collector     = None
    self.__execute()

  def __execute(self):
    #
    self.caller(self.initialize,            "Initialization, inputs")
    self.caller(self.process_inputs,        "Processing inputs")
    self.caller(self.atom_selection,        "Atom selection")
    self.caller(self.get_restraints,        "Geometry Restraints")
    self.caller(self.setup_riding_h,        "Setup riding H")
    self.caller(self.minimization,          "Minimization")
    self.caller(self.write_pdb_file,        "Write PDB file")
    self.caller(self.write_geo_file,        "Write GEO file")
    self.caller(self.show_model_statistics, "Model statistics")
    #
    self.show_times()

  def master_params(self):
    return master_params()

  def caller(self, func, prefix):
    timer = user_plus_sys_time()
    func(prefix = prefix)
    t = timer.elapsed()
    self.total_time += t
    self.time_strings.append("  %s: %s"%(prefix, str("%8.3f"%t).strip()))

  def show_times(self):
    broadcast(m="Detailed timing", log = self.log)
    max_len = 0
    for ts in self.time_strings:
      lts = len(ts)
      if(lts > max_len): max_len = lts
    fmt = "  %-"+str(lts)+"s"
    for ts in self.time_strings:
      sts = ts.split()
      l = " ".join(sts[:len(sts)-1])
      print(fmt%l, sts[len(sts)-1], file=self.log)
    print("  Sum of individual times: %s"%\
      str("%8.3f"%self.total_time).strip(), file=self.log)

  def format_usage_message(self):
    format_usage_message(log=self.log)

  def setup_output_file_names(self):
    # for pdb
    ofn = self.params.output_file_name_prefix
    directory = self.params.directory
    base_name = ""
    if self.use_directory_prefix and directory is not None:
      base_name = directory
    suffix = "_" + self._pdb_suffix  + ".pdb"
    if self.params.output_file_name_prefix is None:
      in_fn = os.path.basename(self.pdb_file_names[0])
      ind = max(0, in_fn.rfind("."))
      ofn = in_fn + suffix
      if ind > 0:
        ofn = in_fn[:ind]+suffix
    else:
      ofn = self.params.output_file_name_prefix+".pdb"
    self.result_model_fname = os.path.join(base_name, ofn)
    self.result_states_fname = self.result_model_fname[:].replace(".pdb","_all_states.pdb")
    self.final_geo_fname = self.result_model_fname[:].replace(".pdb",".geo")

  def initialize(self, prefix):
    if (self.log is None) : self.log = sys.stdout
    if(len(self.args)==0):
      self.format_usage_message()
    parsed = self.master_params()
    self.inputs = mmtbx.utils.process_command_line_args(args = self.args,
      master_params = parsed)
    self.params = self.inputs.params.extract()
    if(self.params.silent): self.log = StringIO()
    broadcast(m=prefix, log = self.log)
    self.inputs.params.show(prefix="  ", out=self.log)
    if(len(self.args)==0): sys.exit(0)

  def process_inputs(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.pdb_file_names = list(self.inputs.pdb_file_names)
    if(self.params.file_name is not None):
      self.pdb_file_names.append(self.params.file_name)

    #=================================================
    cs = self.inputs.crystal_symmetry
    is_non_crystallographic_unit_cell = False
    import iotbx.pdb
    pdb_combined = combine_unique_pdb_files(file_names = self.pdb_file_names)
    pdb_inp = iotbx.pdb.input(lines=pdb_combined.raw_records, source_info=None)
    if(cs is None):
      cs=pdb_inp.crystal_symmetry()
    if(cs is None):
      is_non_crystallographic_unit_cell = True
      box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
        sites_cart   = pdb_inp.atoms().extract_xyz(),
        buffer_layer = 10)
      cs = box.crystal_symmetry()
    cif_objects = list(self.inputs.cif_objects)
    if (len(self.params.restraints) > 0):
      import iotbx.cif
      for file_name in self.params.restraints :
        cif_object = iotbx.cif.reader(file_path=file_name, strict=False).model()
        cif_objects.append((file_name, cif_object))
    if (self.params.restraints_directory is not None):
      restraint_files = os.listdir(self.params.restraints_directory)
      for file_name in restraint_files :
        if (file_name.endswith(".cif")):
          full_path = os.path.join(self.params.restraints_directory, file_name)
          cif_object = iotbx.cif.reader(file_path=full_path,
            strict=False).model()
          cif_objects.append((full_path, cif_object))
    self.model = mmtbx.model.manager(
        model_input = pdb_inp,
        crystal_symmetry = cs,
        restraint_objects = cif_objects,
        stop_for_unknowns = self.params.stop_for_unknowns,
        log = self.log)
    self.model.process(
      pdb_interpretation_params = self.params,
      make_restraints           = True)
    self.ncs_obj = self.model.get_ncs_obj()
    self.output_crystal_symmetry = not is_non_crystallographic_unit_cell
    self.sites_cart_start = self.model.get_xray_structure().sites_cart().deep_copy()
    if(self.params.show_states):
      self.states_collector = mmtbx.utils.states(
        xray_structure = self.model.get_xray_structure(),
        pdb_hierarchy  = self.model.get_hierarchy())
    self.setup_output_file_names()

  def atom_selection(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.selection = self.model.selection(string = self.params.selection)
    print("  selected %s atoms out of total %s"%(
      str(self.selection.count(True)),str(self.selection.size())), file=self.log)

  def get_restraints(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.model.get_restraints_manager()

  def setup_riding_h(self, prefix):
    if not self.params.minimization.riding_h: return
    broadcast(m=prefix, log = self.log)
    self.model.setup_riding_h_manager(idealize=True)

  def minimization(self, prefix): # XXX USE alternate_nonbonded_off_on etc
    broadcast(m=prefix, log = self.log)
    use_amber = False
    if self.ncs_obj is not None:
      print("Using NCS constraints:", file=self.log)
      self.ncs_obj.show(format='phil', log=self.log)
    ncs_restraints_group_list = []
    if self.ncs_obj is not None:
      ncs_restraints_group_list = self.ncs_obj.get_ncs_restraints_group_list()
    run_minimization(
      selection              = self.selection,
      restraints_manager     = self.model.get_restraints_manager(),
      riding_h_manager       = self.model.get_riding_h_manager(),
      params                 = self.params.minimization,
      pdb_hierarchy          = self.model.get_hierarchy(),
      cdl                    = self.params.pdb_interpretation.restraints_library.cdl,
      rdl                    = self.params.pdb_interpretation.restraints_library.rdl,
      correct_hydrogens      = self.params.pdb_interpretation.correct_hydrogens,
      fix_rotamer_outliers   = self.params.fix_rotamer_outliers,
      allow_allowed_rotamers = self.params.allow_allowed_rotamers,
      states_collector       = self.states_collector,
      log                    = self.log,
      ncs_restraints_group_list = ncs_restraints_group_list,
      mon_lib_srv            = self.model.get_mon_lib_srv())
    self.model.set_sites_cart_from_hierarchy()

  def write_pdb_file(self, prefix):
    broadcast(m=prefix, log = self.log)
    # self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    print(self.min_max_mean_shift(), file=self.log)

    if self.model.can_be_output_as_pdb():
      print("  output file name:", self.result_model_fname, file=self.log)
      r = self.model.model_as_pdb(output_cs=self.output_crystal_symmetry)
      f = open(self.result_model_fname, 'w')
      f.write(r)
      f.close()

    cif_fname = self.result_model_fname.replace('.pdb', '.cif')
    if not os.path.isfile(cif_fname):
      print("  output file name:", cif_fname, file=self.log)
      r = self.model.model_as_mmcif(output_cs=self.output_crystal_symmetry)
      with open(cif_fname, 'w') as f:
        f.write(r)

    if(self.states_collector):
      self.states_collector.write(
        file_name=self.result_states_fname)

  def min_max_mean_shift(self):
    return "min,max,mean shift from start: %6.3f %6.3f %6.3f"%flex.sqrt((
      self.sites_cart_start - self.model.get_xray_structure().sites_cart()).dot()
      ).min_max_mean().as_tuple()

  def write_geo_file(self, prefix):
    if self.params.write_geo_file:
      broadcast(m=prefix, log = self.log)
      # no output of NCS stuff here
      restr_txt = self.model.restraints_as_geo()
      f = open(self.final_geo_fname, "w")
      f.write("# Geometry restraints after refinement\n")
      f.write(restr_txt)
      f.close()

  def show_model_statistics(self, prefix):
    if self.params.write_geo_file:
      broadcast(m=prefix, log = self.log)
      s = self.model.geometry_statistics()
      s.show(log = self.log, uppercase=False)

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    filename = run(args=self.args, log=sys.stdout,
                   use_directory_prefix=False).result_model_fname
    return os.path.join(self.output_dir, filename)

def validate_params(params):
  if (params.file_name is None):
    raise Sorry("Please specify a model file to minimize.")
  if (params.restraints_directory is not None):
    if (not os.path.isdir(params.restraints_directory)):
      raise Sorry("The path '%s' does not exist or is not a directory." %
        params.restraints_directory)
  return True

def finish_job(result):
  output_files = []
  if (result is not None):
    output_files.append((result, "Minimized model"))
  return output_files, []

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  log = sys.stdout
  o = run(sys.argv[1:], log=log)
  tt = timer.elapsed()
  print("Overall runtime: %-8.3f" % tt, file=o.log)
  assert abs(tt-o.total_time) < 0.1 # guard against unaccounted times

