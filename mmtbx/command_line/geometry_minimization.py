# LIBTBX_SET_DISPATCHER_NAME phenix.geometry_minimization

from __future__ import division
import mmtbx.refinement.geometry_minimization
import mmtbx.utils
from iotbx.pdb import combine_unique_pdb_files
import iotbx.phil
from cctbx.array_family import flex
from libtbx.utils import user_plus_sys_time, Sorry
from libtbx import runtime_utils
import os
import sys
from cStringIO import StringIO
from mmtbx.monomer_library import pdb_interpretation

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
stop_for_unknowns = True
  .type = bool
  .short_caption = Stop for unknown residues
  .style = noauto
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
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
  move
    .help = Define what to include into refinement target
    .short_caption = Geometry terms
    .style = box auto_align columns:3 noauto
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
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def format_usage_message(log):
  print >> log, "-"*79
  msg = """\
phenix.geometry_minimization: regularize model geometry

Usage examples:
  phenix.geometry_minimization model.pdb
  phenix.geometry_minimization model.pdb ligands.cif
"""
  print >> log, msg
  print >> log, "-"*79

def process_input_files(inputs, params, log):
  pdb_file_names = []
  pdb_file_names = list(inputs.pdb_file_names)
  if (params.file_name is not None) :
    pdb_file_names.append(params.file_name)
  cs = inputs.crystal_symmetry
  is_non_crystallographic_unit_cell = False
  if(cs is None):
    is_non_crystallographic_unit_cell = True
    import iotbx.pdb
    pdb_combined = combine_unique_pdb_files(file_names = pdb_file_names)
    cs = iotbx.pdb.input(source_info = None, lines = flex.std_string(
      pdb_combined.raw_records)).xray_structure_simple().\
        cubic_unit_cell_around_centered_scatterers(
        buffer_size = 10).crystal_symmetry()
  cif_objects = list(inputs.cif_objects)
  if (len(params.restraints) > 0) :
    import iotbx.cif
    for file_name in params.restraints :
      cif_object = iotbx.cif.reader(file_path=file_name, strict=False).model()
      cif_objects.append((file_name, cif_object))
  if (params.restraints_directory is not None) :
    restraint_files = os.listdir(params.restraints_directory)
    for file_name in restraint_files :
      if (file_name.endswith(".cif")) :
        full_path = os.path.join(params.restraints_directory, file_name)
        cif_object = iotbx.cif.reader(file_path=full_path,
          strict=False).model()
        cif_objects.append((full_path, cif_object))
  processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
    crystal_symmetry          = cs,
    pdb_interpretation_params = params.pdb_interpretation,
    stop_for_unknowns         = params.stop_for_unknowns,
    log                       = log,
    cif_objects               = cif_objects,
    use_neutron_distances     = params.pdb_interpretation.use_neutron_distances)
  processed_pdb_file, junk = processed_pdb_files_srv.\
    process_pdb_files(pdb_file_names = pdb_file_names) # XXX remove junk
  processed_pdb_file.is_non_crystallographic_unit_cell = \
    is_non_crystallographic_unit_cell # XXX bad hack
  return processed_pdb_file

def get_geometry_restraints_manager(processed_pdb_file, xray_structure,
    params=None, log=sys.stdout):
  """This function should be transfered to be a member of processed_pdb_file
     class. It should be used in all places where geometry_restraints_manager
     is needed. """
  has_hd = None
  if(xray_structure is not None):
    sctr_keys = xray_structure.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
  reference_torsion_proxies = None
  if params is None:
    params = master_params().fetch().extract()
  # disabled temporarily due to architecture changes
  """
  id_params = params.secondary_structure.idealization
  from mmtbx.secondary_structure.build import substitute_ss
  if id_params.enabled:
    print >> log, "Substituting secondary structure elements with ideal ones:"
    annot = None
    if len(params.secondary_structure.helix) +\
       len(params.secondary_structure.sheet) >0:
      annot = iotbx.pdb_secondary_structure.annotation.from_phil(
        phil_helices=params.secondary_structure.helix,
        phil_sheets=params.secondary_structure.sheet,
        pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
        log=log)
    else:
      annot = processed_pdb_file.all_chain_proxies.pdb_inp.\
          extract_secondary_structure()
    if annot is not None:
      reference_torsion_proxies = substitute_ss(
          real_h=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
          xray_structure=xray_structure,
          ss_annotation=annot,
          sigma_on_reference_non_ss=id_params.sigma_on_reference_non_ss,
          sigma_on_reference_helix=id_params.sigma_on_reference_helix,
          sigma_on_reference_sheet=id_params.sigma_on_reference_sheet,
          sigma_on_torsion_ss=id_params.sigma_on_torsion_ss,
          sigma_on_torsion_nonss=id_params.sigma_on_torsion_nonss,
          sigma_on_ramachandran=id_params.sigma_on_ramachandran,
          sigma_on_cbeta=id_params.sigma_on_cbeta,
          n_macro=id_params.n_macro,
          n_iter=id_params.n_iter,
          log=log,
          verbose=False)
      xray_structure.set_sites_cart(
          processed_pdb_file.all_chain_proxies.pdb_hierarchy.\
              atoms().extract_xyz())
      print >> log, "Substituting secondary structure elements is done"
    else:
      print >> log, "No secondary structure definition found in phil"+\
          "or HELIX/SHEET records. No substitution done."
  """
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    params_edits                 = params.geometry_restraints.edits,
    plain_pairs_radius           = 5,
    assume_hydrogens_all_missing = not has_hd)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)
  if (reference_torsion_proxies is not None
      and id_params.restrain_torsion_angles):
    geometry.add_dihedrals_in_place(
        additional_dihedral_proxies=reference_torsion_proxies)
  if(xray_structure is not None):
    restraints_manager.crystal_symmetry = xray_structure.crystal_symmetry()
  return restraints_manager

def run_minimization(
      selection,
      restraints_manager,
      pdb_hierarchy,
      params,
      cdl,
      correct_hydrogens,
      states_collector,
      fix_rotamer_outliers,
      log):
  o = mmtbx.refinement.geometry_minimization.run2(
    restraints_manager             = restraints_manager,
    pdb_hierarchy                  = pdb_hierarchy,
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
    states_collector               = states_collector,
    correct_hydrogens              = correct_hydrogens,
    fix_rotamer_outliers           = fix_rotamer_outliers,
    log                            = log)

def run_minimization_amber (
      selection,
      restraints_manager,
      pdb_hierarchy,
      params,
      log,
      prmtop,
      ambcrd,
      use_sander):
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
    use_sander                     = use_sander)

class run(object):
  _pdb_suffix = "minimized"
  def __init__(self, args, log, use_directory_prefix=True):
    self.log                  = log
    self.params               = None
    self.inputs               = None
    self.args                 = args
    self.processed_pdb_file   = None
    self.xray_structure       = None
    self.pdb_hierarchy        = None
    self.selection            = None
    self.restrain_selection   = None
    self.grm                  = None
    self.time_strings         = []
    self.total_time           = 0
    self.output_file_name     = None
    self.pdb_file_names       = []
    self.use_directory_prefix = use_directory_prefix
    self.sites_cart_start     = None
    self.states_collector     = None
    self.__execute()

  def __execute(self):
    #
    self.caller(self.initialize,           "Initialization, inputs")
    self.caller(self.process_inputs,       "Processing inputs")
    self.caller(self.atom_selection,       "Atom selection")
    self.caller(self.get_restraints,       "Geometry Restraints")
    self.caller(self.minimization,         "Minimization")
    self.caller(self.write_pdb_file,       "Write PDB file")
    self.caller(self.write_geo_file,       "Write GEO file")
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
      print >> self.log, fmt%l, sts[len(sts)-1]
    print >> self.log, "  Sum of individual times: %s"%\
      str("%8.3f"%self.total_time).strip()

  def format_usage_message (self) :
    format_usage_message(log=self.log)

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
    self.processed_pdb_file = process_input_files(inputs=self.inputs,
      params=self.params, log=self.log)
    self.output_crystal_symmetry = \
      not self.processed_pdb_file.is_non_crystallographic_unit_cell
    self.xray_structure = self.processed_pdb_file.xray_structure()
    self.sites_cart_start = self.xray_structure.sites_cart().deep_copy()
    self.pdb_hierarchy = self.processed_pdb_file.all_chain_proxies.pdb_hierarchy
    if(self.params.show_states):
      self.states_collector = mmtbx.utils.states(
        xray_structure = self.xray_structure,
        pdb_hierarchy  = self.pdb_hierarchy)

  def atom_selection(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.selection = mmtbx.utils.atom_selection(
      all_chain_proxies = self.processed_pdb_file.all_chain_proxies,
      string = self.params.selection)
    print >> self.log, "  selected %s atoms out of total %s"%(
      str(self.selection.count(True)),str(self.selection.size()))

  def get_restraints(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.grm = get_geometry_restraints_manager(
      processed_pdb_file = self.processed_pdb_file,
      xray_structure     = self.xray_structure,
      params             = self.params,
      log                = self.log)

  def minimization(self, prefix): # XXX USE alternate_nonbonded_off_on etc
    broadcast(m=prefix, log = self.log)
    use_amber = False
    if hasattr(self.params, "amber"):
      use_amber = self.params.amber.use_amber
    if(use_amber):
      if not self.params.amber.topology_file_name:
        raise Sorry("Need to supply topology file using amber.topology_file_name=<filename>")
      if not self.params.amber.coordinate_file_name:
        raise Sorry("Need to supply topology file using amber.coordinate_file_name=<filename>")
      run_minimization_amber(
        selection = self.selection,
        restraints_manager = self.grm,
        params = self.params.minimization,
        pdb_hierarchy = self.pdb_hierarchy,
        log = self.log,
        prmtop = self.params.amber.topology_file_name,
        ambcrd = self.params.amber.coordinate_file_name,
        use_sander = self.params.amber.use_sander)
    else:
      run_minimization(
        selection            = self.selection,
        restraints_manager   = self.grm,
        params               = self.params.minimization,
        pdb_hierarchy        = self.pdb_hierarchy,
        cdl                  = self.params.pdb_interpretation.cdl,
        correct_hydrogens    = self.params.pdb_interpretation.correct_hydrogens,
        fix_rotamer_outliers = self.params.fix_rotamer_outliers,
        states_collector= self.states_collector,
        log                  = self.log)
    self.xray_structure.set_sites_cart(
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz())

  def write_pdb_file(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    ofn = self.params.output_file_name_prefix
    directory = self.params.directory
    suffix = "_" + self._pdb_suffix  + ".pdb"
    if(ofn is None):
      pfn = os.path.basename(self.pdb_file_names[0])
      ind = max(0,pfn.rfind("."))
      ofn = pfn+suffix if ind==0 else pfn[:ind]+suffix
    else: ofn = self.params.output_file_name_prefix+".pdb"
    if (self.use_directory_prefix) and (directory is not None) :
      ofn = os.path.join(directory, ofn)
    print >> self.log, "  output file name:", ofn
    print >> self.log, self.min_max_mean_shift()
    if (self.output_crystal_symmetry) :
      self.pdb_hierarchy.write_pdb_file(file_name = ofn, crystal_symmetry =
        self.xray_structure.crystal_symmetry())
    else :
      self.pdb_hierarchy.write_pdb_file(file_name = ofn)
    if(self.states_collector):
      self.states_collector.write(
        file_name=ofn[:].replace(".pdb","_all_states.pdb"))
    self.output_file_name = os.path.abspath(ofn)

  def min_max_mean_shift(self):
    return "min,max,mean shift from start: %6.3f %6.3f %6.3f"%flex.sqrt((
      self.sites_cart_start - self.xray_structure.sites_cart()).dot()
      ).min_max_mean().as_tuple()

  def write_geo_file(self, prefix):
    if(self.params.write_geo_file and self.grm is not None):
      broadcast(m=prefix, log = self.log)
      ofn = os.path.basename(self.output_file_name).replace(".pdb",".geo")
      directory = self.params.directory
      if (self.use_directory_prefix) and (directory is not None) :
        ofn = os.path.join(directory, ofn)
      # no output of NCS stuff here
      self.grm.write_geo_file(
          file_name=ofn,
          header="# Geometry restraints after refinement\n",
          xray_structure=self.xray_structure)

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=self.args, log=sys.stdout,
      use_directory_prefix=False).output_file_name

def validate_params (params) :
  if (params.file_name is None) :
    raise Sorry("Please specify a model file to minimize.")
  if (params.restraints_directory is not None) :
    if (not os.path.isdir(params.restraints_directory)) :
      raise Sorry("The path '%s' does not exist or is not a directory." %
        params.restraints_directory)
  return True

def finish_job (result) :
  output_files = []
  if (result is not None) :
    output_files.append((result, "Minimized model"))
  return output_files, []

if(__name__ == "__main__"):
  timer = user_plus_sys_time()
  log = sys.stdout
  o = run(sys.argv[1:], log=log)
  tt = timer.elapsed()
  print >> o.log, "Overall runtime: %-8.3f" % tt
  assert abs(tt-o.total_time) < 0.1 # guard against unaccounted times
