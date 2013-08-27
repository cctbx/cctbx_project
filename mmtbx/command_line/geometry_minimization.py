# LIBTBX_SET_DISPATCHER_NAME phenix.geometry_minimization

from __future__ import division
import mmtbx.refinement.geometry_minimization
from mmtbx import monomer_library
import mmtbx.utils
from iotbx.pdb import combine_unique_pdb_files
import iotbx.phil
from cctbx.array_family import flex
from libtbx.utils import user_plus_sys_time, Sorry
from libtbx import runtime_utils
import os
import mmtbx.secondary_structure
import sys
from cStringIO import StringIO

base_params_str = """\
silent = False
  .type = bool
write_geo_file = True
  .type = bool
file_name = None
  .type = path
  .short_caption = Model file
  .style = file_type:pdb bold input_file
restraints = None
  .type = path
  .multiple = True
  .short_caption = Restraints
  .style = file_type:cif bold input_file
restraints_directory = None
  .type = path
  .style = directory
use_neutron_distances = False
  .type = bool
  .help = Use neutron X-H distances (which are longer than X-ray ones)
  .short_caption = Use neutron X-H distances
output_file_name_prefix = None
  .type = str
  .input_size = 400
  .style = bold
directory = None
  .type = path
  .short_caption = Output directory
  .style = output_dir
include scope libtbx.phil.interface.tracking_params
reference_restraints {
  restrain_starting_coord_selection = None
    .type = str
    .help = Atom selection string: restraint selected to starting position
    .short_caption = Restrain selection
    .input_size = 400
  coordinate_sigma = 0.5
    .type = float
    .help = sigma value for coordinates restrained to starting positions
}
stop_for_unknowns = True
  .type = bool
  .short_caption = Stop for unknown residues
  .style = noauto
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
secondary_structure_restraints = False
  .type = bool
  .short_caption = Secondary structure restraints
secondary_structure
  .alias = refinement.secondary_structure
{
  include scope mmtbx.secondary_structure.sec_str_master_phil
}
use_c_beta_deviation_restraints=True
  .type = bool
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
  }
}
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
  if(cs is None):
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
    use_neutron_distances     = params.use_neutron_distances,
    mon_lib_srv               = monomer_library.server.server(),
    ener_lib                  = monomer_library.server.ener_lib())
  processed_pdb_file, junk = processed_pdb_files_srv.\
    process_pdb_files(pdb_file_names = pdb_file_names) # XXX remove junk
  return processed_pdb_file

def get_geometry_restraints_manager(processed_pdb_file, xray_structure, params,
    log=sys.stdout):
  has_hd = None
  if(xray_structure is not None):
    sctr_keys = xray_structure.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
  hbond_params = None
  if(params.secondary_structure_restraints):
    sec_str = mmtbx.secondary_structure.process_structure(
      params             = params.secondary_structure,
      processed_pdb_file = processed_pdb_file,
      tmp_dir            = os.getcwd(),
      log                = log,
      assume_hydrogens_all_missing=(not has_hd))
    sec_str.initialize(log=log)
    build_proxies = sec_str.create_hbond_proxies(
      log          = log,
      hbond_params = None)
    hbond_params = build_proxies.proxies
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    params_edits                 = params.geometry_restraints.edits,
    plain_pairs_radius           = 5,
    hydrogen_bond_proxies        = hbond_params,
    assume_hydrogens_all_missing = not has_hd)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)
  if(xray_structure is not None):
    restraints_manager.crystal_symmetry = xray_structure.crystal_symmetry()
  return restraints_manager

def run_minimization(
      sites_cart,
      selection,
      restraints_manager,
      pdb_hierarchy,
      params,
      cdl,
      log):
  o = mmtbx.refinement.geometry_minimization.run2(
    sites_cart                     = sites_cart,
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
    generic_restraints             = False,
    rmsd_bonds_termination_cutoff  = params.rmsd_bonds_termination_cutoff,
    rmsd_angles_termination_cutoff = params.rmsd_angles_termination_cutoff,
    alternate_nonbonded_off_on     = params.alternate_nonbonded_off_on,
    cdl=cdl,
    log                            = log)

class run(object):
  _pdb_suffix = "minimized"
  def __init__(self, args, log, use_directory_prefix=True):
    self.log                = log
    self.params             = None
    self.inputs             = None
    self.args               = args
    self.processed_pdb_file = None
    self.xray_structure     = None
    self.pdb_hierarchy      = None
    self.selection          = None
    self.restrain_selection = None
    self.grm                = None
    self.sites_cart         = None
    self.time_strings       = []
    self.total_time         = 0
    self.output_file_name   = None
    self.pdb_file_names     = []
    self.use_directory_prefix = use_directory_prefix
    self.__execute()

  def __execute (self) :
    #
    self.caller(func = self.initialize,     prefix="Initialization, inputs")
    self.caller(func = self.process_inputs, prefix="Processing inputs")
    self.caller(func = self.atom_selection, prefix="Atom selection")
    self.caller(func = self.get_restraints, prefix="Geometry Restraints")
    self.caller(func = self.minimization,   prefix="Minimization")
    self.caller(func = self.write_pdb_file, prefix="Write PDB file")
    self.caller(func = self.write_geo_file, prefix="Write GEO file")
    #
    self.show_times()

  def master_params (self) :
    return master_params()

  def caller(self, func, prefix):
    timer = user_plus_sys_time()
    func(prefix = prefix)
    t = timer.elapsed()
    self.total_time += t
    fmt = "  %s: %s"%(prefix, str("%8.3f"%t).strip())
    self.time_strings.append(fmt)

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
    self.xray_structure = self.processed_pdb_file.xray_structure()
    self.pdb_hierarchy = self.processed_pdb_file.all_chain_proxies.pdb_hierarchy

  def atom_selection(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.selection = mmtbx.utils.atom_selection(
      all_chain_proxies = self.processed_pdb_file.all_chain_proxies,
      string = self.params.selection)
    print >> self.log, "  selected %s atoms out of total %s"%(
      str(self.selection.count(True)),str(self.selection.size()))
    self.generate_restrain_selection()

  def generate_restrain_selection(self):
    if(self.params.reference_restraints.\
       restrain_starting_coord_selection is not None):
      self.restrain_selection = mmtbx.utils.atom_selection(
        all_chain_proxies = self.processed_pdb_file.all_chain_proxies,
        string =
          self.params.reference_restraints.restrain_starting_coord_selection)

  def addcbetar(self):
    mmtbx.torsion_restraints.utils.add_c_beta_restraints(
      geometry      = self.grm.geometry,
      pdb_hierarchy = self.pdb_hierarchy,
      log           = self.log)

  def get_restraints(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.grm = get_geometry_restraints_manager(
      processed_pdb_file = self.processed_pdb_file,
      xray_structure     = self.xray_structure,
      params             = self.params,
      log                = self.log)
    if(self.params.use_c_beta_deviation_restraints): self.addcbetar()
    # reference
    if(self.restrain_selection is not None):
      restrain_sites_cart = self.processed_pdb_file.all_chain_proxies.\
        sites_cart.deep_copy().select(self.restrain_selection)
      self.grm.geometry.generic_restraints_manager.reference_manager.\
        add_coordinate_restraints(
          sites_cart=restrain_sites_cart,
          selection=self.restrain_selection,
          sigma=self.params.reference_restraints.coordinate_sigma)
      # sanity check
      assert self.grm.geometry.generic_restraints_manager.flags.reference is True
      assert self.grm.geometry.generic_restraints_manager.reference_manager.\
        reference_coordinate_proxies is not None
      assert len(self.grm.geometry.generic_restraints_manager.reference_manager.
        reference_coordinate_proxies) == len(restrain_sites_cart)

  def minimization(self, prefix): # XXX USE alternate_nonbonded_off_on etc
    broadcast(m=prefix, log = self.log)
    self.sites_cart = self.xray_structure.sites_cart()
    run_minimization(sites_cart = self.sites_cart, selection = self.selection,
      restraints_manager = self.grm, params = self.params.minimization,
      pdb_hierarchy = self.pdb_hierarchy,
      cdl=self.params.pdb_interpretation.cdl,
      log = self.log)

  def write_pdb_file(self, prefix):
    broadcast(m=prefix, log = self.log)
    self.xray_structure.set_sites_cart(sites_cart = self.sites_cart)
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
    self.pdb_hierarchy.write_pdb_file(file_name = ofn, crystal_symmetry =
      self.xray_structure.crystal_symmetry())
    self.output_file_name = os.path.abspath(ofn)

  def write_geo_file(self, prefix):
    if(self.params.write_geo_file and self.grm is not None):
      broadcast(m=prefix, log = self.log)
      ofn = os.path.basename(self.output_file_name).replace(".pdb",".geo")
      directory = self.params.directory
      if (self.use_directory_prefix) and (directory is not None) :
        ofn = os.path.join(directory, ofn)
      f=file(ofn,"wb")
      print >> self.log, "  output file name:", ofn
      print >> f, "# Geometry restraints after refinement"
      print >> f
      xray_structure = self.xray_structure
      sites_cart = xray_structure.sites_cart()
      site_labels = xray_structure.scatterers().extract_labels()
      self.grm.geometry.show_sorted(
        sites_cart=sites_cart,
        site_labels=site_labels,
        f=f)
      f.close()

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
