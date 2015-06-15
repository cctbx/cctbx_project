from __future__ import division
import sys
import iotbx.pdb
import mmtbx.utils
from mmtbx import monomer_library
from scitbx.array_family import flex
from libtbx.utils import Sorry
import iotbx.phil
import mmtbx.refinement.real_space.expload_and_refine

master_phil = iotbx.phil.parse("""

  input_files {
    map_coeffs_file = None
      .type = path
      .help = File with map coefficients
      .short_caption = Map coefficients
      .style = bold file_type:hkl input_file process_hkl child:fobs:data_labels\
        child:space_group:space_group child:unit_cell:unit_cell

    map_coeffs_labels = None
      .type = str
      .input_size = 160
      .help = Optional label specifying which columns of of map coefficients \
          to use
      .short_caption = Map coeffs label
      .style = bold renderer:draw_fobs_label_widget

    map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file

    pdb_in = None
      .type = path
      .help = Input PDB file to minimize
      .short_caption = Input PDB file

  }
  output_files {

    pdb_out = minimize_ca.pdb
      .type = path
      .help = Output PDB file with CA positions
      .short_caption = Output PDB file

    prefix = tst_00
      .type = str
      .help = Prefix for output files
      .short_caption = Prefix for output files
  }
  crystal_info {
     resolution = None
       .type = float
       .help = High-resolution limit. Data will be truncated at this\
               resolution. If a map is supplied, it will be Fourier \
               filtered at this resolution. Required if input is a map and \
                only_original_map is not set.
       .short_caption = High-resolution limit
       .style = resolution
     space_group = None
       .type = space_group
       .short_caption = Space Group
       .help = Space group (normally read from the data file)
     unit_cell = None
       .type = unit_cell
       .short_caption = Unit Cell
       .help = Unit Cell (normally read from the data file)
  }
  sharpening {
     b_sharpen = None
       .type = float
       .help = Ignored. \
               B-factor for sharpening (positive is sharpen, negative is blur) \
               Note: if resolution is specified map will be Fourier filtered \
               to that resolution even if it is not sharpened.
       .short_caption = B-factor for sharpening
  }
  minimization {
     strategy = *ca_only all_atoms
       .type = choice
       .help = Ignored. \
          Strategy.  CA_only uses just CA atoms, all_atoms uses all
       .short_caption = CA only or all atoms

     number_of_macro_cycles=5
       .type = int
       .short_caption = Number of overall cycles of minimization
       .help = Number of overall (macro) cycles of minimization

     target_bond_rmsd = 0.02
       .type = float
       .short_caption = Target bond rmsd
       .help = Target bond rmsd

     target_angle_rmsd = 2.0
       .type = float
       .short_caption = Target angle RMSD
       .help = Target angle RMSD

     number_of_trials = 20
       .type = int
       .short_caption = Number of trials
       .help = Number of trials

     start_xyz_error = 10.0
       .type = float
       .short_caption = Starting coordinate error
       .help = Starting coordinate error

  }
  control {
      verbose = False
        .type = bool
        .help = Verbose output
        .short_caption = Verbose output
  }
""", process_includes=True)
master_params = master_phil

def get_params(args,out=sys.stdout):
  command_line = iotbx.phil.process_command_line_with_files(
    reflection_file_def="input_files.map_coeffs_file",
    map_file_def="input_files.map_file",
    pdb_file_def="input_files.pdb_in",
    args=args,
    master_phil=master_phil)
  params = command_line.work.extract()
  print >>out,"\nMinimize_ca ... optimize a coarse-grain model in EM or low-resolution X-ray map"
  master_phil.format(python_object=params).show(out=out)
  return params

def get_params_edits(n_ca=None,edit_file_name='bond_edits.eff'):

  from mmtbx.monomer_library.pdb_interpretation import geometry_restraints_edits_str
  edit_phil = iotbx.phil.parse(
  """refinement.geometry_restraints.edits
      { %s }
  """ %geometry_restraints_edits_str
  )
  # Generate restraint file
  f=open(edit_file_name,'w')
  print >>f,"""
  refinement.geometry_restraints.edits {
    excessive_bond_distance_limit = 10
  """

  for i in xrange(1,n_ca+1):
    if i+1 <=n_ca:
     print >>f,"""
      bond {

      action = *add delete change
      atom_selection_1 = chain E and resseq %d  and name CA
      atom_selection_2 = chain E and resseq %d and name CA
      distance_ideal = 3.8
      sigma = 0.5
      slack = None
      symmetry_operation = None
      }
    """ %(i,i+1)

    if i+2 <=n_ca:
     print >>f,"""
      bond {

      action = *add delete change
      atom_selection_1 = chain E and resseq %d  and name CA
      atom_selection_2 = chain E and resseq %d and name CA
      distance_ideal = 5.5
      sigma = 0.5
      slack = None
      symmetry_operation = None
      }
    """ %(i,i+2)

  print >>f,"""
    }
  """
  f.close()
  command_line = iotbx.phil.process_command_line_with_files(
      args=[edit_file_name],master_phil=edit_phil)
  params = command_line.work.extract()
  edit_phil.format(python_object=params).show()
  return params.refinement.geometry_restraints.edits

def ccp4_map(crystal_symmetry, file_name, map_data):
  from iotbx import ccp4_map
  ccp4_map.write_ccp4_map(
      file_name=file_name,
      unit_cell=crystal_symmetry.unit_cell(),
      space_group=crystal_symmetry.space_group(),
      #gridding_first=(0,0,0),# This causes a bug (map gets shifted)
      #gridding_last=n_real,  # This causes a bug (map gets shifted)
      map_data=map_data,
      labels=flex.std_string([""]))

def run(args,map_data=None,map_coeffs=None,pdb_inp=None,pdb_string=None,
    params_edits=None,out=sys.stdout):

  # Get the parameters
  params=get_params(args=args,out=out)

  # Get the starting model and reset the sequence if necessary
  if pdb_inp is None:
    if not pdb_string:
      if params.input_files.pdb_in:
        pdb_string=open(params.input_files.pdb_in).read()
      else:
        raise Sorry("Need an input PDB file")
    pdb_inp=iotbx.pdb.input(source_info=None, lines = pdb_string)
  hierarchy=pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  xrs=pdb_inp.xray_structure_simple()
  if not pdb_string:
    pdb_string=hierarchy.as_pdb_string()
  #
  asc = hierarchy.atom_selection_cache()
  sel_ca    = asc.selection(string="name CA")
  sel_other = ~sel_ca
  #

  # Initialize states accumulator
  states = mmtbx.utils.states(pdb_hierarchy=hierarchy, xray_structure=xrs)
  states.add(sites_cart = xrs.sites_cart())

  # Build geometry restraints
  pi_params = monomer_library.pdb_interpretation.master_params
  pi_params = pi_params.extract()
  pi_params.max_reasonable_bond_distance=500
  pi_params.nonbonded_weight=2000
  #pi_params.peptide_link.ramachandran_restraints=True
  #pi_params.peptide_link.rama_potential="emsley"
  if 0:
    pi_params.nonbonded_distance_cutoff=10
    pi_params.vdw_1_4_factor=2.0
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = monomer_library.server.server(),
    ener_lib                 = monomer_library.server.ener_lib(),
    params                   = pi_params,
    raw_records              = pdb_string,
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies                = False,
    plain_pairs_radius           = 5,
    params_edits = params_edits,
    assume_hydrogens_all_missing = True)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)

  if params.control.verbose:
    geometry.write_geo_file()

  ear = mmtbx.refinement.real_space.expload_and_refine.run(
    xray_structure          = xrs,
    pdb_hierarchy           = hierarchy,
    map_data                = map_data,
    restraints_manager      = restraints_manager,
    number_of_macro_cycles  = params.minimization.number_of_macro_cycles,
    target_bond_rmsd        = params.minimization.target_bond_rmsd,
    target_angle_rmsd       = params.minimization.target_angle_rmsd,
    number_of_trials        = params.minimization.number_of_trials,
    xyz_shake               = params.minimization.start_xyz_error,
    states                  = states)

  return ear.pdb_hierarchy, ear.xray_structure, ear.states

if (__name__ == "__main__"):
  args=sys.argv[1:]
  run(args=args)
