from __future__ import absolute_import, division, print_function
from iotbx import pdb
import iotbx.phil
import iotbx.ncs
from mmtbx.monomer_library import server
from mmtbx.monomer_library import cif_types
from mmtbx.monomer_library import rna_sugar_pucker_analysis
from mmtbx.monomer_library import conformation_dependent_restraints
from mmtbx.geometry_restraints import ramachandran
import mmtbx.geometry_restraints
import mmtbx.geometry_restraints.torsion_restraints.utils
# from mmtbx.secondary_structure.build import model_idealization_master_phil_str
from cctbx import geometry_restraints
import cctbx.geometry_restraints.manager
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.python_utils import dicts
from libtbx.str_utils import show_string
from libtbx.utils import flat_list, Sorry, user_plus_sys_time, plural_s
from libtbx.utils import format_exception
from libtbx import Auto, group_args, slots_getstate_setstate
from past.builtins import cmp
from six.moves import cStringIO as StringIO
import string
import sys, os
import time
import math

from cctbx.geometry_restraints.linking_class import linking_class
from six.moves import zip, range
origin_ids = linking_class()

# see iotbx/pdb/common_residue_names.h; additionally here only: U I S
ad_hoc_single_atom_residue_element_types = """\
ZN CA MG CL NA MN K FE CU CD HG NI CO BR XE SR CS PT BA TL PB SM AU RB YB LI
KR MO LU CR OS GD TB LA F AR AG HO GA CE W SE RU RE PR IR EU AL V TE SB PD
U I S
""".split()

class ad_hoc_single_atom_residue(object):

  def __init__(self, residue_name, atom_name, atom_element):
    atom_element = atom_element.strip().upper()
    if (atom_element in ad_hoc_single_atom_residue_element_types):
      self.scattering_type = atom_element
      self.energy_type = atom_element
      return
    atom_name = atom_name.strip().upper()
    if (    len(atom_element) == 0
        and atom_name == residue_name.strip()
        and atom_name in ad_hoc_single_atom_residue_element_types):
      self.scattering_type = atom_name
      self.energy_type = atom_name
      return
    if (    residue_name == "NH3"
        and atom_element == "N"
        or (len(atom_element) == 0 and atom_name.startswith("N"))):
      self.scattering_type = "N"
      self.energy_type = "N"
      return
    if (    residue_name == "CH4"
        and atom_element == "C"
        or (len(atom_element) == 0 and atom_name.startswith("C"))):
      self.scattering_type = "C"
      self.energy_type = "C"
      return
    self.scattering_type = None
    self.energy_type = None

dihedral_function_type_params_str = """\
  dihedral_function_type = *determined_by_sign_of_periodicity \
                            all_sinusoidal \
                            all_harmonic
    .type=choice
    .optional=False
"""

clash_guard_params_str = """\
  clash_guard
    .short_caption = Clash guard
    .style = noauto box auto_align
    .expert_level=2
  {
    nonbonded_distance_threshold = 0.5
      .type = float
    max_number_of_distances_below_threshold = 100
      .type = int
    max_fraction_of_distances_below_threshold = 0.1
      .type = float
  }
"""

test_cdl_params = """\
  conformation_dependent_restraints_ideal = True
    .type = bool
  conformation_dependent_restraints_esd = True
    .type = bool
  conformation_dependent_restraints_default = True
    .type = bool
  cdl_interpolation = False
    .type = bool
  cdl_weight = 1.
    .type = float
"""

altloc_weighting_params = """\
  altloc_weighting
    .expert_level = 2
    .style = box noauto auto_align
  {
    weight = False
      .type = bool
    bonds = True
      .type = bool
    angles = True
      .type = bool
    factor = 1.
      .type = float
    sqrt = False
      .type = bool
    min_occupancy = 0.5
      .type = float
  }
"""
restraints_library_str = """
  restraints_library
    .short_caption = Restraints library selection
    .style = box auto_align
  {
    cdl = True
      .type = bool
      .short_caption = Use Conformation-Dependent Library
      .help = Use Conformation Dependent Library (CDL) \
        for geometry restraints
      .style = bold
    mcl = True
      .type = bool
      .short_caption = Use Metal Coordination Library (MCL)
      .help = Use Metal Coordination Library (MCL) \
        for tetrahedral Zn++ and iron-sulfur clusters SF4, FES, F3S, ...
      .style = bold
    cis_pro_eh99 = False
      .type = bool
      .style = hidden
    omega_cdl = False
      .type = bool
      .short_caption = Use Omega Conformation-Dependent Library
      .help = Use Omega Conformation Dependent Library (omega-CDL) \
        for geometry restraints
      .style = hidden
    cdl_svl = False
      .type = bool
      .short_caption = Use improved SVL values for CDL classes
      .style = hidden
    rdl = False
      .type = bool
      .style = hidden
    hpdl = False
      .type = bool
      .style = hidden
  }
"""
ideal_ligands = ['SF4', 'F3S', 'DVT']
ideal_ligands_str = ' '.join(ideal_ligands)
master_params_str = """\
  %(restraints_library_str)s
  sort_atoms = True
    .type = bool
    .short_caption = Sort atoms in input pdb so they would be in the same order
  superpose_ideal_ligand = *None all %(ideal_ligands_str)s
    .type = choice(multi=True)
    .short_caption = Substitute correctly oriented SF4 metal cluster
  flip_symmetric_amino_acids = False
    .type = bool
    .short_caption = Flip symmetric amino acids to conform to IUPAC convention
    .style = noauto
  disable_uc_volume_vs_n_atoms_check = False
    .type = bool
    .short_caption = Disable check of unit cell volume to be compatible with the \
                     number of atoms
  allow_polymer_cross_special_position = False
    .type = bool
  correct_hydrogens = True
    .type = bool
    .short_caption = Correct the hydrogen positions trapped in chirals etc
  # include scope mmtbx.secondary_structure.build.model_idealization_master_phil_str
  include scope mmtbx.secondary_structure.sec_str_master_phil
  c_beta_restraints=True
    .type = bool
    .short_caption = Use C-beta deviation restraints
  reference_coordinate_restraints
    .short_caption = Reference coordinate restraints
    .caption = Harmonic restraints on the starting coordinates
    .help = Restrains coordinates in Cartesian space to stay near their \
      starting positions.  This is intended for use in generating \
      simulated annealing omit maps, to prevent refined atoms from \
      collapsing in on the region missing atoms.  For conserving geometry \
      quality at low resolution, the more flexible reference model \
      restraints should be used instead.
    .style = box auto_align
  {
    enabled = False
      .type = bool
    exclude_outliers = True
      .caption = Exclude ramachandran plot and rotamer outliers from reference \
        coordinate selection
      .type = bool
    selection = all
      .type = atom_selection
    sigma = 0.2
      .type = float
    limit = 1.0
      .type = float
    top_out = False
      .type = bool
  }
  automatic_linking
    .style = box auto_align noauto
    .short_caption = Automatic covalent linking
    .caption = You may choose to run with the defaults or use the buttons \
      below to select all or none of the linking options. Also, you may enable \
      individual link types separately and adjust the cutoff values to \
      fine-tune the link candidates.
  {
    link_all = False
      .type = bool
      .short_caption = Automatic linking
      .help = If True, bond restraints will be generated for any appropriate \
        ligand-protein or ligand-nucleic acid covalent bonds. This includes \
        sugars, amino acid modifications, and other prosthetic groups.
      .style = noauto #renderer:draw_automatic_linking_control
    link_none = False
      .type = bool
      .short_caption = Overrides all automatic linking
      .style = noauto
    link_metals = False
      .type = bool
    link_residues = False
      .type = bool
      .short_caption = Link amino acids in "special" ways such as cyclic and \
        side-chain to side-chain links
    link_amino_acid_rna_dna = False
      .type = bool
      .short_caption = Link any RNA/DNA to amino acid possibities
    link_carbohydrates = True
      .type = bool
      .short_caption = Link carbohydrates to protein and other carbohydrates
    link_ligands = True
      .type = bool
      .short_caption = Link ligands to protein
    link_small_molecules = False
      .type = bool
      .short_caption = Link small molecules such as SO4, PO4 to protein
    metal_coordination_cutoff = 3.5
      .type = float
      .short_caption = Maximum distance for automatic linking of metals
    amino_acid_bond_cutoff = 1.9
      .type = float
      .short_caption = Distance cutoff for automatic linking of aminoacids
    inter_residue_bond_cutoff = 2.2
      .type = float
      .short_caption = Distance cutoff for automatic linking of other residues
    buffer_for_second_row_elements = 0.5
      .type = float
      .short_caption = Distance to add to intra_residue_bond_cutoff if one \
        or more element is in the second (or more) row
    carbohydrate_bond_cutoff = 1.99
      .type = float
    ligand_bond_cutoff = 1.99
      .type = float
    small_molecule_bond_cutoff = 1.98
      .type = float
  }
  include_in_automatic_linking
    .optional = True
    .multiple = True
    .short_caption = exclude
    .style = noauto auto_align
  {
    selection_1 = None
      .type = atom_selection
    selection_2 = None
      .type = atom_selection
    bond_cutoff = 4.5
      .type = float
  }
  exclude_from_automatic_linking
    .optional = True
    .multiple = True
    .short_caption = exclude
    .style = noauto auto_align
  {
    selection_1 = None
      .type = atom_selection
    selection_2 = None
      .type = atom_selection
  }
  use_neutron_distances = False
    .type = bool
    .short_caption = Use the nuclear distances for X-H/D
    .help = Use neutron X-H distances (which are longer than X-ray ones)
  apply_cis_trans_specification
    .optional = True
    .multiple = True
    .short_caption = Modify default cis-trans specification
    .style = noauto auto_align
  {
    cis_trans_mod = cis *trans
      .type = choice
    residue_selection = None
      .type = atom_selection
      .help = Residues containing C-alpha atom of omega dihedral
  }
  apply_cif_restraints
    .optional = True
    .multiple = True
    .short_caption = Use CIF restraints
    .style = noauto auto_align
  {
    restraints_file_name= None
      .type = path
    residue_selection = None
      .type = atom_selection
  }
  apply_cif_modification
    .optional = True
    .multiple = True
    .short_caption = Modify CIF
    .style = noauto auto_align
  {
    data_mod = None
      .type = str
    residue_selection = None
      .type = atom_selection
  }
  apply_cif_link
    .optional = True
    .multiple = True
    .short_caption = Add link to CIF
    .style = noauto auto_align
  {
    data_link = None
      .type = str
    residue_selection_1 = None
      .type = atom_selection
    residue_selection_2 = None
      .type = atom_selection
  }
  disulfide_bond_exclusions_selection_string = None
    .type = str
  exclusion_distance_cutoff = 3
    .type = float
    .help = If SG of CYS forming SS bond is closer than this distance to an \
            atom that it may coordinate then this SG is excluded from SS bond.
  link_distance_cutoff = 3
    .type=float
    .optional=False
    .help = Length of link between the linked residues
  disulfide_distance_cutoff = 3
    .type=float
    .optional=False
  add_angle_and_dihedral_restraints_for_disulfides = True
    .type = bool
    .optional = False
  %(dihedral_function_type_params_str)s
  chir_volume_esd = 0.2
    .type=float
    .optional=False
    .short_caption = Chiral volume E.S.D.
  peptide_link
    .short_caption = Peptide link settings
    .style = box
  {
    ramachandran_restraints = False
      .type = bool
      .help = !!! OBSOLETED. Kept for backward compatibility only !!! \
        Restrains peptide backbone to fall within allowed regions of \
        Ramachandran plot.  Although it does not eliminate outliers, it can \
        significantly improve the percent favored and percent outliers at \
        low resolution.  Probably not useful (and maybe even harmful) at \
        resolutions much higher than 3.5A.
      .expert_level = 2
      .short_caption = Ramachandran restraints
    cis_threshold = 45
      .type = float
      .optional = False
      .short_caption = Threshold (degrees) for cis-peptides
    apply_all_trans = False
      .type = bool
      .style = hidden
    discard_omega = False
      .type = bool
    discard_psi_phi = True
      .type = bool
      .optional = False
      .short_caption = Ignore monomer library Phi/Psi restraints
    apply_peptide_plane = False
      .type = bool
    omega_esd_override_value = None
      .type = float
      .short_caption = Omega-ESD override value
    include scope mmtbx.geometry_restraints.ramachandran.old_master_phil
  }
  include scope mmtbx.geometry_restraints.ramachandran.master_phil
  max_reasonable_bond_distance = 50.0
    .type=float
  nonbonded_distance_cutoff = None
    .type=float
  default_vdw_distance = 1
    .type=float
    .optional=False
    .short_caption = Default VDW distance
  min_vdw_distance = 1
    .type=float
    .optional=False
    .short_caption = Minimum VDW distance
  nonbonded_buffer = 1
    .type=float
    .optional=False
    .help = **EXPERIMENTAL, developers only**
    .expert_level = 3
  nonbonded_weight = None
    .type = float
    .short_caption = Nonbonded weight
    .help = Weighting of nonbonded restraints term.  By default, this will be \
      set to 16 if explicit hydrogens are used (this was the default in \
      earlier versions of Phenix), or 100 if hydrogens are missing.
  const_shrink_donor_acceptor = 0.6
    .type=float
    .optional=False
    .expert_level=3
    .help = **EXPERIMENTAL, developers only**
  vdw_1_4_factor = 0.8
    .type=float
    .optional=False
    .short_caption = Scale factor for 1-4 VDW interactions
  min_distance_sym_equiv = 0.5
    .type=float
    .short_caption = Min. symmetry-equivalent distance
  custom_nonbonded_symmetry_exclusions = None
    .optional = True
    .type = atom_selection
    .multiple = True
  translate_cns_dna_rna_residue_names = None
    .type=bool
    .optional=False
    .short_caption = Translate CNS DNA/RNA names
  proceed_with_excessive_length_bonds = False
    .type=bool
  rna_sugar_pucker_analysis
    .short_caption = RNA sugar pucker analysis
    .style = box noauto auto_align menu_item parent_submenu:advanced
  {
    include scope mmtbx.monomer_library.rna_sugar_pucker_analysis.master_phil
  }
  show_histogram_slots
    .style = box auto_align noauto
    .expert_level = 2
  {
    bond_lengths = 5
      .type=int
    nonbonded_interaction_distances = 5
      .type=int
    bond_angle_deviations_from_ideal = 5
      .type=int
    dihedral_angle_deviations_from_ideal = 5
      .type=int
    chiral_volume_deviations_from_ideal = 5
      .type=int
  }
  show_max_items
    .expert_level = 2
    .style = box auto_align noauto
  {
    not_linked = 5
      .type=int
    bond_restraints_sorted_by_residual = 5
      .type=int
    nonbonded_interactions_sorted_by_model_distance = 5
      .type=int
    bond_angle_restraints_sorted_by_residual = 5
      .type=int
    dihedral_angle_restraints_sorted_by_residual = 3
      .type=int
    chirality_restraints_sorted_by_residual = 3
      .type=int
    planarity_restraints_sorted_by_residual = 3
      .type=int
    residues_with_excluded_nonbonded_symmetry_interactions = 12
      .type=int
    fatal_problem_max_lines = 10
      .type=int
  }
  include scope iotbx.ncs.ncs_group_phil_str
  include scope iotbx.ncs.ncs_search_options
  %(clash_guard_params_str)s
""" % vars()

master_params = iotbx.phil.parse(
  input_string=master_params_str,
  process_includes=True)

geometry_restraints_edits_str = """\
excessive_bond_distance_limit = 10
  .type = float
bond
  .optional = True
  .multiple = True
  .short_caption = Bond
  .style = auto_align
{
  action = *add delete change
    .type = choice
  atom_selection_1 = None
    .type = atom_selection
    .input_size = 400
  atom_selection_2 = None
    .type = atom_selection
    .input_size = 400
  symmetry_operation = None
    .help = "The bond is between atom_1 and symmetry_operation * atom_2,"
            " with atom_1 and atom_2 given in fractional coordinates."
            " Example: symmetry_operation = -x-1,-y,z"
    .type = str
  distance_ideal = None
    .type = float
  sigma = None
    .type = float
  slack = None
    .type = float
  limit = -1.0
    .type = float
  top_out = False
    .type = bool
}
angle
  .optional = True
  .multiple = True
  .short_caption = Angle
  .style = auto_align
{
  action = *add delete change
    .type = choice
  atom_selection_1 = None
    .type = atom_selection
    .input_size = 400
  atom_selection_2 = None
    .type = atom_selection
    .input_size = 400
  atom_selection_3 = None
    .type = atom_selection
    .input_size = 400
  angle_ideal = None
    .type = float
  sigma = None
    .type = float
}
dihedral
  .optional = True
  .multiple = True
  .short_caption = Dihedral
  .style = hidden
{
  action = *add delete change
    .type = choice
  atom_selection_1 = None
    .type = atom_selection
    .input_size = 400
  atom_selection_2 = None
    .type = atom_selection
    .input_size = 400
  atom_selection_3 = None
    .type = atom_selection
    .input_size = 400
  atom_selection_4 = None
    .type = atom_selection
    .input_size = 400
  angle_ideal = None
    .type = float
  alt_angle_ideals = None
    .type = floats
  sigma = None
    .type = float
  periodicity = 1
    .type = int
}
planarity
  .optional = True
  .multiple = True
  .short_caption = Planarity
  .style = auto_align
{
  action = *add delete change
    .type = choice
  atom_selection = None
    .type = atom_selection
    .input_size = 400
  sigma = None
    .type = float
}
parallelity
  .optional = True
  .multiple = True
  .short_caption = Parallelity
  .style = auto_align
{
  action = *add delete change
    .type = choice
  atom_selection_1 = None
    .type = atom_selection
    .input_size = 400
  atom_selection_2 = None
    .type = atom_selection
    .input_size = 400
  sigma = 0.027
    .type = float
  target_angle_deg = 0
    .type = float
}
"""
obsoleted_scale_restraints = """
scale_restraints
  .multiple = True
  .optional = True
  .help = Apply a scale factor to restraints for specific atom selections, \
    to tighten geometry without changing the overall scale of the geometry \
    target.
{
  atom_selection = None
    .type = atom_selection
  scale = 1.0
    .type = float
  apply_to = *bond *angle *dihedral *chirality
    .type = choice(multi=True)
}
"""

def validate_geometry_edits_params(params):
  for k, bond in enumerate(params.bond):
    if (None in [bond.atom_selection_1, bond.atom_selection_2]):
      raise Sorry(("A custom bond definition (#%d in the list) is incomplete; "+
        "two atom selections are required.") % k)
    elif (bond.distance_ideal is None):
      raise Sorry(("The ideal distance for custom bond #%d is not defined. "+
        "(Atom selections: %s, %s)") % (k, bond.atom_selection_1,
        bond.atom_selection_2))
    elif (bond.sigma is None):
      raise Sorry(("The sigma for custom bond #%d is not defined. "+
        "(Atom selections: %s, %s)") % (k, bond.atom_selection_1,
        bond.atom_selection_2))
  for k, angle in enumerate(params.angle):
    if (None in [angle.atom_selection_1, angle.atom_selection_2,
                 angle.atom_selection_3]):
      raise Sorry(("A custom angle definition (#%d in the list) is "+
        "incomplete; two atom selections are required.") % k)
    elif (angle.angle_ideal is None):
      raise Sorry(("The ideal angle for custom angle #%d is not defined. "+
        "(Atom selections: %s, %s, %s)") % (k, angle.atom_selection_1,
        angle.atom_selection_2, angle.atom_selection_3))
    elif (angle.sigma is None):
      raise Sorry(("The sigma for custom angle #%d is not defined. "+
        "(Atom selections: %s, %s, %s)") % (k, angle.atom_selection_1,
        angle.atom_selection_2, angle.atom_selection_3))
  for k, plane in enumerate(params.planarity):
    if (plane.atom_selection is None):
      raise Sorry("The atom selection for custom plane #%d is not defined."% k)
    elif (plane.sigma is None):
      raise Sorry(("The sigma for custom plane #%d is not defined.  (Atom "+
        "selection: '%s')") % (k, plane.atom_selection))
  for k, parallelity in enumerate(params.parallelity):
    if (None in [parallelity.atom_selection_1, parallelity.atom_selection_2]):
      raise Sorry(("A custom parallelity definition (#%d in the list) is "+
          "incomplete; two atom selections are required.") % k)
  return True

geometry_restraints_remove_str = """\
angles=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
dihedrals=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
chiralities=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
planarities=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
parallelities=None
  .optional=True
  .type=str
  .multiple=True
  .input_size=400
"""

grand_master_phil_str = """\
pdb_interpretation
  .alias = refinement.pdb_interpretation
  .short_caption = Model interpretation
{
  %(master_params_str)s
}
geometry_restraints
  .alias = refinement.geometry_restraints
{
  edits
    .short_caption = Custom geometry restraints
  {
    %(geometry_restraints_edits_str)s
  }
  remove {
    %(geometry_restraints_remove_str)s
  }
}
""" % vars()

def flush_log(log):
  if (log is not None):
    flush = getattr(log, "flush", None)
    if (flush is not None): flush()

def all_atoms_are_in_main_conf(atoms):
  for atom in atoms:
    if (atom.parent().altloc != ""): return False
  return True

def residue_id_str(residue, suppress_segid=0):
  try :
    return residue.id_str(suppress_segid=suppress_segid)
  except ValueError as e :
    raise Sorry(str(e))

class counters(object):

  def __init__(self, label):
    self.label = label
    self.corrupt_monomer_library_definitions = 0
    self.already_assigned_to_first_conformer = 0
    self.unresolved_non_hydrogen = 0
    self.unresolved_hydrogen = 0
    self.undefined = 0
    self.resolved = 0
    self.discarded_because_of_special_positions = 0

class special_position_dict():
  def __init__(self, special_position_indices):
    self.spi = special_position_indices
    self.iseq_mapping = {}

  def involves_special_positions(self, i_seqs):
    if self.spi is None: return False
    for i_seq in i_seqs:
      mapping = self.iseq_mapping.get(i_seq)
      if mapping:
        return True
      elif mapping is None:
        self.iseq_mapping[i_seq] = i_seq in self.spi
    return False

def involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs):
  if (broken_bond_i_seq_pairs is None): return False
  i_seqs = sorted(i_seqs)
  for i in range(len(i_seqs)-1):
    for j in range(i+1,len(i_seqs)):
      if ((i_seqs[i],i_seqs[j]) in broken_bond_i_seq_pairs):
        return True
  return False

class source_info_server(object):

  def __init__(self, m_i, m_j):
    self.m_i, self.m_j = m_i, m_j

  def labels(self):
    if (self.m_j is None):
      return "residue: %s" % self.m_i.residue_altloc()
    return "residues: %s + %s" % (
      self.m_i.residue_altloc(),
      self.m_j.residue_altloc())

  def n_expected_atoms(self):
    if (self.m_j is None):
      return len(self.m_i.expected_atoms)
    return len(self.m_i.expected_atoms) \
         + len(self.m_j.expected_atoms)

def _show_atom_labels(pdb_atoms, i_seqs, out=None, prefix="", max_lines=None):
  if (out is None): out = sys.stdout
  for i_line,i_seq in enumerate(i_seqs):
    if (i_line == max_lines and len(i_seqs) > max_lines+1):
      print(prefix + "... (remaining %d not shown)" % (
        len(i_seqs)-max_lines), file=out)
      break
    print(prefix + pdb_atoms[i_seq].quote(), file=out)

def format_exception_message(
      m_i,
      m_j,
      i_seqs,
      base_message,
      source_labels=None,
      show_residue_names=True,
      lines=[]):
  s = StringIO()
  print(base_message, file=s)
  for line in lines:
    print(" ", line, file=s)
  if (source_labels is not None):
    for i,label in enumerate(source_labels):
      print("  %d. definition from: %s" % (i+1, label), file=s)
  if (show_residue_names):
    print("  " + source_info_server(m_i, m_j).labels(), file=s)
  print("  atom%s:" % plural_s(len(i_seqs))[1], file=s)
  _show_atom_labels(
    pdb_atoms=m_i.pdb_atoms, i_seqs=i_seqs, out=s, prefix="    ", max_lines=10)
  return s.getvalue()[:-1]

def discard_conflicting_pdb_element_column(mm):
  trusted_library_definitions = [
    "GLY", "VAL", "ALA", "LEU", "ILE", "PRO", "MET", "PHE", "TRP", "SER",
    "THR", "TYR", "CYS", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
    "MSE"]
  if (mm.residue_name in trusted_library_definitions): return True
  if (mm.residue_name[3:] == "%COO"
      and mm.residue_name[:3] in trusted_library_definitions):
    return True
  return False

class type_symbol_registry_base(object):

  def __init__(self, type_label, symbols, strict_conflict_handling):
    assert type_label in ["scattering", "nonbonded energy"]
    self.type_label = type_label
    self.symbols = symbols
    self.strict_conflict_handling = strict_conflict_handling
    self.n_resolved_conflicts = 0
    self.source_labels = flex.std_string(symbols.size())
    self.source_n_expected_atoms = flex.int(symbols.size(), -1)
    self.charges = flex.int(symbols.size(), 0)

  def discard_tables(self):
    self.source_labels = None
    self.source_n_expected_atoms = None

  def assign_directly(self, i_seq, symbol):
    self.symbols[i_seq] = symbol

  def assign_charge(self, i_seq, charge=0):
    self.charges[i_seq] = charge

  def assign_from_monomer_mapping(self, conf_altloc, mm):
    atom_dict = mm.monomer_atom_dict
    for atom_id,atom in mm.expected_atoms.items():
      i_seq = atom.i_seq
      if (self.type_label == "scattering"):
        symbol = atom_dict[atom_id].type_symbol
        if (symbol == "H" and self.symbols[i_seq] == "D"):
          symbol = "D"
      else:
        symbol = atom_dict[atom_id].type_energy
      if (symbol is None): continue
      source_label = mm.residue_altloc()
      source_n_expected_atoms = len(mm.expected_atoms)
      prev_symbol = self.symbols[i_seq]
      prev_source_label = self.source_labels[i_seq]
      prev_source_n_expected_atoms = self.source_n_expected_atoms[i_seq]
      assign = False
      raise_conflict = False
      if (prev_symbol == ""
          or (self.type_label == "scattering") and prev_symbol in ["Q","X"]):
        assign = True
      elif (prev_symbol.upper() != symbol.upper()):
        if (self.strict_conflict_handling):
          raise_conflict = True
        elif (self.type_label == "nonbonded energy"):
          if (prev_source_n_expected_atoms == source_n_expected_atoms):
            raise_conflict = True
          else:
            self.n_resolved_conflicts += 1
            if (prev_source_n_expected_atoms < source_n_expected_atoms):
              assign = True
        elif (self.type_label == "scattering"):
          if (prev_source_label == ""
              and discard_conflicting_pdb_element_column(mm)):
            self.n_resolved_conflicts += 1
            assign = True
          else:
            raise_conflict = True
      assert not (assign and raise_conflict)
      if (assign):
        self.symbols[i_seq] = symbol
        self.source_labels[i_seq] = source_label
        self.source_n_expected_atoms[i_seq] = source_n_expected_atoms
        if (self.type_label == "nonbonded energy"):
          charge_str = atom.charge_tidy(strip=True)
          if (charge_str is not None) and (len(charge_str) != 0):
            charge = 0
            if ("-" in charge_str):
              charge = - int(charge_str.replace("-", ""))
            elif ("+" in charge_str):
              charge = int(charge_str.replace("+", ""))
            self.charges[i_seq] = charge
          elif ((len(mm.expected_atoms) == 1) and
                (mm.residue_name.strip() == atom.name.strip()) and
                (atom.name.strip() == atom_dict[atom_id].type_symbol.strip())):
            # this is an ion, even if the charge isn't set in the PDB file -
            # the charge is set to an unrealistic number to avoid ambiguity
            self.charges[i_seq] = 99
      elif (raise_conflict):
        source = "with residue name %s" % source_label
        if (prev_source_label == ""):
          assert self.type_label == "scattering"
          prev_source = "from pdb element column"
        else:
          if (prev_source_label == source_label):
            source = "also " + source
          prev_source = "with residue name %s" % prev_source_label
        raise Sorry(format_exception_message(
          m_i=mm,
          m_j=None,
          i_seqs=[i_seq],
          base_message="Conflicting %s type symbols:" % self.type_label,
          show_residue_names=False,
          lines=['initial symbol: "%s" (%s)' % (prev_symbol, prev_source),
                 '    new symbol: "%s" (%s)' % (symbol, source)]))

  def n_unknown_type_symbols(self):
    return self.symbols.count("")

  def report_unknown_message(self):
    return "Number of atoms with unknown %s type symbols" % self.type_label

  def report(self, pdb_atoms, log, prefix, max_lines=10):
    n_unknown = self.n_unknown_type_symbols()
    if (n_unknown > 0):
      print("%s%s: %d" % (
        prefix, self.report_unknown_message(), n_unknown), file=log)
      i_seqs = (self.symbols == "").iselection()
      _show_atom_labels(
        pdb_atoms=pdb_atoms, i_seqs=i_seqs,
        out=log, prefix=prefix+"  ", max_lines=max_lines)
    if (self.n_resolved_conflicts > 0):
      print("%sNumber of resolved %s type symbol conflicts: %d" % (
        prefix, self.type_label, self.n_resolved_conflicts), file=log)

  def get_unknown_atoms(self, pdb_atoms, return_iseqs=False):
    n_unknown = self.n_unknown_type_symbols()
    if (n_unknown > 0):
      i_seqs = (self.symbols == "").iselection()
      if return_iseqs: return i_seqs
      return pdb_atoms.select(i_seqs)

class scattering_type_registry(type_symbol_registry_base):

  def __init__(self, scattering_types, strict_conflict_handling):
    type_symbol_registry_base.__init__(self,
      type_label="scattering",
      symbols=scattering_types,
      strict_conflict_handling=strict_conflict_handling)

class nonbonded_energy_type_registry(type_symbol_registry_base):

  def __init__(self, n_seq, strict_conflict_handling):
    type_symbol_registry_base.__init__(self,
      type_label="nonbonded energy",
      symbols=flex.std_string(n_seq),
      strict_conflict_handling=strict_conflict_handling)

class monomer_mapping_summary(slots_getstate_setstate):

  __slots__ = [
    "conf_altloc",
    "residue_name",
    "expected_atoms",
    "unexpected_atoms",
    "duplicate_atoms",
    "ignored_atoms",
    "classification",
    "incomplete_info",
    "is_terminus",
    "is_unusual"]

  def __init__(self, **keyword_args):
    for key in monomer_mapping_summary.__slots__:
      setattr(self, key, keyword_args.get(key, None))

  def all_associated_i_seqs(self):
    return flex.size_t(
        [a.i_seq for a in self.expected_atoms]
      + [a.i_seq for a in self.unexpected_atoms]
      + [a.i_seq for a in self.duplicate_atoms]
      + [a.i_seq for a in self.ignored_atoms])

  def summary(self):
    return self

  def __getstate__(self):
    sdict = super(monomer_mapping_summary, self).__getstate__()
    # The atom arrays have no pickling implemented so convert them to lists of
    # indices of the atoms in the common root hierarchy. Since a hierarchy can
    # be pickled we pass this on together with the sdict that holds the indices arrays.
    import iotbx_pdb_hierarchy_ext
    hroot = None
    if len(self.expected_atoms) and type(self.expected_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
      hroot = self.expected_atoms[0].parent().parent().parent().parent().parent()
    if len(self.unexpected_atoms) and type(self.unexpected_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
      hroot = self.unexpected_atoms[0].parent().parent().parent().parent().parent()
    if len(self.duplicate_atoms) and type(self.duplicate_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
      hroot = self.duplicate_atoms[0].parent().parent().parent().parent().parent()
    if len(self.ignored_atoms) and type(self.ignored_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
      hroot = self.ignored_atoms[0].parent().parent().parent().parent().parent()
    if hroot:
      # assuming atoms all come from the same hierarchy
      indexer = dict( ( a, i) for ( i, a ) in enumerate( hroot.atoms() ) )
      if len(self.expected_atoms) and type(self.expected_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
        sdict["expected_atoms"] = [indexer[ a ] for a in self.expected_atoms]
      if len(self.unexpected_atoms) and type(self.unexpected_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
        sdict["unexpected_atoms"] = [indexer[ a ]  for a in self.unexpected_atoms]
      if len(self.duplicate_atoms) and type(self.duplicate_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
        sdict["duplicate_atoms"] = [indexer[ a ]  for a in self.duplicate_atoms]
      if len(self.ignored_atoms) and type(self.ignored_atoms[0]) == iotbx_pdb_hierarchy_ext.atom:
        sdict["ignored_atoms"] = [indexer[ a ]  for a in self.ignored_atoms]
    return sdict, hroot

  def __setstate__(self, state):
    # Restore the atoms by using the indices arrays to pick the individual atoms
    # from their common root hierarchy
    sdict, hroot = state
    sdict["expected_atoms"] = []
    sdict["unexpected_atoms"] = []
    sdict["duplicate_atoms"] = []
    sdict["ignored_atoms"] = []
    if hroot:
      sdict["expected_atoms"] = [ hroot.atoms()[i] for i in sdict["expected_atoms"] ]
      sdict["unexpected_atoms"] = [ hroot.atoms()[i] for i in sdict["unexpected_atoms"] ]
      sdict["duplicate_atoms"] = [ hroot.atoms()[i] for i in sdict["duplicate_atoms"] ]
      sdict["ignored_atoms"] = [ hroot.atoms()[i] for i in sdict["ignored_atoms"] ]
    for name,value in sdict.items(): setattr(self, name, value)


class monomer_mapping(slots_getstate_setstate):

  __slots__ = [
    "active_atoms",
    "angle_counters",
    "atom_name_interpretation",
    "atom_names_given",
    "bond_counters",
    "chem_mod_ids",
    "chirality_counters",
    "classification",
    "conf_altloc",
    "dihedral_counters",
    "duplicate_atoms",
    "expected_atoms",
    "i_conformer",
    "ignored_atoms",
    "incomplete_info",
    "is_first_conformer_in_chain",
    "is_rna2p",
    "is_rna_dna",
    "is_terminus",
    "is_unusual",
    "lib_link",
    "missing_hydrogen_atoms",
    "missing_non_hydrogen_atoms",
    "mon_lib_names",
    "mon_lib_srv",
    "monomer",
    "monomer_atom_dict",
    "pdb_atoms",
    "pdb_residue",
    "pdb_residue_id_str",
    "planarity_counters",
    "parallelity_counters",
    "residue_name",
    "unexpected_atoms",
    "chainid",
    ]

  def __init__(self,
        pdb_atoms,
        mon_lib_srv,
        translate_cns_dna_rna_residue_names,
        rna_sugar_pucker_analysis_params,
        apply_cif_modifications,
        apply_cif_links_mm_pdbres_dict,
        i_model,
        i_conformer,
        is_first_conformer_in_chain,
        conf_altloc,
        pdb_residue,
        next_pdb_residue,
        chainid,
        specific_residue_restraints=None):
    self.chainid = chainid
    self.pdb_atoms = pdb_atoms
    self.mon_lib_srv = mon_lib_srv
    self.i_conformer = i_conformer
    self.is_first_conformer_in_chain = is_first_conformer_in_chain
    self.conf_altloc = conf_altloc
    self.pdb_residue = pdb_residue
    self.pdb_residue_id_str = residue_id_str(pdb_residue, suppress_segid=-1)
    self.residue_name = pdb_residue.resname
    atom_id_str_pdbres_list = self._collect_atom_names()
    self.monomer, self.atom_name_interpretation \
      = self.mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(
          residue_name=self.residue_name,
          atom_names=self.atom_names_given,
          translate_cns_dna_rna_residue_names
            =translate_cns_dna_rna_residue_names,
          specific_residue_restraints=specific_residue_restraints)
    if (self.atom_name_interpretation is None):
      self.mon_lib_names = None
    else:
      self.mon_lib_names = self.atom_name_interpretation.mon_lib_names()
    if (self.monomer is None):
      self.expected_atoms = {}
      self.unexpected_atoms = {}
      self.duplicate_atoms = {}
      for atom,atom_name in zip(self.active_atoms, self.atom_names_given):
        self.unexpected_atoms[atom_name] = atom
      self.incomplete_info = None
      self.is_terminus = None
    else:
      self.chem_mod_ids = set()
      self._rna_sugar_pucker_analysis(
        params=rna_sugar_pucker_analysis_params,
        next_pdb_residue=next_pdb_residue)
      self._get_mappings()
      for id_str in atom_id_str_pdbres_list:
        for apply_data_mod in apply_cif_modifications.get(id_str, []):
          self.apply_mod(
            mod_mod_id=self.mon_lib_srv.mod_mod_id_dict[apply_data_mod])
      self._set_incomplete_info()
      self.is_terminus = None
      self.monomer.set_classification()
      if (self.atom_name_interpretation is not None):
        d_aa_rn = getattr(
          self.atom_name_interpretation, "d_aa_residue_name", None)
        if (d_aa_rn is not None):
          mmid = self.mon_lib_srv.mod_mod_id_dict
          self.apply_mod(mod_mod_id=mmid["PEPT-D"])
          if (d_aa_rn == "DIL"):
            self.apply_mod(mod_mod_id=mmid["DIL_chir_02_both"])
          elif (d_aa_rn == "DTH"):
            self.apply_mod(mod_mod_id=mmid["DTH_chir_02_both"])
      cif_object = self.monomer.cif_object
      if (self.incomplete_info is None):
        self.resolve_unexpected()
      # strange behaviour of mods destroys cif_object
      self.monomer.cif_object = cif_object
    if (self.pdb_residue_id_str in apply_cif_links_mm_pdbres_dict):
      apply_cif_links_mm_pdbres_dict[self.pdb_residue_id_str].setdefault(
        self.i_conformer, []).append(self)

  def __setattr1__(self, attr, value):
    if attr == "lib_link":
      if value:
        for plane in value.plane_list:
          print(dir(plane))
          plane.show()
    slots_getstate_setstate.__setattr__(self, attr, value)

  def _collect_atom_names(self):
    self.ignored_atoms = {}
    self.active_atoms = []
    self.atom_names_given = []
    atom_id_str_pdbres_set = set()
    for atom in self.pdb_residue.atoms():
      atom_id_str_pdbres_set.add(atom.id_str(pdbres=True))
      if (atom.element.strip() == "Q"):
        self.ignored_atoms.setdefault(atom.name, []).append(atom)
      else:
        self.active_atoms.append(atom)
        self.atom_names_given.append(atom.name.replace(" ",""))
    return sorted(atom_id_str_pdbres_set)

  def _rna_sugar_pucker_analysis(self, params, next_pdb_residue):
    self.is_rna_dna = False
    self.is_rna2p = None
    if (self.monomer.is_peptide()): return
    from iotbx.pdb.rna_dna_detection import residue_analysis
    resname = self.pdb_residue.resname
    ra1 = residue_analysis(
      residue_atoms=self.pdb_residue.atoms(),
      distance_tolerance=params.bond_detection_distance_tolerance)
    if (ra1.problems is not None): return
    self.is_rna_dna = True
    if (not ra1.is_rna): return
    residue_2_p_atom = None
    if (next_pdb_residue is not None):
      residue_2_p_atom = next_pdb_residue.find_atom_by(name=" P  ")
    ana = rna_sugar_pucker_analysis.evaluate(
      params=params,
      residue_1_deoxy_ribo_atom_dict=ra1.deoxy_ribo_atom_dict,
      residue_1_c1p_outbound_atom=ra1.c1p_outbound_atom,
      residue_2_p_atom=residue_2_p_atom)
    self.is_rna2p = ana.is_2p
    mod_resname = mmtbx.geometry_restraints.torsion_restraints.utils.\
        modernize_rna_resname(resname=resname)
    if (self.is_rna2p):
      if mod_resname in ['  A', '  G']:
        primary_mod_id = "rna2p_pur"
      elif mod_resname in ['  C', '  U']:
        primary_mod_id = "rna2p_pyr"
      else:
        primary_mod_id = "rna2p"
    else:
      if mod_resname in ['  A', '  G']:
        primary_mod_id = "rna3p_pur"
      elif mod_resname in ['  C', '  U']:
        primary_mod_id = "rna3p_pyr"
      else:
        primary_mod_id = "rna3p"
    self.monomer, chem_mod_ids = self.mon_lib_srv.get_comp_comp_id_mod(
      comp_comp_id=self.monomer,
      mod_ids=(primary_mod_id,))
    self._track_mods(chem_mod_ids=chem_mod_ids)

  def _get_mappings(self):
    self.monomer_atom_dict = atom_dict = self.monomer.atom_dict()
    deuterium_aliases = None
    processed_atom_names = {}
    self.expected_atoms = {}
    self.unexpected_atoms = {}
    self.duplicate_atoms = {}
    if (self.atom_name_interpretation is not None):
      replace_primes = False
    elif (self.is_rna_dna or self.monomer.is_rna_dna()):
      replace_primes = True
    else:
      n_primes = 0
      n_stars = 0
      for atom_name in self.atom_names_given:
        if (atom_name.find("'") >= 0): n_primes += 1
        if (atom_name.find("*") >= 0): n_stars += 1
      replace_primes = (n_primes != 0 and n_stars == 0)
    _ = "".join([_ for _ in atom_dict.keys()])
    handle_case_insensitive = (_.upper() == _ or _.lower() == _)
    rna_dna_bb_cif_by_ref = None
    for i_atom,atom in enumerate(self.active_atoms):
      atom_name_given = self.atom_names_given[i_atom]
      if (handle_case_insensitive):
        atom_name_given = atom_name_given.upper()
      if (self.mon_lib_names is None):
        atom_name = atom_name_given
      else:
        atom_name = self.mon_lib_names[i_atom]
        if (atom_name is None):
          atom_name = atom_name_given
      if (len(atom_name) != 0 and atom_name not in atom_dict):
        auto_synomyms = []
        if (atom_name[0] in string.digits):
          auto_synomyms.append(atom_name[1:] + atom_name[0])
        elif (atom_name[-1] in string.digits):
          auto_synomyms.append(atom_name[-1] + atom_name[0:-1])
        if (replace_primes):
          atom_name = atom_name.replace("'", "*")
          if (atom_name != atom_name_given):
            auto_synomyms.append(atom_name)
            if (atom_name[0] in string.digits):
              auto_synomyms.append(atom_name[1:] + atom_name[0])
            elif (atom_name[-1] in string.digits):
              auto_synomyms.append(atom_name[-1] + atom_name[0:-1])
        for atom_name in auto_synomyms:
          if (atom_name in atom_dict): break
        else:
          auto_synomyms.insert(0, atom_name_given)
          if (deuterium_aliases is None):
            deuterium_aliases = self.monomer.hydrogen_deuterium_aliases()
          for atom_name in auto_synomyms:
            atom_name = deuterium_aliases.get(atom_name)
            if (atom_name is not None): break
          else:
            for atom_name in auto_synomyms:
              atom_name = self.mon_lib_srv.comp_synonym_atom_list_dict.get(
                self.monomer.chem_comp.id, {}).get(atom_name, None)
              if (atom_name is not None): break
            else:
              atom_name = atom_name_given
      if (    len(atom_name) != 0
          and atom_name not in atom_dict
          and ((self.is_rna_dna) or (self.monomer.normalized_rna_dna))):
        aliases = pdb.rna_dna_atom_names_backbone_aliases
        if (rna_dna_bb_cif_by_ref is None):
          rna_dna_bb_cif_by_ref = {}
          for cif_name in atom_dict:
            ref_name = aliases.get(cif_name)
            if (ref_name is not None):
              rna_dna_bb_cif_by_ref[ref_name] = cif_name
        ref_name = aliases.get(atom_name)
        cif_name = rna_dna_bb_cif_by_ref.get(ref_name)
        if (cif_name is not None):
          atom_name = cif_name
      prev_atom = processed_atom_names.get(atom_name)
      if (prev_atom is None):
        processed_atom_names[atom_name] = atom
        if (atom_name in atom_dict):
          self.expected_atoms[atom_name] = atom
        else:
          self.unexpected_atoms[atom_name] = atom
      else:
        self.duplicate_atoms.setdefault(atom_name, []).append(atom)
    if (    self.monomer.is_peptide()
        and self.atom_name_interpretation is None):
      self._rename_ot1_ot2("OXT" in atom_dict)
      self._auto_alias_h_h1()
    self._set_missing_atoms()

  def _rename_ot1_ot2(self, oxt_in_atom_dict):
    if ("O" not in self.expected_atoms):
      i_seq = self.unexpected_atoms.get("OT1", None)
      if (i_seq is not None):
        self.expected_atoms["O"] = i_seq
        del self.unexpected_atoms["OT1"]
    if (oxt_in_atom_dict):
      oxt_dict = self.expected_atoms
    else:
      oxt_dict = self.unexpected_atoms
    if ("OXT" not in oxt_dict):
      i_seq = self.unexpected_atoms.get("OT2", None)
      if (i_seq is not None):
        oxt_dict["OXT"] = i_seq
        del self.unexpected_atoms["OT2"]

  def _auto_alias_h_h1(self):
    if (self.monomer_atom_dict.get("H1") is None):
      return
    e = self.expected_atoms
    if ("H1" in e or "D1" in e):
      return
    u = self.unexpected_atoms
    h = u.get("H")
    d = u.get("D")
    key = "H"
    if (h is None):
      key = "D"
      h = d
    elif (d is not None):
      h = None
    if (h is not None):
      e["H1"] = h
      del u[key]

  def _set_missing_atoms(self):
    self.missing_non_hydrogen_atoms = {}
    self.missing_hydrogen_atoms = {}
    for atom in self.monomer.atom_list:
      if (atom.atom_id not in self.expected_atoms):
        if (atom.type_symbol != "H"):
          self.missing_non_hydrogen_atoms[atom.atom_id] = atom
        else:
          self.missing_hydrogen_atoms[atom.atom_id] = atom

  def _get_incomplete_info(self):
    if (    len(self.unexpected_atoms) == 0
        and len(self.missing_non_hydrogen_atoms) > 0):
      if (self.monomer.is_peptide()):
        atom_ids = list(self.expected_atoms.keys())
        atom_ids.sort()
        atom_ids = " ".join(atom_ids)
        if (atom_ids == "CA"): return "c_alpha_only"
        if (atom_ids == "C CA N"): return "n_c_alpha_c_only"
        if (atom_ids == "C CA N O"): return "backbone_only"
        if (atom_ids == "C CA CB N O"): return "truncation_to_alanine"
      elif (self.is_rna_dna or self.monomer.is_rna_dna()):
        atom_ids = " ".join(self.expected_atoms.keys())
        if (atom_ids == "P"): return "p_only"
    return None

  def _set_incomplete_info(self):
    self.incomplete_info = self._get_incomplete_info()

  def resolve_unexpected(self):
    mod_dict = self.mon_lib_srv.mod_mod_id_dict
    mod_mod_ids = []
    ani = self.atom_name_interpretation
    u = self.unexpected_atoms
    if (self.monomer.classification == "peptide"):
      if (ani is not None):
        u_mon_lib = {}
        for given_name,mon_lib_name in zip(ani.atom_names,
                                           self.mon_lib_names):
          i_seq = u.get(given_name)
          if (i_seq is None): continue
          # special case for terminating breaks with HC hydrogen
          if given_name in ["HC"] and "OC" not in ani.atom_names:
            u_mon_lib[given_name]=i_seq
          elif (mon_lib_name is None):
            u_mon_lib[given_name] = i_seq
          else:
            u_mon_lib[mon_lib_name] = i_seq
        u = u_mon_lib
      if ("HXT" in u):
        mod_mod_ids.append(mod_dict["COOH"])
      elif ("OXT" in u):
        mod_mod_ids.append(mod_dict["COO"])
      elif ("HC" in u):
        mod_mod_ids.append(mod_dict["CF-COH"])
      if (self.monomer.chem_comp.id == "GLU"):
        if ("HE2" in u):
          mod_mod_ids.append(mod_dict["ACID-GLU"])
      elif (self.monomer.chem_comp.id == "ASP"):
        if ("HD2" in u):
          mod_mod_ids.append(mod_dict["ACID-ASP"])
      def raise_missing_notpro(id, h):
        raise RuntimeError("""\
A modified version of the monomer library is required to correctly
handle N-terminal hydrogens. The mod_%sNOTPRO modification is missing.
This is a copy of mod_%s, but without the %s-N-CD angle.
Please contact cctbx@cci.lbl.gov for more information.""" % (id, id, h))
      if (ani is not None):
        nitrogen_hydrogens = []
        for name in u.keys():
          # name is a mon_lib_name
          if (name in ["H1", "H2", "H3"]):
            nitrogen_hydrogens.append(name)
        nitrogen_hydrogen_translation = None
        if (len(nitrogen_hydrogens) == 3):
          if (self.monomer_atom_dict.get("H") is not None):
            mod_mod_ids.append(mod_dict["NH3"])
        elif (len(nitrogen_hydrogens) == 2):
          if (self.monomer.chem_comp.id == "PRO"):
            mod_mod_id = mod_dict["NH2"]
          else:
            mod_mod_id = mod_dict.get("NH2NOTPRO")
            if (mod_mod_id is None):
              raise_missing_notpro("NH2", "HN2")
          mod_mod_ids.append(mod_mod_id)
          nitrogen_hydrogen_translation = ["HN1", "HN2"]
        elif (len(nitrogen_hydrogens) == 1):
          if (self.monomer.chem_comp.id == "PRO"):
            mod_mod_id = mod_dict["NH1"]
          else:
            mod_mod_id = mod_dict["NH1NOTPRO"]
            if (mod_mod_id is None):
              raise_missing_notpro("NH1", "HN")
          mod_mod_ids.append(mod_mod_id)
          nitrogen_hydrogen_translation = ["HN"]
        if (nitrogen_hydrogen_translation is not None):
          j = 0
          for i,mon_lib_name in enumerate(self.mon_lib_names):
            if (not mon_lib_name in u): continue
            if (mon_lib_name in ["H1", "H2", "H3"]):
              self.mon_lib_names[i] = nitrogen_hydrogen_translation[j]
              j += 1
          assert j == len(nitrogen_hydrogen_translation)
      else:
        if (      ("H1" in u or "1H" in u)
              and ("H2" in u or "2H" in u)
              and ("H3" in u or "3H" in u)):
          if (self.monomer_atom_dict.get("H") is not None):
            mod_mod_ids.append(mod_dict["NH3"])
        elif (    ("HN1" in u or "1HN" in u)
              and ("HN2" in u or "2HN" in u)):
          mod_mod_ids.append(mod_dict["NH2"])
        elif (    "HN" in u):
          mod_mod_ids.append(mod_dict["NH1"])
        elif (    ("H1" in u or "1H" in u)
              or  ("H2" in u or "2H" in u)
              or  ("H3" in u or "3H" in u)):
          if (self.monomer_atom_dict.get("H") is not None):
            mod_mod_ids.append(mod_dict["NH3"])
        elif (    ("HN1" in u or "1HN" in u)
              or  ("HN2" in u or "2HN" in u)):
          mod_mod_ids.append(mod_dict["NH2"])
    elif (self.monomer.classification in ["RNA", "DNA"]):
      if (ani is not None):
        if (ani.have_op3_or_hop3):
          mod_mod_ids.append(mod_dict["p5*END"])
        elif (not ani.have_phosphate):
          mod_mod_ids.append(mod_dict["5*END"])
        if (ani.have_ho3prime):
          mod_mod_ids.append(mod_dict["3*END"])
      else:
        if ("O3T" in u):
          mod_mod_ids.append(mod_dict["p5*END"])
        else:
          e = self.expected_atoms
          if (    not "P"   in e
              and not "OP1" in e
              and not "OP2" in e):
            mod_mod_ids.append(mod_dict["5*END"])
        if ("HO3*" in u):
          mod_mod_ids.append(mod_dict["3*END"])
    for mod_mod_id in mod_mod_ids:
      self.apply_mod(mod_mod_id=mod_mod_id)

  def _track_mods(self, chem_mod_ids):
    for chem_mod_id in chem_mod_ids:
      self.chem_mod_ids.add(chem_mod_id)
      self.residue_name += "%" + chem_mod_id

  def apply_mod(self, mod_mod_id):
    if (mod_mod_id.chem_mod.id in self.chem_mod_ids):
      return
        # mod previously applied already, e.g. two links to same carbohydrate
    try:
      mod_mon = self.monomer.apply_mod(mod_mod_id)
    except Exception as e:
      import traceback
      msg = traceback.format_exc().splitlines()
      msg.extend([
        "apply_mod failure:",
        "  %s" % residue_id_str(self.pdb_residue),
        "  comp id: %s" % self.monomer.chem_comp.id,
        "  mod id: %s" % mod_mod_id.chem_mod.id])
      raise Sorry("\n".join(msg))
    self._track_mods(chem_mod_ids=[mod_mod_id.chem_mod.id])
    mod_mon.classification = self.monomer.classification
    self.monomer = mod_mon
    if (    mod_mod_id.chem_mod.name is not None
        and mod_mod_id.chem_mod.name.lower().find("terminus") >= 0):
      self.is_terminus = True # AD HOC manipulation
    self._get_mappings()

  def residue_altloc(self):
    result = self.residue_name
    if (self.conf_altloc != ""):
      result += ', conformer "%s"' % self.conf_altloc
    return result

  def _is_unusual(self):
    m = self.monomer
    if (m is None): return True
    if (m.is_peptide()): return False
    if (m.is_rna_dna()): return False
    if (m.is_water()): return False
    return True

  def summary(self):
    if (self.monomer is None):
      classification = None
    else:
      classification = self.monomer.classification
    # TODO: verify all these .values calls are dict method otherwise remove the list()
    return monomer_mapping_summary(
      conf_altloc=self.conf_altloc,
      residue_name=self.residue_name,
      expected_atoms=list(self.expected_atoms.values()),
      unexpected_atoms=list(self.unexpected_atoms.values()),
      duplicate_atoms=flat_list(list(self.duplicate_atoms.values())),
      ignored_atoms=list(self.ignored_atoms.values()),
      classification=classification,
      incomplete_info=self.incomplete_info,
      is_terminus=self.is_terminus,
      is_unusual=self._is_unusual())

  def add_bond_proxies(self,
                       bond_simple_proxy_registry,
                       use_neutron_distances=False,
                       ):
    self.bond_counters = add_bond_proxies(
      counters=counters(label="bond"),
      m_i=self,
      m_j=self,
      bond_list=self.monomer.bond_list,
      bond_simple_proxy_registry=bond_simple_proxy_registry,
      use_neutron_distances=use_neutron_distances).counters

  def add_angle_proxies(self, special_position_dict, angle_proxy_registry):
    self.angle_counters = add_angle_proxies(
      counters=counters(label="angle"),
      m_i=self,
      m_j=None,
      angle_list=self.monomer.angle_list,
      angle_proxy_registry=angle_proxy_registry,
      special_position_dict=special_position_dict).counters

  def add_dihedral_proxies(self,
        dihedral_function_type,
        special_position_dict,
        dihedral_proxy_registry):
    self.dihedral_counters = add_dihedral_proxies(
      counters=counters(label="dihedral"),
      m_i=self,
      m_j=None,
      tor_list=self.monomer.tor_list,
      dihedral_function_type=dihedral_function_type,
      peptide_link_params=None,
      dihedral_proxy_registry=dihedral_proxy_registry,
      special_position_dict=special_position_dict).counters

  def add_chirality_proxies(self, special_position_dict,
                                  chirality_proxy_registry,
                                  chir_volume_esd):
    self.chirality_counters = add_chirality_proxies(
      counters=counters(label="chirality"),
      m_i=self,
      m_j=None,
      chir_list=self.monomer.chir_list,
      chirality_proxy_registry=chirality_proxy_registry,
      special_position_dict=special_position_dict,
      chir_volume_esd=chir_volume_esd).counters

  def add_planarity_proxies(self, special_position_dict,
                                  planarity_proxy_registry):
    self.planarity_counters = add_planarity_proxies(
      counters=counters(label="planarity"),
      m_i=self,
      m_j=None,
      plane_list=self.monomer.get_planes(),
      planarity_proxy_registry=planarity_proxy_registry,
      special_position_dict=special_position_dict).counters

  def add_parallelity_proxies(self, special_position_dict,
                                    parallelity_proxy_registry):
    self.parallelity_counters = add_parallelity_proxies(
      counters=counters(label="parallelity"),
      m_i=self,
      m_j=self,
      plane_list=self.monomer.get_planes(),
      planarity_proxy_registry=planarity_proxy_registry,
      special_position_dict=special_position_dict).counters

class link_match_one(object):

  def __init__(self, chem_link_comp_id, chem_link_group_comp,
                     comp_id, comp_group):
    if (comp_group in cif_types.peptide_comp_groups):
      comp_group = "peptide"
    elif (comp_group in cif_types.dna_rna_comp_groups):
      comp_group = "DNA/RNA"
    if (   chem_link_comp_id in [None, ""]
        or (chem_link_comp_id is not None
            and comp_id.lower() == chem_link_comp_id.lower())):
      self.is_comp_id_match = True
      if (chem_link_comp_id is None):
        self.len_comp_id_match = 0
      else:
        self.len_comp_id_match = len(chem_link_comp_id)
    else:
      self.is_comp_id_match = False
      self.len_comp_id_match = -1
    if (   chem_link_group_comp in [None, ""]
        or (comp_group is not None
            and comp_group.lower() == chem_link_group_comp.lower())):
      self.is_group_match = True
      if (chem_link_group_comp is None):
        self.len_group_match = 0
      else:
        self.len_group_match = len(chem_link_group_comp)
    elif ( comp_group in [None, ""] ):
      # linking non-standard amino acids
      if self.is_comp_id_match and self.len_comp_id_match>0:
        self.is_group_match = True
        self.len_group_match = len(chem_link_group_comp)
      else:
        self.is_group_match = False
        self.len_group_match = -1
    else:
      self.is_group_match = False
      self.len_group_match = 0

  def is_match(self):
    return self.is_comp_id_match and self.is_group_match

class link_match(object):

  def __init__(self, link_link_id, comp_id_1, comp_group_1,
                                   comp_id_2, comp_group_2):
    self.link_link_id = None
    chem_link = link_link_id.chem_link
    match_1 = link_match_one(
      chem_link.comp_id_1, chem_link.group_comp_1,
      comp_id_1, comp_group_1)
    if (not match_1.is_match()): return
    match_2 = link_match_one(
      chem_link.comp_id_2, chem_link.group_comp_2,
      comp_id_2, comp_group_2)
    if (not match_2.is_match()): return
    self.link_link_id = link_link_id
    self.len_comp_id_match_1 = match_1.len_comp_id_match
    self.len_group_match_1 = match_1.len_group_match
    self.len_comp_id_match_2 = match_2.len_comp_id_match
    self.len_group_match_2 = match_2.len_group_match

  def is_proper_match(self):
    return (
         self.len_comp_id_match_1 > 0
      or self.len_group_match_1 > 0
      or self.len_comp_id_match_2 > 0
      or self.len_group_match_2 > 0)

  def __lt__(self, other):
    if (self.n_unresolved_bonds < other.n_unresolved_bonds): return 1
    if (self.n_unresolved_angles < other.n_unresolved_angles): return 1
    if (self.len_comp_id_match_1 > other.len_comp_id_match_1): return 1
    if (self.len_comp_id_match_2 > other.len_comp_id_match_2): return 1
    if (self.len_group_match_1 > other.len_group_match_1): return 1
    if (self.len_group_match_2 > other.len_group_match_2): return 1
    return 0

  def __gt__(self, other):
    if (self.n_unresolved_bonds > other.n_unresolved_bonds): return  1
    if (self.n_unresolved_angles > other.n_unresolved_angles): return  1
    if (self.len_comp_id_match_1 < other.len_comp_id_match_1): return  1
    if (self.len_comp_id_match_2 < other.len_comp_id_match_2): return  1
    if (self.len_group_match_1 < other.len_group_match_1): return  1
    if (self.len_group_match_2 < other.len_group_match_2): return  1
    return 0

  def __eq__(self, other):
    if __lt__(other) == 0 and __gt__(other) == 0:
      return 1
    return 0

  def __cmp__(self, other):
    return self.__lt__(other) - self.__gt__(other)

def get_lib_link_peptide(mon_lib_srv, m_i, m_j, include_peptide_plane=False):
  link_id = "TRANS"
  if include_peptide_plane:
    link_id = "PEPTIDE-PLANE"
  elif (m_j.expected_atoms.get("CN", None) is not None):
    link_id = "NM" + link_id
  elif (m_j.monomer.chem_comp.id == "PRO"):
    link_id = "P" + link_id
  return mon_lib_srv.link_link_id_dict[link_id]

def get_lib_link(mon_lib_srv,
                 m_i,
                 m_j,
                 include_peptide_plane=False,
                 verbose=False):
  if (m_i.monomer.is_water() or m_j.monomer.is_water()): return None
  if (m_i.monomer.is_peptide() and m_j.monomer.is_peptide()):
    if verbose: print('peptide-peptide')
    return get_lib_link_peptide(mon_lib_srv,
                                m_i,
                                m_j,
                                include_peptide_plane=include_peptide_plane)
  elif (    (m_i.is_rna_dna or m_i.monomer.is_rna_dna())
        and (m_j.is_rna_dna or m_j.monomer.is_rna_dna())):
    if (m_i.is_rna2p):
      if verbose: print('rna2p')
      return mon_lib_srv.link_link_id_dict["rna2p"]
    if verbose: print('rna3p')
    return mon_lib_srv.link_link_id_dict["rna3p"]
  comp_id_1 = m_i.monomer.chem_comp.id
  comp_id_2 = m_j.monomer.chem_comp.id
  comp_1 = mon_lib_srv.get_comp_comp_id_direct(comp_id_1)
  comp_2 = mon_lib_srv.get_comp_comp_id_direct(comp_id_2)
  group_1 = comp_1.chem_comp.group
  group_2 = comp_2.chem_comp.group
  matches = []
  for link_link_id in mon_lib_srv.link_link_id_list:
    chem_link = link_link_id.chem_link
    if (chem_link.name in cif_types.non_chain_links): continue
    if (    chem_link.comp_id_1 == ""
        and chem_link.mod_id_1 == ""
        and chem_link.group_comp_1 == ""
        and chem_link.comp_id_2 == ""
        and chem_link.mod_id_2 == ""
        and chem_link.group_comp_2 == ""): continue
    match = link_match(link_link_id, comp_id_1, group_1, comp_id_2, group_2)
    if (match.link_link_id is not None):
      def get_atom(restr, comp_id, atom_id):
        ad = None
        if (comp_id == 1): ad = m_i.monomer_atom_dict
        if (comp_id == 2): ad = m_j.monomer_atom_dict
        if (ad is None):
          raise Sorry("""\
Corrupt CIF link definition:
  source info: %s
  link id: %s
  link name: %s
  %s atom atom_id: %s
  %s atom comp_id: %d (must be 1 or 2)""" % (
            str(match.link_link_id.source_info),
            chem_link.id,
            chem_link.name,
            restr, atom_id,
            restr, comp_id))
        atom = ad.get(atom_id, None)
        return atom
      match.n_unresolved_bonds = 0
      for bond in match.link_link_id.bond_list:
        atoms = [
          get_atom("bond", bond.atom_1_comp_id, bond.atom_id_1),
          get_atom("bond", bond.atom_2_comp_id, bond.atom_id_2)]
        if (None in atoms):
          match.n_unresolved_bonds += 1
      match.n_unresolved_angles = 0
      for angle in match.link_link_id.angle_list:
        atoms = [
          get_atom("angle", angle.atom_1_comp_id, angle.atom_id_1),
          get_atom("angle", angle.atom_2_comp_id, angle.atom_id_2),
          get_atom("angle", angle.atom_3_comp_id, angle.atom_id_3)]
        if (None in atoms):
          match.n_unresolved_angles += 1
      matches.append(match)
  if (len(matches) == 0): return None
  matches.sort()
  best_matches = []
  def _show_match(match):
    print('match '*10)
    print(match.link_link_id.source_info)
    match.link_link_id.chem_link.show()
    for attr in [
      "n_unresolved_bonds",
      "n_unresolved_angles",
      "len_comp_id_match_1",
      "len_comp_id_match_2",
      "len_group_match_1",
      "len_group_match_2",
      ]:
      print(attr, getattr(match, attr, None))
    print('_'*80)
  for m in matches:
    if verbose: _show_match(m)
    if (cmp(m, matches[0]) != 0): break
    best_matches.append(m)
  match = best_matches[0]
  if (not match.is_proper_match()):
    return None
  return match.link_link_id

def get_restraints_loading_flags(params):
  rc = {}
  if params:
    rc["use_neutron_distances"] = params.use_neutron_distances
  return rc

def evaluate_registry_process_result(
      proxy_label,
      m_i, m_j, i_seqs,
      registry_process_result,
      lines=[]):
  if (registry_process_result.is_conflicting):
    raise Sorry(format_exception_message(
      m_i=m_i,
      m_j=m_j,
      i_seqs=i_seqs,
      base_message="Conflicting %s restraints:" % proxy_label,
      source_labels=registry_process_result.conflict_source_labels,
      show_residue_names=False,
      lines=lines))
  pdb_atoms = m_i.pdb_atoms
  atoms = [pdb_atoms[i_seq] for i_seq in i_seqs]
  if (not registry_process_result.is_new
      and not all_atoms_are_in_main_conf(atoms=atoms)):
    raise Sorry(format_exception_message(
      m_i=m_i,
      m_j=m_j,
      i_seqs=i_seqs,
      base_message="Duplicate %s restraints:" % proxy_label,
      lines=lines))

class add_bond_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        bond_list,
        bond_simple_proxy_registry,
        sites_cart=None,
        distance_cutoff=None,
        use_neutron_distances=False,
        ):
    if (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    self.counters = counters
    self.broken_bond_i_seq_pairs = set()
    value = "value_dist"
    for bond in bond_list:
      if use_neutron_distances:
        value = "value_dist"
        if getattr(bond, "value_dist_neutron", None):
          value = "value_dist_neutron"
      if (   bond.atom_id_1 not in m_i.monomer_atom_dict
          or bond.atom_id_2 not in m_j.monomer_atom_dict):
        #
        # replace primes with stars to see if that will work!!!
        #
        if (    bond.atom_id_1.replace("'", "*") in m_i.monomer_atom_dict
            and bond.atom_id_2.replace("'", "*") in m_j.monomer_atom_dict):
          bond.atom_id_1 = bond.atom_id_1.replace("'", "*")
          bond.atom_id_2 = bond.atom_id_2.replace("'", "*")
        else:
          counters.corrupt_monomer_library_definitions += 1
          continue
      atoms = (m_i.expected_atoms.get(bond.atom_id_1, None),
               m_j.expected_atoms.get(bond.atom_id_2, None))
      if (None in atoms):
        if (   m_i.monomer_atom_dict[bond.atom_id_1].type_symbol == "H"
            or m_j.monomer_atom_dict[bond.atom_id_2].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif ( getattr(bond, value) is None
                  or bond.value_dist_esd in [None, 0]):
        counters.undefined += 1
      else:
        counters.resolved += 1
        i_seqs = [atom.i_seq for atom in atoms]
        proxy = geometry_restraints.bond_simple_proxy(
          i_seqs=i_seqs,
          distance_ideal=getattr(bond, value),
          weight=1/bond.value_dist_esd**2)
        is_large_distance = False
        if (sites_cart is not None):
          r = geometry_restraints.bond(sites_cart=sites_cart, proxy=proxy)
          if (r.distance_model > distance_cutoff):
            is_large_distance = True
            self.broken_bond_i_seq_pairs.add(tuple(sorted(i_seqs)))
        if (not is_large_distance):
          if (    not m_i.is_first_conformer_in_chain
              and all_atoms_are_in_main_conf(atoms=atoms)):
            counters.already_assigned_to_first_conformer += 1
          else:
            registry_process_result = bond_simple_proxy_registry.process(
              source_info=source_info_server(m_i=m_i, m_j=m_j),
              proxy=proxy)
            evaluate_registry_process_result(
              proxy_label="bond_simple", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
              registry_process_result=registry_process_result)

class add_angle_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        angle_list,
        angle_proxy_registry,
        special_position_dict,
        broken_bond_i_seq_pairs=None):
    self.counters = counters
    if (m_j is None):
      m_1,m_2,m_3 = m_i,m_i,m_i
    elif (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    for angle in angle_list:
      if (m_j is not None):
        m_1,m_2,m_3 = [(m_i, m_j)[comp_id-1] for comp_id in (
          angle.atom_1_comp_id, angle.atom_2_comp_id, angle.atom_3_comp_id)]
      if (   angle.atom_id_1 not in m_1.monomer_atom_dict
          or angle.atom_id_2 not in m_2.monomer_atom_dict
          or angle.atom_id_3 not in m_3.monomer_atom_dict):
        if (    angle.atom_id_1.replace("'", "*") in m_1.monomer_atom_dict
            and angle.atom_id_2.replace("'", "*") in m_2.monomer_atom_dict
            and angle.atom_id_3.replace("'", "*") in m_3.monomer_atom_dict):
          angle.atom_id_1 = angle.atom_id_1.replace("'", "*")
          angle.atom_id_2 = angle.atom_id_2.replace("'", "*")
          angle.atom_id_3 = angle.atom_id_3.replace("'", "*")
        else:
          counters.corrupt_monomer_library_definitions += 1
          continue
      atoms = (m_1.expected_atoms.get(angle.atom_id_1, None),
               m_2.expected_atoms.get(angle.atom_id_2, None),
               m_3.expected_atoms.get(angle.atom_id_3, None))
      if (None in atoms):
        if (   m_1.monomer_atom_dict[angle.atom_id_1].type_symbol == "H"
            or m_2.monomer_atom_dict[angle.atom_id_2].type_symbol == "H"
            or m_3.monomer_atom_dict[angle.atom_id_3].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif (   angle.value_angle is None
            or angle.value_angle_esd in [None, 0]):
        counters.undefined += 1
      else:
        counters.resolved += 1
        i_seqs = [atom.i_seq for atom in atoms]
        if (special_position_dict.involves_special_positions(i_seqs)):
          counters.discarded_because_of_special_positions += 1
        elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
          pass
        else:
          registry_process_result = angle_proxy_registry.process(
            source_info=source_info_server(m_i=m_i, m_j=m_j),
            proxy=geometry_restraints.angle_proxy(
              i_seqs=i_seqs,
              angle_ideal=angle.value_angle,
              weight=1/angle.value_angle_esd**2))
          evaluate_registry_process_result(
            proxy_label="angle", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
            registry_process_result=registry_process_result)

class add_dihedral_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        tor_list,
        dihedral_function_type,
        peptide_link_params,
        dihedral_proxy_registry,
        special_position_dict,
        sites_cart=None,
        chem_link_id=None,
        broken_bond_i_seq_pairs=None,
        cis_trans_specifications=None,
        ):
    self.counters = counters
    self.chem_link_id = chem_link_id
    if (chem_link_id not in ["TRANS", "PTRANS", "NMTRANS"]):
      sites_cart = None
    if (m_j is None):
      m_1,m_2,m_3,m_4 = m_i,m_i,m_i,m_i
    elif (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    for tor in tor_list:
      if (m_j is not None):
        m_1,m_2,m_3,m_4 = [(m_i, m_j)[comp_id-1] for comp_id in (
          tor.atom_1_comp_id,
          tor.atom_2_comp_id,
          tor.atom_3_comp_id,
          tor.atom_4_comp_id)]
      if (   tor.atom_id_1 not in m_1.monomer_atom_dict
          or tor.atom_id_2 not in m_2.monomer_atom_dict
          or tor.atom_id_3 not in m_3.monomer_atom_dict
          or tor.atom_id_4 not in m_4.monomer_atom_dict):
        if (    tor.atom_id_1.replace("'", "*") in m_1.monomer_atom_dict
            and tor.atom_id_2.replace("'", "*") in m_2.monomer_atom_dict
            and tor.atom_id_3.replace("'", "*") in m_3.monomer_atom_dict
            and tor.atom_id_4.replace("'", "*") in m_4.monomer_atom_dict):
          tor.atom_id_1 = tor.atom_id_1.replace("'", "*")
          tor.atom_id_2 = tor.atom_id_2.replace("'", "*")
          tor.atom_id_3 = tor.atom_id_3.replace("'", "*")
          tor.atom_id_4 = tor.atom_id_4.replace("'", "*")
        else:
          counters.corrupt_monomer_library_definitions += 1
          continue
      atoms = (m_1.expected_atoms.get(tor.atom_id_1, None),
               m_2.expected_atoms.get(tor.atom_id_2, None),
               m_3.expected_atoms.get(tor.atom_id_3, None),
               m_4.expected_atoms.get(tor.atom_id_4, None))
      if (None in atoms):
        if (   m_1.monomer_atom_dict[tor.atom_id_1].type_symbol == "H"
            or m_2.monomer_atom_dict[tor.atom_id_2].type_symbol == "H"
            or m_3.monomer_atom_dict[tor.atom_id_3].type_symbol == "H"
            or m_4.monomer_atom_dict[tor.atom_id_4].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif (   tor.value_angle is None
            or tor.value_angle_esd in [None, 0]):
        counters.undefined += 1
      else:
        counters.resolved += 1
        i_seqs = [atom.i_seq for atom in atoms]
        trans_cis_ids = [
          "TRANS", "PTRANS", "NMTRANS",
          "CIS",   "PCIS",   "NMCIS"]
        if (special_position_dict.involves_special_positions(i_seqs)):
          counters.discarded_because_of_special_positions += 1
        elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
          pass
        elif (    tor.id in ["psi", "phi"]
              and self.chem_link_id in trans_cis_ids
              and peptide_link_params.discard_psi_phi):
          pass
        elif (    tor.id == "omega"
              and self.chem_link_id in trans_cis_ids
              and peptide_link_params.discard_omega):
          pass
        else:
          if (dihedral_function_type == "determined_by_sign_of_periodicity"):
            periodicity = tor.period
          elif (dihedral_function_type == "all_sinusoidal"):
            periodicity = max(1, tor.period)
          elif (dihedral_function_type == "all_harmonic"):
            periodicity = -abs(tor.period)
          else:
            raise RuntimeError(
              "Unknown dihedral_function_type: %s"
                % str(dihedral_function_type))
          try:
            if len(tor.alt_value_angle) == 0:
              alt_value_angle = None
            else:
              alt_value_angle = [float(t) for t in tor.alt_value_angle.split(",")]
          except Exception:
            alt_value_angle = None
          proxy = geometry_restraints.dihedral_proxy(
            i_seqs=i_seqs,
            angle_ideal=tor.value_angle,
            weight=1/tor.value_angle_esd**2,
            periodicity=periodicity, alt_angle_ideals=alt_value_angle)
          if (sites_cart is not None and tor.id == "omega"):
            assert abs(tor.value_angle - 180) < 1.e-6
            if (peptide_link_params.omega_esd_override_value is not None):
              assert peptide_link_params.omega_esd_override_value > 0
              proxy.weight = 1/peptide_link_params.omega_esd_override_value**2
            r = geometry_restraints.dihedral(
              sites_cart=sites_cart,
              proxy=proxy)
            if ( not peptide_link_params.apply_all_trans and
                 abs(r.delta) > 180-peptide_link_params.cis_threshold):
              self.chem_link_id = self.chem_link_id.replace("TRANS", "CIS")
              proxy.angle_ideal = 0
            if cis_trans_specifications:
              for ca_i_seq in cis_trans_specifications:
                if ca_i_seq[0] == i_seqs[3]: # specify the trailing CA
                  if self.chem_link_id in ["PTRANS", "PCIS"]:
                    cis_trans_specifications[ca_i_seq] = "P%s" % (
                      cis_trans_specifications[ca_i_seq],
                      )
                  elif self.chem_link_id in ["NMTRANS", "NMCIS"]:
                    cis_trans_specifications[ca_i_seq] = "NM%s" % (
                      cis_trans_specifications[ca_i_seq],
                      )
                  if self.chem_link_id!=cis_trans_specifications[ca_i_seq].upper():
                    if self.chem_link_id.find("TRANS")>-1:
                      self.chem_link_id=self.chem_link_id.replace('TRANS', 'CIS')
                      proxy.angle_ideal=0
                    else:
                      self.chem_link_id=self.chem_link_id.replace('CIS','TRANS')
                      proxy.angle_ideal=180

          registry_process_result = dihedral_proxy_registry.process(
            source_info=source_info_server(m_i=m_i, m_j=m_j),
            proxy=proxy)
          evaluate_registry_process_result(
            proxy_label="dihedral", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
            registry_process_result=registry_process_result,
            lines=["tor id: " + str(tor.id)])

class add_chirality_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        chir_list,
        chirality_proxy_registry,
        special_position_dict,
        chir_volume_esd,
        lib_link=None,
        broken_bond_i_seq_pairs=None):
    self.counters = counters
    self.counters.unsupported_volume_sign = dicts.with_default_value(0)
    if (m_j is None):
      m_c,m_1,m_2,m_3 = m_i,m_i,m_i,m_i
    elif (m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    for chir in chir_list:
      if (m_j is not None):
        m_c,m_1,m_2,m_3 = [(m_i, m_j)[comp_id-1] for comp_id in (
          chir.atom_centre_comp_id,
          chir.atom_1_comp_id,
          chir.atom_2_comp_id,
          chir.atom_3_comp_id)]
      volume_sign = chir.volume_sign
      if (volume_sign is not None):
        volume_sign = volume_sign[:4].lower()
      if (volume_sign not in ["posi", "nega", "both"]):
        counters.unsupported_volume_sign[volume_sign] += 1
        continue
      if (   chir.atom_id_centre not in m_c.monomer_atom_dict
          or chir.atom_id_1 not in m_1.monomer_atom_dict
          or chir.atom_id_2 not in m_2.monomer_atom_dict
          or chir.atom_id_3 not in m_3.monomer_atom_dict):
        if (    chir.atom_id_1.replace("'", "*") in m_1.monomer_atom_dict
            and chir.atom_id_2.replace("'", "*") in m_2.monomer_atom_dict
            and chir.atom_id_3.replace("'", "*") in m_3.monomer_atom_dict
            and chir.atom_id_centre.replace("'", "*") in m_c.monomer_atom_dict):
          chir.atom_id_1 = chir.atom_id_1.replace("'", "*")
          chir.atom_id_2 = chir.atom_id_2.replace("'", "*")
          chir.atom_id_3 = chir.atom_id_3.replace("'", "*")
          chir.atom_id_centre = chir.atom_id_centre.replace("'", "*")
        else:
          counters.corrupt_monomer_library_definitions += 1
          continue
      atoms = (m_c.expected_atoms.get(chir.atom_id_centre, None),
               m_1.expected_atoms.get(chir.atom_id_1, None),
               m_2.expected_atoms.get(chir.atom_id_2, None),
               m_3.expected_atoms.get(chir.atom_id_3, None))
      if (None in atoms):
        if (   m_c.monomer_atom_dict[chir.atom_id_centre].type_symbol == "H"
            or m_1.monomer_atom_dict[chir.atom_id_1].type_symbol == "H"
            or m_2.monomer_atom_dict[chir.atom_id_2].type_symbol == "H"
            or m_3.monomer_atom_dict[chir.atom_id_3].type_symbol == "H"):
          counters.unresolved_hydrogen += 1
        else:
          counters.unresolved_non_hydrogen += 1
      elif (   volume_sign is None
            or chir_volume_esd in [None, 0]):
        counters.undefined += 1
      else:
        if (m_j is None):
          volume_ideal = m_i.monomer.get_chir_volume_ideal(chir)
        else:
          volume_ideal = lib_link.get_chir_volume_ideal(
            m_i.monomer, m_j.monomer, chir)
        if (volume_ideal is None):
          counters.undefined += 1
        else:
          counters.resolved += 1
          i_seqs = [atom.i_seq for atom in atoms]
          if (special_position_dict.involves_special_positions(i_seqs)):
            counters.discarded_because_of_special_positions += 1
          elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
            pass
          else:
            registry_process_result = chirality_proxy_registry.process(
              source_info=source_info_server(m_i=m_i, m_j=m_j),
              proxy=geometry_restraints.chirality_proxy(
                i_seqs=i_seqs,
                volume_ideal=volume_ideal,
                both_signs=(volume_sign == "both"),
                weight=1/chir_volume_esd**2))
            evaluate_registry_process_result(
              proxy_label="chirality", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
              registry_process_result=registry_process_result)

class add_planarity_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        plane_list,
        planarity_proxy_registry,
        special_position_dict,
        peptide_link_params=None,
        broken_bond_i_seq_pairs=None):
    self.counters = counters
    self.counters.less_than_four_sites = dicts.with_default_value(0)
    if (    m_j is not None
        and m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    for plane in plane_list:
      this_plane_has_unresolved_non_hydrogen = False
      i_seqs = []
      weights = []
      for plane_atom in plane.plane_atoms:
        if (m_j is None):
          m_x = m_i
        else:
          assert plane_atom.atom_comp_id in (1,2)
          m_x = (m_i, m_j)[plane_atom.atom_comp_id-1]
        if (plane_atom.atom_id not in m_x.monomer_atom_dict):
          if (   plane_atom.atom_id.replace("'", "*") in m_x.monomer_atom_dict):
            plane_atom.atom_id = plane_atom.atom_id.replace("'", "*")
          else:
            counters.corrupt_monomer_library_definitions += 1
            continue
        atom = m_x.expected_atoms.get(plane_atom.atom_id, None)
        if (atom is None):
          if (m_x.monomer_atom_dict[plane_atom.atom_id].type_symbol == "H"):
            counters.unresolved_hydrogen += 1
          else:
            counters.unresolved_non_hydrogen += 1
            this_plane_has_unresolved_non_hydrogen = True
        elif (plane_atom.dist_esd in [None, 0]):
          counters.undefined += 1
        else:
          counters.resolved += 1
          i_seq = atom.i_seq
          if (special_position_dict.involves_special_positions([i_seq])):
            counters.discarded_because_of_special_positions += 1
          else:
            i_seqs.append(i_seq)
            weights.append(1/plane_atom.dist_esd**2)
      if (len(i_seqs) < 4):
        if (this_plane_has_unresolved_non_hydrogen):
          counters.less_than_four_sites[plane.plane_id] += 1
      elif (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs)):
        pass
      else:
        registry_process_result = planarity_proxy_registry.process(
          source_info=source_info_server(m_i=m_i, m_j=m_j),
          proxy=geometry_restraints.planarity_proxy(
            i_seqs=flex.size_t(i_seqs),
            weights=flex.double(weights)))
        evaluate_registry_process_result(
          proxy_label="planarity", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
          registry_process_result=registry_process_result,
          lines=["plane id: " + str(plane.plane_id)])

class add_parallelity_proxies(object):

  def __init__(self,
        counters,
        m_i,
        m_j,
        parallelity_proxy_registry,
        special_position_dict,
        broken_bond_i_seq_pairs=None,
        weight=0.05):
    # probably is not used anymore
    if weight <=0:
      raise Sorry("Weight for parallelity restraint should be > 0.")
    self.counters = counters
    if (    m_j is not None
        and m_i.i_conformer != 0 and m_j.i_conformer != 0):
      assert m_i.i_conformer == m_j.i_conformer
    counters.resolved += 1
    # making i_seqs, j_seqs,weight
    i_seqs = []
    j_seqs = []
    for seqs, m in [(i_seqs, m_i), (j_seqs, m_j)]:
      for p in m.monomer.get_planes():
        for plane_atom in p.plane_atoms:
          atom = m.expected_atoms.get(plane_atom.atom_id, None)
          if atom is not None and atom.i_seq not in seqs:
            seqs.append(atom.i_seq)
    if (involves_broken_bonds(broken_bond_i_seq_pairs, i_seqs+j_seqs)):
      pass
    elif len(i_seqs) < 3 or len(j_seqs) < 3:
      pass
    else:
      registry_process_result = parallelity_proxy_registry.process(
        source_info=source_info_server(m_i=m_i, m_j=m_j),
        proxy=geometry_restraints.parallelity_proxy(
          i_seqs=flex.size_t(i_seqs),
          j_seqs=flex.size_t(j_seqs),
          weight=weight))
      evaluate_registry_process_result(
        proxy_label="parallelity", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
        registry_process_result=registry_process_result,
          lines=["plane id: " + "???"])


# XXX TODO synonymes
def ener_lib_as_nonbonded_params(
      ener_lib,
      assume_hydrogens_all_missing,
      factor_1_4_interactions,
      default_distance,
      minimum_distance,
      const_shrink_donor_acceptor,
      use_lib_vdw=False):
  params = geometry_restraints.nonbonded_params(
    factor_1_4_interactions=factor_1_4_interactions,
    const_shrink_1_4_interactions=0,
    default_distance=default_distance,
    minimum_distance=minimum_distance,
    const_shrink_donor_acceptor=const_shrink_donor_acceptor)
  if (use_lib_vdw):
    tables = {"": [], "h": []}
    for vdw in ener_lib.lib_vdw:
      assert vdw.H_flag in ["", "h"]
      if (vdw.H_flag == ""):
        tables[""].append(vdw)
      else:
        tables["h"].append(vdw)
    if (assume_hydrogens_all_missing):
      reverse_prefs = ["", "h"]
    else:
      reverse_prefs = ["h", ""]
    for code in reverse_prefs:
      for vdw in tables[code]:
        atom_types = [vdw.atom_type_1, vdw.atom_type_2]
        atom_types.sort()
        params.distance_table.setdefault(
          atom_types[0])[atom_types[1]] = vdw.radius_min
  if (assume_hydrogens_all_missing):
    pref1, pref2 = ["vdwh_radius", "vdw_radius"]
  else:
    pref1, pref2 = ["vdw_radius", "vdwh_radius"]
  for atom_type,energy_lib_atom in ener_lib.lib_atom.items():
    if (len(atom_type) == 0): continue
    r = getattr(energy_lib_atom, pref1)
    if (r is None):
      r = getattr(energy_lib_atom, pref2)
    if (r is not None):
      params.radius_table[atom_type] = r
    r_ionic = getattr(energy_lib_atom, "ion_radius")
    if (r_ionic is not None):
      params.ionic_radius_table[atom_type] = r_ionic
    # N = 0, D = 1, A = 2, B = 3, H = 4
    if getattr(energy_lib_atom, "hb_type") == 'N':
      params.donor_acceptor_table[atom_type] = 0
    elif getattr(energy_lib_atom, "hb_type") == 'D':
      params.donor_acceptor_table[atom_type] = 1
    elif getattr(energy_lib_atom, "hb_type") == 'A':
      params.donor_acceptor_table[atom_type] = 2
    elif getattr(energy_lib_atom, "hb_type") == 'B':
      params.donor_acceptor_table[atom_type] = 3
    elif getattr(energy_lib_atom, "hb_type") == 'H':
      params.donor_acceptor_table[atom_type] = 4
  return params

def is_same_model_as_before(model_type_indices, i_model, models):
  m_i = models[i_model]
  for j_model in range(0, i_model):
    if (model_type_indices[j_model] != j_model): continue
    if (m_i.is_identical_hierarchy(other=models[j_model])):
      model_type_indices[i_model] = j_model
      return True
  model_type_indices[i_model] = i_model
  return False

class build_chain_proxies(object):

  def __init__(self,
        mon_lib_srv,
        ener_lib,
        translate_cns_dna_rna_residue_names,
        rna_sugar_pucker_analysis_params,
        apply_cif_modifications,
        apply_cif_links_mm_pdbres_dict,
        link_distance_cutoff,
        not_linked_show_max,
        dihedral_function_type,
        chir_volume_esd,
        peptide_link_params,
        pdb_hierarchy,
        pdb_atoms,
        sites_cart,
        special_position_dict,
        keep_monomer_mappings,
        all_monomer_mappings,
        scattering_type_registry,
        nonbonded_energy_type_registry,
        geometry_proxy_registries,
        cystein_sulphur_i_seqs,
        cystein_monomer_mappings,
        is_unique_model,
        i_model,
        i_conformer,
        is_first_conformer_in_chain,
        conformer,
        conformation_dependent_restraints_list,
        cis_trans_specifications,
        apply_restraints_specifications,
        log,
        restraints_loading_flags=None,
        fatal_problem_max_lines=10,
               ):
    if restraints_loading_flags is None: restraints_loading_flags={}
    self._cif = cif_output_holder()
    self.pdb_link_records = {}
    self.conformation_dependent_restraints_list = \
      conformation_dependent_restraints_list
    unknown_residues = dicts.with_default_value(0)
    ad_hoc_single_atom_residues = dicts.with_default_value(0)
    unusual_residues = dicts.with_default_value(0)
    inner_chain_residues_flagged_as_termini = []
    n_expected_atoms = 0
    unexpected_atoms = dicts.with_default_value(0)
    ignored_atoms = dicts.with_default_value(0)
    duplicate_atoms = dicts.with_default_value(0)
    classifications = dicts.with_default_value(0)
    modifications_used = dicts.with_default_value(0)
    incomplete_infos = dicts.with_default_value(0)
    link_ids = dicts.with_default_value(0)
    mm_pairs_not_linked = []
    n_unresolved_chain_links = 0
    n_chain_breaks = 0
    n_unresolved_chain_link_angles = 0
    n_unresolved_chain_link_dihedrals = 0
    n_unresolved_chain_link_chiralities = 0
    n_unresolved_chain_link_planarities = 0
    n_unresolved_chain_link_parallelities = 0
    corrupt_monomer_library_definitions = dicts.with_default_value(0)
    n_bond_proxies_already_assigned_to_first_conformer = 0
    n_unresolved_non_hydrogen_bonds = 0
    n_unresolved_non_hydrogen_angles = 0
    n_angles_discarded_because_of_special_positions = 0
    n_unresolved_non_hydrogen_dihedrals = 0
    n_dihedrals_discarded_because_of_special_positions = 0
    unsupported_chir_volume_sign = dicts.with_default_value(0)
    n_unresolved_non_hydrogen_chiralities = 0
    n_chiralities_discarded_because_of_special_positions = 0
    planarities_with_less_than_four_sites = dicts.with_default_value(0)
    n_unresolved_non_hydrogen_planarities = 0
    n_planarities_discarded_because_of_special_positions = 0
    mm = None
    prev_mm = None
    prev_prev_mm = None
    pdb_residues = conformer.residues()
    for i_residue,residue in enumerate(pdb_residues):
      def _get_next_residue():
        j = i_residue + 1
        if (j == len(pdb_residues)): return None
        return pdb_residues[j]
      # specific_residue_restraints
      specific_residue_restraints = None
      if apply_restraints_specifications:
        atoms = residue.atoms()
        alt_locs = {}
        for atom in atoms: alt_locs.setdefault(atom.parent().altloc, 0)
        residue_i_seqs=atoms.extract_i_seq()
        min_i_seq, max_i_seq = min(residue_i_seqs), max(residue_i_seqs)
        for selection, item in apply_restraints_specifications.items():
          #print min_i_seq,max_i_seq,list(selection)
          if min_i_seq in selection and max_i_seq in selection:
            if selection.all_eq(residue_i_seqs):
              print('%sResidue %s was targeted for' % (' '*8,
                                                               residue.id_str(),
                ), file=log)
              print('%srestraints from file: "%s"' % (' '*10,
                                                              item[1],
                ), file=log)
              if len(alt_locs)!=1:
                print('%sbut ignored because residue has complex alt. loc.' % (
                  ' '*10), file=log)
                continue
              specific_residue_restraints=item[1]
              apply_restraints_specifications[selection]="OK"
              break
      #
      mm = monomer_mapping(
        pdb_atoms=pdb_atoms,
        mon_lib_srv=mon_lib_srv,
        translate_cns_dna_rna_residue_names
          =translate_cns_dna_rna_residue_names,
        rna_sugar_pucker_analysis_params=rna_sugar_pucker_analysis_params,
        apply_cif_modifications=apply_cif_modifications,
        apply_cif_links_mm_pdbres_dict=apply_cif_links_mm_pdbres_dict,
        i_model=i_model,
        i_conformer=i_conformer,
        is_first_conformer_in_chain=is_first_conformer_in_chain,
        conf_altloc=conformer.altloc,
        pdb_residue=residue,
        next_pdb_residue=_get_next_residue(),
        chainid=residue.parent().parent().id,
        specific_residue_restraints=specific_residue_restraints)
      if mm.monomer and mm.monomer.cif_object:
        self._cif.chem_comps.append(mm.monomer.chem_comp)
        if specific_residue_restraints:
          self._cif.cif["comp_specific_%s" % residue.resname.strip()] = mm.monomer.cif_object
        else:
          self._cif.cif["comp_%s" % residue.resname.strip()] = mm.monomer.cif_object
      if (mm.monomer is None):
        def use_scattering_type_if_available_to_define_nonbonded_type():
          if (   residue.atoms_size() != 1
              or len(mm.active_atoms) != 1): return False
          atom = mm.active_atoms[0]
          ad_hoc = ad_hoc_single_atom_residue(
            residue_name=residue.resname,
            atom_name=atom.name,
            atom_element=atom.element)
          if (ad_hoc.scattering_type is None): return False
          entry = ener_lib.lib_atom.get(ad_hoc.energy_type, None)
          if (entry is None): return False
          i_seq = atom.i_seq
          scattering_type_registry.assign_directly(
            i_seq=i_seq, symbol=ad_hoc.scattering_type)
          nonbonded_energy_type_registry.assign_directly(
            i_seq=i_seq, symbol=ad_hoc.energy_type)
          ad_hoc_single_atom_residues[mm.residue_name] += 1
          return True
        if (not use_scattering_type_if_available_to_define_nonbonded_type()):
          unknown_residues[mm.residue_name] += 1
        n_chain_breaks += 1
      elif (prev_mm is not None and not residue.link_to_previous):
        n_chain_breaks += 1
      else:
        if (prev_mm is not None and prev_mm.monomer is not None):
          prev_mm.lib_link = get_lib_link(
            mon_lib_srv=mon_lib_srv,
            m_i=prev_mm,
            m_j=mm,
            )
          if (prev_mm.lib_link is None):
            link_ids[None] += 1
            mm_pairs_not_linked.append((prev_mm, mm))
          else:
            def within_linking_cutoff():
              result = True
              if(prev_mm.lib_link is None): return result
              for bond in prev_mm.lib_link.bond_list:
                atoms = (prev_mm.expected_atoms.get(bond.atom_id_1, None),
                         mm     .expected_atoms.get(bond.atom_id_2, None))
                # should not rely on sequence for ALL
                if None in atoms: return 999.
                i_seqs = [atom.i_seq for atom in atoms]
                s1 = sites_cart[i_seqs[0]]
                s2 = sites_cart[i_seqs[1]]
                import math
                link_distance = math.sqrt(
                  (s1[0]-s2[0])**2+(s1[1]-s2[1])**2+(s1[2]-s2[2])**2)
                if(link_distance > link_distance_cutoff): result = False
                return result
            mod_id = prev_mm.lib_link.chem_link.mod_id_1
            if (mod_id not in [None, ""] and within_linking_cutoff()):
              mod_mod_id = mon_lib_srv.mod_mod_id_dict[mod_id]
              prev_mm.apply_mod(mod_mod_id=mod_mod_id)
              prev_mm.resolve_unexpected()
            mod_id = prev_mm.lib_link.chem_link.mod_id_2
            if (mod_id not in [None, ""] and within_linking_cutoff()):
              mod_mod_id = mon_lib_srv.mod_mod_id_dict[mod_id]
              mm.apply_mod(mod_mod_id=mod_mod_id)
              mm.resolve_unexpected()
            link_resolution = add_bond_proxies(
              counters=counters(label="link_bond"),
              m_i=prev_mm,
              m_j=mm,
              bond_list=prev_mm.lib_link.bond_list,
              bond_simple_proxy_registry=geometry_proxy_registries.bond_simple,
              sites_cart=sites_cart,
              distance_cutoff=link_distance_cutoff)
            n_bond_proxies_already_assigned_to_first_conformer += \
              link_resolution.counters.already_assigned_to_first_conformer
            n_unresolved_chain_links \
              += link_resolution.counters.unresolved_non_hydrogen
            broken_bond_i_seq_pairs = link_resolution.broken_bond_i_seq_pairs
            n_chain_breaks += len(broken_bond_i_seq_pairs)
            link_resolution = add_angle_proxies(
              counters=counters(label="link_angle"),
              m_i=prev_mm,
              m_j=mm,
              angle_list=prev_mm.lib_link.angle_list,
              angle_proxy_registry=geometry_proxy_registries.angle,
              special_position_dict=special_position_dict,
              broken_bond_i_seq_pairs=broken_bond_i_seq_pairs)
            n_unresolved_chain_link_angles \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_dihedral_proxies(
              counters=counters(label="link_dihedral"),
              m_i=prev_mm,
              m_j=mm,
              tor_list=prev_mm.lib_link.tor_list,
              dihedral_function_type=dihedral_function_type,
              peptide_link_params=peptide_link_params,
              dihedral_proxy_registry=geometry_proxy_registries.dihedral,
              special_position_dict=special_position_dict,
              sites_cart=sites_cart,
              chem_link_id=prev_mm.lib_link.chem_link.id,
              broken_bond_i_seq_pairs=broken_bond_i_seq_pairs,
              cis_trans_specifications=cis_trans_specifications,
              )
            n_unresolved_chain_link_dihedrals \
              += link_resolution.counters.unresolved_non_hydrogen
            link_ids[link_resolution.chem_link_id] += 1
            link_resolution = add_chirality_proxies(
              counters=counters(label="link_chirality"),
              m_i=prev_mm,
              m_j=mm,
              chir_list=prev_mm.lib_link.chir_list,
              chirality_proxy_registry=geometry_proxy_registries.chirality,
              special_position_dict=special_position_dict,
              chir_volume_esd=chir_volume_esd,
              lib_link=prev_mm.lib_link,
              broken_bond_i_seq_pairs=broken_bond_i_seq_pairs)
            n_unresolved_chain_link_chiralities \
              += link_resolution.counters.unresolved_non_hydrogen
            def _add_planarity_proxies():
              link_resolution = add_planarity_proxies(
                counters=counters(label="link_planarity"),
                m_i=prev_mm,
                m_j=mm,
                plane_list=prev_mm.lib_link.get_planes(),
                planarity_proxy_registry=geometry_proxy_registries.planarity,
                special_position_dict=special_position_dict,
                broken_bond_i_seq_pairs=broken_bond_i_seq_pairs)
            _add_planarity_proxies()
            n_unresolved_chain_link_planarities \
              += link_resolution.counters.unresolved_non_hydrogen
            if peptide_link_params.apply_peptide_plane:
              prev_mm.lib_link = get_lib_link(
                mon_lib_srv=mon_lib_srv,
                m_i=prev_mm,
                m_j=mm,
                include_peptide_plane=True,
                )
              _add_planarity_proxies()
              link_ids["peptide plane"] += 1
              n_unresolved_chain_link_planarities \
                += link_resolution.counters.unresolved_non_hydrogen

      if (mm.monomer is not None):
        if (mm._is_unusual()):
          unusual_residues[mm.residue_name] += 1
        if (    mm.is_terminus == True
            and i_residue > 0
            and i_residue < conformer.residues_size()-1):
          inner_chain_residues_flagged_as_termini.append(
            residue_id_str(residue))
        n_expected_atoms += len(mm.expected_atoms)
        for atom_name in mm.unexpected_atoms.keys():
          unexpected_atoms[mm.residue_name+","+atom_name] += 1
        for atom_name,i_seqs in mm.ignored_atoms.items():
          ignored_atoms[mm.residue_name+","+atom_name] += len(i_seqs)
        for atom_name,i_seqs in mm.duplicate_atoms.items():
          duplicate_atoms[mm.residue_name+","+atom_name] += len(i_seqs)
        if (mm.incomplete_info is not None):
          incomplete_infos[mm.incomplete_info] += 1
        if (mm.monomer.classification is not None):
          classifications[mm.monomer.classification] += 1
        for chem_mod_id in mm.chem_mod_ids:
          modifications_used[chem_mod_id] += 1
        scattering_type_registry.assign_from_monomer_mapping(
          conf_altloc=conformer.altloc, mm=mm)
        nonbonded_energy_type_registry.assign_from_monomer_mapping(
          conf_altloc=conformer.altloc, mm=mm)
        mm.add_bond_proxies(
          bond_simple_proxy_registry=geometry_proxy_registries.bond_simple,
          use_neutron_distances= \
            restraints_loading_flags.get("use_neutron_distances", None))
        n_bond_proxies_already_assigned_to_first_conformer += \
          mm.bond_counters.already_assigned_to_first_conformer
        if (mm.bond_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.bond_counters.corrupt_monomer_library_definitions
        n_unresolved_non_hydrogen_bonds \
          += mm.bond_counters.unresolved_non_hydrogen
        mm.add_angle_proxies(
          special_position_dict=special_position_dict,
          angle_proxy_registry=geometry_proxy_registries.angle)
        if (mm.angle_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.angle_counters.corrupt_monomer_library_definitions
        n_unresolved_non_hydrogen_angles \
          += mm.angle_counters.unresolved_non_hydrogen
        n_angles_discarded_because_of_special_positions \
          += mm.angle_counters.discarded_because_of_special_positions
        mm.add_dihedral_proxies(
          dihedral_function_type=dihedral_function_type,
          special_position_dict=special_position_dict,
          dihedral_proxy_registry=geometry_proxy_registries.dihedral)
        if (mm.dihedral_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.dihedral_counters.corrupt_monomer_library_definitions
        n_unresolved_non_hydrogen_dihedrals \
          += mm.dihedral_counters.unresolved_non_hydrogen
        n_dihedrals_discarded_because_of_special_positions \
          += mm.dihedral_counters.discarded_because_of_special_positions
        mm.add_chirality_proxies(
          special_position_dict=special_position_dict,
          chirality_proxy_registry=geometry_proxy_registries.chirality,
          chir_volume_esd=chir_volume_esd)
        if (mm.chirality_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.chirality_counters.corrupt_monomer_library_definitions
        for s,n in mm.chirality_counters.unsupported_volume_sign.items():
          unsupported_chir_volume_sign[s] += n
        n_unresolved_non_hydrogen_chiralities \
          += mm.chirality_counters.unresolved_non_hydrogen
        n_chiralities_discarded_because_of_special_positions \
          += mm.chirality_counters.discarded_because_of_special_positions
        mm.add_planarity_proxies(
          special_position_dict=special_position_dict,
          planarity_proxy_registry=geometry_proxy_registries.planarity)
        if (mm.planarity_counters.corrupt_monomer_library_definitions > 0):
          corrupt_monomer_library_definitions[mm.residue_name] \
            += mm.planarity_counters.corrupt_monomer_library_definitions
        for p,n in mm.planarity_counters.less_than_four_sites.items():
          planarities_with_less_than_four_sites[mm.residue_name+":"+p] += n
        n_unresolved_non_hydrogen_planarities \
          += mm.planarity_counters.unresolved_non_hydrogen
        n_planarities_discarded_because_of_special_positions \
          += mm.planarity_counters.discarded_because_of_special_positions
        if (mm.monomer.chem_comp.id == "CYS"):
          sulphur_atom = mm.expected_atoms.get("SG", None)
          # XXX keep track of weights
          if (sulphur_atom is not None
              and sulphur_atom.i_seq not in cystein_sulphur_i_seqs):
            cystein_sulphur_i_seqs.append(sulphur_atom.i_seq)
            cystein_monomer_mappings.append(mm)
      if (conformation_dependent_restraints.is_available):
        cdr = conformation_dependent_restraints \
                .build_conformation_dependent_angle_proxies(
          angle_proxy_registry=geometry_proxy_registries.angle,
          dihedral_proxy_registry=geometry_proxy_registries.dihedral,
          monomer_mappings=(prev_prev_mm, prev_mm, mm),
          connectivity_i_j=True,
          connectivity_j_k=True,
          sites_cart=sites_cart)
        self.conformation_dependent_restraints_list.append(cdr)
      if (keep_monomer_mappings):
        all_monomer_mappings.append(mm)
      else:
        all_monomer_mappings.append(mm.summary())
      prev_prev_mm = prev_mm
      prev_mm = mm
      prev_mm.lib_link = None
    # ========================
    # End of residue loop for i_residue,residue in enumerate(pdb_residues):
    # ========================

    if (is_unique_model and log is not None):
      print("        Number of residues, atoms: %d, %d" % (
        conformer.residues_size(),
        n_expected_atoms + flex.sum(flex.long(list(unexpected_atoms.values())))), file=log)
      if (len(unknown_residues) > 0):
        print("          Unknown residues:", unknown_residues, file=log)
      if (len(ad_hoc_single_atom_residues) > 0):
        print("          Ad-hoc single atom residues:", \
          ad_hoc_single_atom_residues, file=log)
      if (len(unusual_residues) > 0):
        print("          Unusual residues:", unusual_residues, file=log)
      if (len(inner_chain_residues_flagged_as_termini) > 0):
        print("          Inner-chain residues flagged as termini:", \
          inner_chain_residues_flagged_as_termini, file=log)
      if (len(unexpected_atoms) > 0):
        print("          Unexpected atoms:", unexpected_atoms, file=log)
      if (len(ignored_atoms) > 0):
        print("          Ignored atoms:", ignored_atoms, file=log)
      if (len(duplicate_atoms) > 0):
        print("          Duplicate atoms:", duplicate_atoms, file=log)
      if (len(classifications) > 0):
        print("          Classifications:", classifications, file=log)
      if (len(modifications_used) > 0):
        print("          Modifications used:", modifications_used, file=log)
      if (len(incomplete_infos) > 0):
        print("          Incomplete info:", incomplete_infos, file=log)
    if (log is not None):
      if (len(link_ids) > 0):
        print("          Link IDs:", link_ids, file=log)
        if (len(link_ids) != 1):
          if (not_linked_show_max is None):
            show_max = len(mm_pairs_not_linked)
          else:
            show_max = not_linked_show_max
          n_not_shown = max(0, len(mm_pairs_not_linked) - show_max)
          if (n_not_shown == 1):
            show_max += 1
            n_not_shown = 0
          for pair in mm_pairs_not_linked[:show_max]:
            print("            Not linked:", file=log)
            for mm in pair:
              print("              %s" % residue_id_str(mm.pdb_residue), file=log)
          if (n_not_shown != 0):
            print("            ... (remaining %d not shown)" % n_not_shown, file=log)
    if (is_unique_model and log is not None):
      if (n_unresolved_chain_links > 0):
        print("          Unresolved chain links:", \
          n_unresolved_chain_links, file=log)
    if (log is not None):
      if (n_chain_breaks > 0):
        print("          Chain breaks:", n_chain_breaks, file=log)
    if (is_unique_model and log is not None):
      if (n_unresolved_chain_link_angles > 0):
        print("          Unresolved chain link angles:", \
          n_unresolved_chain_link_angles, file=log)
      if (n_unresolved_chain_link_dihedrals > 0):
        print("          Unresolved chain link dihedrals:", \
          n_unresolved_chain_link_dihedrals, file=log)
      if (n_unresolved_chain_link_chiralities > 0):
        print("          Unresolved chain link chiralities:", \
          n_unresolved_chain_link_chiralities, file=log)
      if (n_unresolved_chain_link_planarities > 0):
        print("          Unresolved chain link planarities:", \
          n_unresolved_chain_link_planarities, file=log)
      if (n_unresolved_chain_link_parallelities > 0):
        print("          Unresolved chain link parallelities:", \
          n_unresolved_chain_link_parallelities, file=log)
      if (len(corrupt_monomer_library_definitions) > 0):
        print("          Corrupt monomer library definitions:", \
          corrupt_monomer_library_definitions, file=log)
      if (n_unresolved_non_hydrogen_bonds > 0):
        print("          Unresolved non-hydrogen bonds:", \
          n_unresolved_non_hydrogen_bonds, file=log)
      if (n_unresolved_non_hydrogen_angles > 0):
        print("          Unresolved non-hydrogen angles:", \
          n_unresolved_non_hydrogen_angles, file=log)
    if (log is not None):
      if (n_angles_discarded_because_of_special_positions > 0):
        print("          Angles discarded because of special positions:", \
          n_angles_discarded_because_of_special_positions, file=log)
    if (is_unique_model and log is not None):
      if (n_unresolved_non_hydrogen_dihedrals > 0):
        print("          Unresolved non-hydrogen dihedrals:", \
          n_unresolved_non_hydrogen_dihedrals, file=log)
    if (log is not None):
      if (n_dihedrals_discarded_because_of_special_positions > 0):
        print("          Dihedrals discarded because of special positions:",\
          n_dihedrals_discarded_because_of_special_positions, file=log)
    if (is_unique_model and log is not None):
      if (len(unsupported_chir_volume_sign) > 0):
        print("          Unsupported chir.volume_sign:", \
          unsupported_chir_volume_sign, file=log)
      if (n_unresolved_non_hydrogen_chiralities > 0):
        print("          Unresolved non-hydrogen chiralities:", \
          n_unresolved_non_hydrogen_chiralities, file=log)
    if (log is not None):
      if (n_chiralities_discarded_because_of_special_positions > 0):
        print("          Chiralities discarded because of special positions:", \
          n_chiralities_discarded_because_of_special_positions, file=log)
    if (is_unique_model and log is not None):
      if (len(planarities_with_less_than_four_sites) > 0):
        print("          Planarities with less than four sites:", \
          planarities_with_less_than_four_sites, file=log)
      if (n_unresolved_non_hydrogen_planarities > 0):
        print("          Unresolved non-hydrogen planarities:", \
          n_unresolved_non_hydrogen_planarities, file=log)
    if (log is not None):
      if (n_planarities_discarded_because_of_special_positions > 0):
        print("          planarities discarded because of special positions:", \
          n_planarities_discarded_because_of_special_positions, file=log)
      if (n_bond_proxies_already_assigned_to_first_conformer > 0):
        print("          bond proxies already assigned to first conformer:", \
          n_bond_proxies_already_assigned_to_first_conformer, file=log)

class geometry_restraints_proxy_registries(object):

  def __init__(self, n_seq, strict_conflict_handling):
    self.bond_simple = geometry_restraints.bond_simple_proxy_registry(
      n_seq=n_seq, strict_conflict_handling=strict_conflict_handling)
    self.angle = geometry_restraints.angle_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)
    self.dihedral = geometry_restraints.dihedral_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)
    self.chirality = geometry_restraints.chirality_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)
    self.planarity = geometry_restraints.planarity_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)
    self.parallelity = geometry_restraints.parallelity_proxy_registry(
      strict_conflict_handling=strict_conflict_handling)

  # XXX TODO use counts to modify weights

  def initialize_tables(self):
    self.bond_simple.initialize_table()
    self.angle.initialize_table()
    self.dihedral.initialize_table()
    self.chirality.initialize_table()
    self.planarity.initialize_table()
    self.parallelity.initialize_table()

  def discard_tables(self):
    self.bond_simple.discard_table()
    self.angle.discard_table()
    self.dihedral.discard_table()
    self.chirality.discard_table()
    self.planarity.discard_table()
    self.parallelity.discard_table()

  def report(self, prefix, log):
    if (self.bond_simple.n_resolved_conflicts > 0):
      print(prefix + (
        "Number of resolved bond restraint conflicts: %d"
          % self.bond_simple.n_resolved_conflicts), file=log)
    if (self.angle.n_resolved_conflicts > 0):
      print(prefix + (
        "Number of resolved angle restraint conflicts: %d"
          % self.angle.n_resolved_conflicts), file=log)
    if (self.dihedral.n_resolved_conflicts > 0):
      print(prefix + (
        "Number of resolved dihedral restraint conflicts: %d"
          % self.dihedral.n_resolved_conflicts), file=log)
    if (self.chirality.n_resolved_conflicts > 0):
      print(prefix + (
        "Number of resolved chirality restraint conflicts: %d"
          % self.chirality.n_resolved_conflicts), file=log)
    if (self.planarity.n_resolved_conflicts > 0):
      print(prefix + (
        "Number of resolved planarity restraint conflicts: %d"
          % self.planarity.n_resolved_conflicts), file=log)
    if (self.parallelity.n_resolved_conflicts > 0):
      print(prefix + (
        "Number of resolved parallelity restraint conflicts: %d"
          % self.planarity.n_resolved_conflicts), file=log)

class cif_output_holder:
  def __init__(self):
    import iotbx.cif.model
    self.cif = iotbx.cif.model.cif()
    self.chem_comps = []

  def update(self, other):
    self.cif.update(other.cif)
    for cc in other.chem_comps:
      if cc not in self.chem_comps:
        self.chem_comps.append(cc)

from mmtbx.monomer_library.linking_mixins import linking_mixins

class build_all_chain_proxies(linking_mixins):

  def __init__(self,
        mon_lib_srv,
        ener_lib,
        params=None,
        pdb_inp=None,
        pdb_hierarchy=None,
        atom_selection_string=None,
        special_position_settings=None,
        crystal_symmetry=None,
        force_symmetry=False,
        substitute_non_crystallographic_unit_cell_if_necessary=False,
        strict_conflict_handling=True,
        keep_monomer_mappings=False,
        max_atoms=None,
        log=None,
        carbohydrate_callback=None,
        restraints_loading_flags=None,
                ):
    self._cif = cif_output_holder()
    # Proposal: use origin id in proxies.
    self.pdb_link_records = {}
    # END_MARKED_FOR_DELETION_OLEG
    if restraints_loading_flags is None:
      restraints_loading_flags = get_restraints_loading_flags(params)
    self.mon_lib_srv = mon_lib_srv
    assert special_position_settings is None or crystal_symmetry is None
    if (params is None): params = master_params.extract()
    self.params = params
    timer = user_plus_sys_time()
    self.time_building_chain_proxies = None
    self.pdb_inp = None
    self.pdb_hierarchy = None
    if (pdb_inp is not None):
      self.pdb_inp = pdb_inp
    if (pdb_hierarchy is not None):
      self.pdb_hierarchy = pdb_hierarchy
    else:
      self.pdb_hierarchy = self.pdb_inp.construct_hierarchy(
        sort_atoms=self.params.sort_atoms)
    #
    # optionally modify the model before processing
    #
    if 'all' in self.params.superpose_ideal_ligand:
      self.params.superpose_ideal_ligand = ideal_ligands
    from mmtbx.conformation_dependent_library import mcl
    for residue in self.params.superpose_ideal_ligand:
      if residue in [None, 'None']: continue
      info = mcl.superpose_ideal_residue_coordinates(self.pdb_hierarchy,
                                                     resname=residue,
                                                   )
      if info: print(info, file=log)
    if self.params.flip_symmetric_amino_acids:
      info = self.pdb_hierarchy.flip_symmetric_amino_acids()
      if info:
        print("\n  Symmetric amino acids flipped", file=log)
        print(info, file=log)
    if atom_selection_string is not None:
      sel = self.pdb_hierarchy.atom_selection_cache().selection(atom_selection_string)
      temp_string = self.pdb_hierarchy.select(sel).as_pdb_string()
      self.pdb_inp = pdb.input(source_info=None, lines=temp_string)
      self.pdb_hierarchy = self.pdb_inp.construct_hierarchy(sort_atoms=self.params.sort_atoms)
      if self.params.flip_symmetric_amino_acids:
        self.pdb_hierarchy.flip_symmetric_amino_acids()
    self.pdb_atoms = self.pdb_hierarchy.atoms()
    self.pdb_atoms.reset_i_seq()
    self.counts = self.pdb_hierarchy.overall_counts()
    self.counts.raise_residue_groups_with_multiple_resnames_using_same_altloc_if_necessary()
    self.counts.raise_improper_alt_conf_if_necessary()
    self.counts.raise_chains_with_mix_of_proper_and_improper_alt_conf_if_necessary()
    self.counts.raise_duplicate_atom_labels_if_necessary()
    if (log is not None):
      print("  Monomer Library directory:", file=log)
      print("   ", show_string(mon_lib_srv.root_path), file=log)
      print("  Total number of atoms:", self.pdb_atoms.size(), file=log)
    selection_cache = self.pdb_hierarchy.atom_selection_cache()
    # cis-trans specifications
    cis_trans_specifications = {}
    for cis_trans in self.params.apply_cis_trans_specification:
      if cis_trans.residue_selection is None: continue
      if cis_trans.residue_selection.lower().find("name ca")==-1:
        cis_trans.residue_selection+=" and name CA"
      t_selection = flex.size_t()
      t_selection = self.pdb_hierarchy.atom_selection_cache().selection(cis_trans.residue_selection).iselection()
      if len(t_selection)!=1:
        msg = """
    cis-trans specification selection "%s"
    produced %d atoms. Need to select one C-alpha atom.
    """ % (cis_trans.residue_selection, len(t_selection))
        print(msg, file=log)
        raise Sorry(msg)
      cis_trans_specifications[t_selection]=cis_trans.cis_trans_mod
    if cis_trans_specifications:
      print("  cis-trans peptide specifications", file=log)
      for cis_trans in self.params.apply_cis_trans_specification:
        print('    "%s" - %s' % ( cis_trans.residue_selection,
                                          cis_trans.cis_trans_mod.upper(),
          ), file=log)
    # apply a specific restraints file to a specific monomer
    apply_restraints_specifications = {}
    for acf in self.params.apply_cif_restraints:
      if acf.residue_selection is None: continue
      t_selection = flex.size_t()
      t_selection = self.pdb_hierarchy.atom_selection_cache().selection(acf.residue_selection).iselection()
      apply_restraints_specifications[t_selection]=[
        acf.residue_selection,
        acf.restraints_file_name,
        ]
    if apply_restraints_specifications:
      print("  Apply specific restraints filenames to specific monomers", file=log)
      for acf in self.params.apply_cif_restraints:
        print('    "%s" - %s' % (acf.residue_selection,
                                         acf.restraints_file_name,
          ), file=log)
    self.special_position_settings = None
    self._site_symmetry_table = None
    self.sites_cart = None
    self._sites_cart_exact = None
    self.use_cdl = None
    if (max_atoms is not None
        and self.pdb_atoms.size() > max_atoms):
      if (log is not None):
        print("  More than %d atoms: no processing." % max_atoms, file=log)
        return
    self.sites_cart = self.pdb_atoms.extract_xyz()
    models = self.pdb_hierarchy.models()
    if (log is not None):
      print("  Number of models:", len(models), file=log)
    n_seq = self.pdb_atoms.size()
    def set_model_indices():
      self.model_indices = flex.size_t(n_seq, n_seq)
      for i_model,model in enumerate(models):
        self.model_indices.set_selected(model.atoms().extract_i_seq(), i_model)
      assert self.model_indices.count(n_seq) == 0
    set_model_indices()
    altloc_i_conformer = {}
    def set_conformer_indices():
      self.conformer_indices = flex.size_t(n_seq, 0)
      altloc_indices = self.pdb_hierarchy.altloc_indices()
      if ("" in altloc_indices): p = 0
      else:                      p = 1
      altlocs = sorted(altloc_indices.keys())
      for i,altloc in enumerate(altlocs):
        if (altloc == ""): continue
        self.conformer_indices.set_selected(altloc_indices[altloc], i+p)
        altloc_i_conformer[altloc] = i+p
      altloc_i_conformer[""] = 0
    set_conformer_indices()
    sym_excl_residue_groups = []
    def set_sym_excl_indices():
      self.sym_excl_indices = flex.size_t(n_seq, 0)
      for rg in self.pdb_hierarchy.residue_groups():
        rg_atoms = rg.atoms()
        if (rg_atoms.extract_occ().all_lt(1.0)):
          sym_excl_residue_groups.append(rg)
          self.sym_excl_indices.set_selected(
            rg_atoms.extract_i_seq(), len(sym_excl_residue_groups))
    set_sym_excl_indices()
    def set_donor_acceptor_excl_groups():
      self.donor_acceptor_excl_groups = flex.size_t(n_seq, 0)
      counter = 0
      for rg in self.pdb_hierarchy.residue_groups():
        rg_atoms = rg.atoms()
        self.donor_acceptor_excl_groups.set_selected(
          rg_atoms.extract_i_seq(), counter)
        counter += 1
    set_donor_acceptor_excl_groups()
    if (    special_position_settings is None
        and crystal_symmetry is not None):
      special_position_settings = crystal_symmetry.special_position_settings(
        min_distance_sym_equiv=params.min_distance_sym_equiv)
    if(self.pdb_inp is None):
      self.special_position_settings = \
        iotbx.pdb.construct_special_position_settings(
          crystal_symmetry          = crystal_symmetry,
          special_position_settings = special_position_settings,
          weak_symmetry             = not force_symmetry,
          min_distance_sym_equiv    = params.min_distance_sym_equiv)
      occ = self.pdb_hierarchy.atoms().extract_occ()
    else:
      self.special_position_settings = self.pdb_inp.special_position_settings(
        special_position_settings=special_position_settings,
        min_distance_sym_equiv=params.min_distance_sym_equiv,
        weak_symmetry=not force_symmetry)
      occ = self.pdb_inp.atoms().extract_occ()
    if(self.special_position_settings is not None and
       self.special_position_settings.unit_cell() is not None and
       self.special_position_settings.unit_cell().volume() <
       flex.sum(occ)*5 and
       not self.params.disable_uc_volume_vs_n_atoms_check):
      msg = """Unit cell volume is incompatible with number of atoms.
  Unit cell parameters: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f
  Check CRYST1 record or other sources of crystal symmetry.
"""
      raise Sorry(msg%self.special_position_settings.unit_cell().parameters())
    if (self.special_position_settings is not None
        and (   self.special_position_settings.unit_cell() is None
             or self.special_position_settings.space_group_info() is None)):
      self.special_position_settings = None
    if (    self.special_position_settings is None
        and substitute_non_crystallographic_unit_cell_if_necessary):
      self.special_position_settings = crystal.non_crystallographic_symmetry(
        sites_cart=self.sites_cart).special_position_settings(
          min_distance_sym_equiv=params.min_distance_sym_equiv)
    if (self.special_position_settings is None):
      self.special_position_indices = None
    else:
      self.special_position_indices = \
        self.site_symmetry_table().special_position_indices()

    self.special_position_dict = special_position_dict(
        special_position_indices=self.special_position_indices)
    self.process_apply_cif_modification(mon_lib_srv=mon_lib_srv, log=log)
    self.process_apply_cif_link(mon_lib_srv=mon_lib_srv, log=log)
    self.conformation_dependent_restraints_list = []
    model_type_indices = [-1] * len(models)
    self.all_monomer_mappings = []
    self.scattering_type_registry = scattering_type_registry(
      # XXX should be same as in pdb_inp.xray_structure_simple
      scattering_types=self.pdb_atoms.extract_element(strip=True),
      strict_conflict_handling=strict_conflict_handling)
    self.nonbonded_energy_type_registry = nonbonded_energy_type_registry(
      n_seq=self.pdb_atoms.size(),
      strict_conflict_handling=strict_conflict_handling)
    self.geometry_proxy_registries = geometry_restraints_proxy_registries(
      n_seq=self.pdb_atoms.size(),
      strict_conflict_handling=strict_conflict_handling)
    self.cystein_sulphur_i_seqs = flex.size_t()
    self.cystein_monomer_mappings = []
    n_unique_models = 0
    for i_model,model in enumerate(models):
      if (log is not None):
        print('  Model: "%s"' % model.id, file=log)
      is_unique_model = not is_same_model_as_before(
        model_type_indices, i_model, models)
      if (is_unique_model):
        n_unique_models += 1
      elif (log is not None):
        print("    Same as model", \
          models[model_type_indices[i_model]].id, file=log)
      if (is_unique_model and log is not None):
        print("    Number of chains:", model.chains_size(), file=log)
      self.geometry_proxy_registries.initialize_tables()
      apply_cif_links_mm_pdbres_dict = dict(
        self.empty_apply_cif_links_mm_pdbres_dict)
      for chain in model.chains():
        conformers = chain.conformers()
        if (is_unique_model and log is not None):
          print('    Chain: "%s"' % chain.id, file=log)
          print("      Number of atoms:", chain.atoms_size(), file=log)
          print("      Number of conformers:", len(conformers), file=log)
          flush_log(log)
        for j_conformer,conformer in enumerate(conformers):
          if (is_unique_model and log is not None):
            print('      Conformer: "%s"' % conformer.altloc, file=log)
            flush_log(log)
          i_conformer = altloc_i_conformer[conformer.altloc]
          chain_proxies = build_chain_proxies(
            mon_lib_srv=mon_lib_srv,
            ener_lib=ener_lib,
            translate_cns_dna_rna_residue_names
              =self.params.translate_cns_dna_rna_residue_names,
            rna_sugar_pucker_analysis_params
              =self.params.rna_sugar_pucker_analysis,
            apply_cif_modifications=self.apply_cif_modifications,
            apply_cif_links_mm_pdbres_dict=apply_cif_links_mm_pdbres_dict,
            link_distance_cutoff=self.params.link_distance_cutoff,
            not_linked_show_max=self.params.show_max_items.not_linked,
            dihedral_function_type=self.params.dihedral_function_type,
            chir_volume_esd=self.params.chir_volume_esd,
            peptide_link_params=self.params.peptide_link,
            pdb_hierarchy=self.pdb_hierarchy,
            pdb_atoms=self.pdb_atoms,
            sites_cart=self.sites_cart,
            special_position_dict=self.special_position_dict,
            keep_monomer_mappings=keep_monomer_mappings,
            all_monomer_mappings=self.all_monomer_mappings,
            scattering_type_registry=self.scattering_type_registry,
            nonbonded_energy_type_registry=self.nonbonded_energy_type_registry,
            geometry_proxy_registries=self.geometry_proxy_registries,
            cystein_sulphur_i_seqs=self.cystein_sulphur_i_seqs,
            cystein_monomer_mappings=self.cystein_monomer_mappings,
            is_unique_model=is_unique_model,
            i_model=i_model,
            i_conformer=i_conformer,
            is_first_conformer_in_chain=(j_conformer == 0),
            conformer=conformer,
            conformation_dependent_restraints_list=
              self.conformation_dependent_restraints_list,
            restraints_loading_flags=restraints_loading_flags,
            fatal_problem_max_lines
              =self.params.show_max_items.fatal_problem_max_lines,
            cis_trans_specifications=cis_trans_specifications,
            apply_restraints_specifications=apply_restraints_specifications,
            log=log,
            )
          self._cif.update(chain_proxies._cif)
          assert not chain_proxies.pdb_link_records
          self.conformation_dependent_restraints_list = \
            chain_proxies.conformation_dependent_restraints_list
          del chain_proxies
          flush_log(log)
      if apply_restraints_specifications:
        for selection, item in apply_restraints_specifications.items():
          if item=="OK": continue
          print("%sRestraints for '%s'" % (' '*6,
                                                   item[0],
            ), file=log)
          print('%swere not modified by "%s"' % (' '*8,
                                                         item[1],
            ), file=log)
      #
      # Identify disulfide bond exclusions BEGIN
      self.disulfide_bond_exclusions_selection = flex.size_t()
      if(self.cystein_sulphur_i_seqs.size()>0):
        if(params.disulfide_bond_exclusions_selection_string is not None):
          self.disulfide_bond_exclusions_selection = \
            self.pdb_hierarchy.atom_selection_cache().selection(
              params.disulfide_bond_exclusions_selection_string).iselection()
        else:
          exclusion_list = ["H","D","T","S","O","P","N","C","SE"]
          cystein_sulphur_atoms_exclude = []
          for atom in self.pdb_atoms:
            e_ = atom.element.strip().upper()
            e = atom.determine_chemical_element_simple()
            if(e is not None):
              e = e.strip().upper()
              if(e_.strip() != ""):
                if(e_!=e):
                  raise Sorry("Bad element type: '%s', '%s'."%(e,e_))
              if(not e in exclusion_list):
                for cs_i_seq in self.cystein_sulphur_i_seqs:
                  csa = self.pdb_atoms[cs_i_seq]
                  if(csa.distance(atom) < params.exclusion_distance_cutoff):
                    cystein_sulphur_atoms_exclude.append(csa)
                    self.disulfide_bond_exclusions_selection.append(csa.i_seq)
      if(self.disulfide_bond_exclusions_selection.size()>0):
        if log is not None:
          print(file=log)
          print("List of CYS excluded from plausible disulfide bonds:", file=log)
          print("  (reason: may participate in coordination)", file=log)
        for i_seq in self.disulfide_bond_exclusions_selection:
          a = self.pdb_atoms[i_seq]
          if log is not None: print("  %s"%a.format_atom_record(), file=log)
          dces = a.determine_chemical_element_simple()
          if dces is None:
            raise Sorry("Atom '%s' has unknown chemical element symbol" % a.format_atom_record())
          e = dces.strip().upper()
          if(e!="S"):
            raise Sorry("disulfide_bond_exclusions_selection_string must select CYS sulfur.")
        if log is not None: print(file=log)
      tmp = flex.size_t()
      for i_seq in self.cystein_sulphur_i_seqs:
        if(not i_seq in self.disulfide_bond_exclusions_selection):
          tmp.append(i_seq)
      self.cystein_sulphur_i_seqs = tmp
      # Identify disulfide bond exclusions END
      #
      n_unresolved_apply_cif_link_bonds = 0
      n_unresolved_apply_cif_link_angles = 0
      n_unresolved_apply_cif_link_dihedrals = 0
      n_unresolved_apply_cif_link_chiralities = 0
      n_unresolved_apply_cif_link_planarities = 0
      for apply in self.apply_cif_links:
        if (apply.was_used):
          if(apply.automatic):
            print('  Automatic links duplication of user input', file=log)
          continue
        mms = []
        for pdbres in apply.pdbres_pair:
          mms.append(apply_cif_links_mm_pdbres_dict[pdbres])
        for m_i_list in mms[0].values():
          assert len(m_i_list) == 1
          m_i = m_i_list[0]
          def raise_missing_cif(i_pair):
            raise Sorry(
              "Error processing apply_cif_link:\n"
              + "  data_link: %s\n" % show_string(apply.data_link)
              + "  Missing CIF file for residue: %s"
                  % apply.pdbres_pair[i_pair])
          if (m_i.monomer is None):
            raise_missing_cif(i_pair=0)
          for m_j_list in mms[1].values():
            assert len(m_j_list) == 1
            m_j = m_j_list[0]
            if (m_j.monomer is None):
              raise_missing_cif(i_pair=1)
            if (    m_i.i_conformer != 0
                and m_j.i_conformer != 0
                and m_i.i_conformer != m_j.i_conformer):
              continue
            apply.was_used = True
            def raise_if_corrupt(link_resolution):
              counters = link_resolution.counters
              if (counters.corrupt_monomer_library_definitions != 0):
                raise Sorry(
                  "Error processing %s:\n" % counters.label
                  + "  data_link: %s\n" % show_string(apply.data_link)
                  + "  residue 1: %s\n" % apply.pdbres_pair[0]
                  + "  residue 2: %s\n" % apply.pdbres_pair[1]
                  + "  Possible problems giving rise to this error include:\n"
                  + "    - residue selections in apply_cif_link block swapped\n"
                  + "    - corrupt CIF link definition\n"
                  + "    - corrupt or missing CIF modifications associated with"
                      " this link\n"
                  + "  If none of this applies, send email to:\n"
                  + "    bugs@phenix-online.org")
            # automatic link creation
            if not mon_lib_srv.link_link_id_dict.get(apply.data_link, False):
              if getattr(apply, "possible_peptide_link", False):
                link = self.create_link(apply, m_i, m_j)
                mon_lib_srv.link_link_id_dict[apply.data_link] = link
              elif getattr(apply, "possible_rna_dna_link", False):
                link = self.create_link(apply, m_i, m_j)
                mon_lib_srv.link_link_id_dict[apply.data_link] = link
              else:
                link = self.create_link(apply, m_i, m_j)
                mon_lib_srv.link_link_id_dict[apply.data_link] = link
                # perform this in process_custom_links
              continue
            link = mon_lib_srv.link_link_id_dict[apply.data_link]
            if hasattr(apply, "atom1") and hasattr(apply, "atom2"):
              i_seqs = [apply.atom1.i_seq, apply.atom2.i_seq]
              if self.geometry_proxy_registries.bond_simple.is_proxy_set(
                  i_seqs=i_seqs,
                  ):
                if apply.automatic:
                  print('%sDuplicate links ignored : %s' % (
                    ' '*6,
                    apply.data_link,
                    ), file=log)
                  continue
            link_resolution = add_bond_proxies(
              counters=counters(label="apply_cif_link_bond"),
              m_i=m_i,
              m_j=m_j,
              bond_list=link.bond_list,
              bond_simple_proxy_registry=self.geometry_proxy_registries
                .bond_simple,
              sites_cart=self.sites_cart,
              distance_cutoff=self.params.link_distance_cutoff)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_bonds \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_angle_proxies(
              counters=counters(label="apply_cif_link_angle"),
              m_i=m_i,
              m_j=m_j,
              angle_list=link.angle_list,
              angle_proxy_registry=self.geometry_proxy_registries.angle,
              special_position_dict=self.special_position_dict)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_angles \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_dihedral_proxies(
              counters=counters(label="apply_cif_link_dihedral"),
              m_i=m_i,
              m_j=m_j,
              tor_list=link.tor_list,
              dihedral_function_type=self.params.dihedral_function_type,
              peptide_link_params=self.params.peptide_link,
              dihedral_proxy_registry=self.geometry_proxy_registries.dihedral,
              special_position_dict=self.special_position_dict,
              sites_cart=self.sites_cart,
              chem_link_id=link.chem_link.id)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_dihedrals \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_chirality_proxies(
              counters=counters(label="apply_cif_link_chirality"),
              m_i=m_i,
              m_j=m_j,
              chir_list=link.chir_list,
              chirality_proxy_registry=self.geometry_proxy_registries.chirality,
              special_position_dict=self.special_position_dict,
              chir_volume_esd=self.params.chir_volume_esd,
              lib_link=link)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_chiralities \
              += link_resolution.counters.unresolved_non_hydrogen
            link_resolution = add_planarity_proxies(
              counters=counters(label="apply_cif_link_planarity"),
              m_i=m_i,
              m_j=m_j,
              plane_list=link.get_planes(),
              planarity_proxy_registry=self.geometry_proxy_registries.planarity,
              special_position_dict=self.special_position_dict)
            raise_if_corrupt(link_resolution)
            n_unresolved_apply_cif_link_planarities \
              += link_resolution.counters.unresolved_non_hydrogen
      if (log is not None):
        if (n_unresolved_apply_cif_link_bonds > 0):
          print("          Unresolved apply_cif_link bonds:", \
            n_unresolved_apply_cif_link_bonds, file=log)
        if (n_unresolved_apply_cif_link_angles > 0):
          print("          Unresolved apply_cif_link angles:", \
            n_unresolved_apply_cif_link_angles, file=log)
        if (n_unresolved_apply_cif_link_dihedrals > 0):
          print("          Unresolved apply_cif_link dihedrals:", \
            n_unresolved_apply_cif_link_dihedrals, file=log)
        if (n_unresolved_apply_cif_link_chiralities > 0):
          print("          Unresolved apply_cif_link chiralities:", \
            n_unresolved_apply_cif_link_chiralities, file=log)
        if (n_unresolved_apply_cif_link_planarities > 0):
          print("          Unresolved apply_cif_link planarities:", \
            n_unresolved_apply_cif_link_planarities, file=log)
        flush_log(log)
    for apply in self.apply_cif_links:
      if (not apply.was_used):
        raise RuntimeError(
          "Unused apply_cif_link: %s %s" % (
            apply.data_link, str(apply.pdbres_pair)))
    #
    if carbohydrate_callback:
      print('\n  Calling carbohydrate callback')
      if not hasattr(carbohydrate_callback, "pdb_interpretation_callback"):
        print('    No PDB interpretation callback found, skipped')
      else:
        carbohydrate_callback.pdb_interpretation_callback(self)
    #
    # self.geometry_proxy_registries.discard_tables()
    # self.scattering_type_registry.discard_tables()
    # self.nonbonded_energy_type_registry.discard_tables()
    if (log is not None):
      if (n_unique_models != 1):
        print("  Number of unique models:", n_unique_models, file=log)
      if (len(sym_excl_residue_groups) != 0):
        print("  Residues with excluded nonbonded symmetry interactions:", \
          len(sym_excl_residue_groups), file=log)
        show_residue_groups(
          residue_groups=sym_excl_residue_groups,
          log=log,
          prefix="    ",
          max_items=self.params.show_max_items
            .residues_with_excluded_nonbonded_symmetry_interactions)
      self.geometry_proxy_registries.report(log=log, prefix="  ")
      self.fatal_problems_report(
        prefix="  ",
        log=log,
        max_lines=self.params.show_max_items.fatal_problem_max_lines,
        )
    self.process_custom_nonbonded_symmetry_exclusions(
      log=log,
      curr_sym_excl_index=len(sym_excl_residue_groups))
    self.time_building_chain_proxies = timer.elapsed()
    # Make sure pdb_hierarchy and xray_structure are consistent
    if(self.special_position_settings is not None):
      self.pdb_hierarchy.adopt_xray_structure(self.extract_xray_structure())

  def __getstate__(self):
    indexer = dict( ( a, i) for ( i, a ) in enumerate( self.pdb_hierarchy.atoms() ) )
    # pdb_atoms is an af_shared_atom array which has no pickling implemented
    # so convert it to a list of indices of the atoms in pdb_hierarchy
    import iotbx_pdb_hierarchy_ext
    if type(self.pdb_atoms) == iotbx_pdb_hierarchy_ext.af_shared_atom:
      self.pdb_atoms = [indexer[ a ] for a in self.pdb_atoms]
    sdict = self.__dict__.copy()
    return sdict

  def __setstate__(self, state):
    hroot = state["pdb_hierarchy"]
    # Assuming pdb_atoms was stored as a list of indices of the atoms in
    # pdb_hierarchy convert these back to an af_shared_atom array
    state["pdb_atoms"] = [ hroot.atoms()[i] for i in state["pdb_atoms"] ]
    self.__dict__.update( state )

  def update_internals_due_to_coordinates_change(self, pdb_h):
    self.pdb_hierarchy = pdb_h
    self.pdb_atoms = self.pdb_hierarchy.atoms()
    self.sites_cart = self.pdb_atoms.extract_xyz()
    self._sites_cart_exact = None

  def fatal_problems_report(self, prefix="", log=None, max_lines=10):
    self.scattering_type_registry.report(
      pdb_atoms=self.pdb_atoms, log=log, prefix=prefix, max_lines=max_lines)
    self.nonbonded_energy_type_registry.report(
      pdb_atoms=self.pdb_atoms, log=log, prefix=prefix, max_lines=max_lines)

  def fatal_problems_message(self,
        ignore_unknown_scattering_types=False,
        ignore_unknown_nonbonded_energy_types=False):
    result = ["Fatal problems interpreting model file:"]
    for reg,ignore in [
          (self.scattering_type_registry,
            ignore_unknown_scattering_types),
          (self.nonbonded_energy_type_registry,
             ignore_unknown_nonbonded_energy_types)]:
      if (not ignore):
        n_unknown = reg.n_unknown_type_symbols()
        if (n_unknown != 0):
          result.append("%s: %d" % (reg.report_unknown_message(), n_unknown))
    if (len(result) == 1): return None
    result.extend([
      "  Please edit the model file to resolve the problems and/or supply a",
      "  CIF file with matching restraint definitions, along with",
      "  apply_cif_modification and apply_cif_link parameter definitions",
      "  if necessary."])
    if (self.scattering_type_registry.n_unknown_type_symbols() != 0):
      result.extend([
        "  It is best practice to define the element names in",
        "  columns 77-78 of the PDB file."])
    return "\n  ".join(result)

  def extract_secondary_structure(self, log=None):
    if hasattr(self.pdb_inp, "extract_secondary_structure"):
      return self.pdb_inp.extract_secondary_structure(log=log)
    else:
      return None

  def site_symmetry_table(self):
    if (self._site_symmetry_table is None):
      assert self.special_position_settings is not None
      assert self.sites_cart is not None
      self._site_symmetry_table = \
        self.special_position_settings.site_symmetry_table(
          sites_cart=self.sites_cart,
          unconditional_general_position_flags=(
            self.pdb_atoms.extract_occ() != 1))
    return self._site_symmetry_table

  def sites_cart_exact(self):
    if (self._sites_cart_exact is None):
      self._sites_cart_exact = self.site_symmetry_table().apply_symmetry_sites(
        unit_cell=self.special_position_settings.unit_cell(),
        sites_cart=self.sites_cart)
    return self._sites_cart_exact

  def sel_classification(self, classification):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification == classification):
        result.set_selected(summary.all_associated_i_seqs(), True)
    return result

  def sel_backbone_or_sidechain(self, backbone_flag, sidechain_flag):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification == "peptide"):
        for atom in summary.expected_atoms:
          # XXX hydrogens not included
          if (atom.name.strip() in ["N", "CA", "C", "O"]):
            result[atom.i_seq] = backbone_flag
          else:
            result[atom.i_seq] = sidechain_flag
      elif (summary.classification in ["RNA", "DNA"]):
        for atom in summary.expected_atoms:
          # XXX hydrogens not included
          if (atom.name.strip()
                in ["P", "O1P", "O2P", "O3'", "O5'",
                         "OP1", "OP2", "O3*", "O5*",
                                       "O2*", "O2'",
                    "O4'", "C1'", "C2'", "C3'", "C4'", "C5'",
                    "O4*", "C1*", "C2*", "C3*", "C4*", "C5*"]):
            result[atom.i_seq] = backbone_flag
          else:
            result[atom.i_seq] = sidechain_flag
    return result

  def sel_backbone(self):
    return self.sel_backbone_or_sidechain(True, False)

  def sel_sidechain(self):
    return self.sel_backbone_or_sidechain(False, True)

  def sel_phosphate(self):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification in ["RNA", "DNA"]):
        for atom in summary.expected_atoms:
          if (atom.name.strip()
                in ["P", "O1P", "O2P", "O3'", "O5'",
                         "OP1", "OP2", "O3*", "O5*"]):
            result[atom.i_seq] = True
    return result

  def sel_ribose(self):
    result = flex.bool(self.pdb_atoms.size(), False)
    for summary in self.all_monomer_mappings:
      if (summary.classification in ["RNA", "DNA"]):
        for atom in summary.expected_atoms:
          if (atom.name.strip()
                in ["O4'", "C1'", "C2'", "C3'", "C4'", "C5'", "O2'",
                    "O4*", "C1*", "C2*", "C3*", "C4*", "C5*", "O2*"]):
            result[atom.i_seq] = True
    return result

  def sel_within(self, radius, primary_selection):
    assert radius > 0
    assert self.special_position_settings is not None
    return crystal.neighbors_fast_pair_generator(
      asu_mappings=self.special_position_settings.asu_mappings(
        buffer_thickness=radius,
        sites_cart=self.sites_cart),
      distance_cutoff=radius).neighbors_of(
        primary_selection=primary_selection)

  def _selection_callback(self, word, word_iterator, result_stack):
    lword = word.value.lower()
    if (lword in ["peptide", "protein"]):
      result_stack.append(self.sel_classification(classification="peptide"))
    elif (lword == "rna"):
      result_stack.append(self.sel_classification(classification="RNA"))
    elif (lword == "dna"):
      result_stack.append(self.sel_classification(classification="DNA"))
    elif (lword == "water"):
      result_stack.append(self.sel_classification(classification="water"))
    elif (lword in ["nucleotide", "nuc"]):
      result_stack.append(
          self.sel_classification(classification="RNA")
        | self.sel_classification(classification="DNA"))
    elif (lword == "backbone"):
      result_stack.append(self.sel_backbone())
    elif (lword == "sidechain"):
      result_stack.append(self.sel_sidechain())
    elif (lword == "phosphate"):
      result_stack.append(self.sel_phosphate())
    elif (lword == "ribose"):
      result_stack.append(self.sel_ribose())
    elif (lword == "within"):
      assert word_iterator.pop().value == "("
      radius = float(word_iterator.pop().value)
      assert word_iterator.pop().value == ","
      sel = self.pdb_hierarchy.atom_selection_cache().selection_parser(
        word_iterator=word_iterator,
        callback=self._selection_callback,
        expect_nonmatching_closing_parenthesis=True)
      result_stack.append(self.sel_within(radius=radius,primary_selection=sel))
    else:
      return False
    return True

  def selection(self, string, cache=None, optional=True):
    if (cache is None): cache = self.pdb_hierarchy.atom_selection_cache()
    return cache.selection(
      string=string,
      optional=optional,
      callback=self._selection_callback)

  def iselection(self, string, cache=None, optional=True):
    result = self.selection(string=string, cache=cache, optional=optional)
    if (result is None):
      return None
    return result.iselection()

  def process_apply_cif_modification(self, mon_lib_srv, log):
    self.apply_cif_modifications = {}
    atoms = self.pdb_atoms
    sel_cache = None
    for apply in self.params.apply_cif_modification:
      if (apply.data_mod is None): continue
      print("  apply_cif_modification:", file=log)
      print("    data_mod:", apply.data_mod, file=log)
      mod = mon_lib_srv.mod_mod_id_dict.get(apply.data_mod)
      if (mod is None):
        print(file=log)
        raise Sorry(
          "Missing CIF modification: data_mod_%s\n" % apply.data_mod
          + "  Please check for spelling errors or specify the file name\n"
          + "  with the modification as an additional argument.")
      print("    residue_selection:", apply.residue_selection, file=log)
      if (sel_cache is None):
        sel_cache = self.pdb_hierarchy.atom_selection_cache()
      iselection = self.phil_atom_selection(
        cache=sel_cache,
        scope_extract=apply,
        attr="residue_selection").iselection()
      pdbres_set = set()
      for i_seq in iselection:
        pdbres_set.add(atoms[i_seq].id_str(pdbres=True))
      for pdbres in pdbres_set:
        self.apply_cif_modifications.setdefault(pdbres, []).append(
          apply.data_mod)

  def process_apply_cif_link(self, mon_lib_srv, log):
    self.apply_cif_links = []
    self.empty_apply_cif_links_mm_pdbres_dict = {}
    atoms = self.pdb_atoms
    sel_cache = None
    for apply in self.params.apply_cif_link:
      if (apply.data_link is None): continue
      print("  apply_cif_link:", file=log)
      print("    data_link:", apply.data_link, file=log)
      link = mon_lib_srv.link_link_id_dict.get(apply.data_link)
      if (link is None):
        print(file=log)
        raise Sorry(
          "Missing CIF link: data_link_%s\n" % apply.data_link
          + "  Please check for spelling errors or specify the file name\n"
          + "  with the link as an additional argument.")
      mod_ids = []
      for mod_attr in ["mod_id_1", "mod_id_2"]:
        mod_id = getattr(link.chem_link, mod_attr)
        if (mod_id == ""): mod_id = None
        mod_ids.append(mod_id)
        if (mod_id is not None):
          print("      %s:" % mod_attr, mod_id, file=log)
          mod = mon_lib_srv.mod_mod_id_dict.get(mod_id)
          if (mod is None):
            print(file=log)
            raise Sorry(
              "Missing CIF modification: data_mod_%s\n" % mod_id
              + "  Please check for spelling errors or specify the file name\n"
              + "  with the modification as an additional argument.")
      sel_attrs = ["residue_selection_"+n for n in ["1", "2"]]
      pdbres_pair = []
      for attr in sel_attrs:
        print("    %s:" % attr, getattr(apply, attr), file=log)
        if (sel_cache is None):
          sel_cache = self.pdb_hierarchy.atom_selection_cache()
        iselection = self.phil_atom_selection(
          cache=sel_cache,
          scope_extract=apply,
          attr=attr).iselection()
        pdbres_set = set()
        pdbres_set_no_segid = set()
        for i_seq in iselection:
          atom = atoms[i_seq]
          pdbres_set.add(atom.id_str(pdbres=True))
          pdbres_set_no_segid.add(atom.id_str(pdbres=True,suppress_segid=True))
        if (len(pdbres_set) != 1):
          # XXX models?
          if (len(pdbres_set_no_segid) != 1):
            raise Sorry("Not exactly one residue selected.")
          raise Sorry(
            "Selected residue has multiple segid. This is not supported.")
        pdbres_pair.append(list(pdbres_set)[0])
      for pdbres,mod_id in zip(pdbres_pair, mod_ids):
        if (mod_id is not None):
          self.apply_cif_modifications.setdefault(pdbres, []).append(mod_id)
      self.apply_cif_links.append(group_args(
        pdbres_pair=pdbres_pair,
        data_link=apply.data_link,
        was_used=False))
      for pdbres in pdbres_pair:
        self.empty_apply_cif_links_mm_pdbres_dict[pdbres] = {}

  def create_link(self,
                  apply,
                  m_i,
                  m_j,
                  order=1,
                  angles=True,
                  verbose=False,
                  ):
    assert 0
    from mmtbx.monomer_library import linking_utils
    from math import sqrt
    from mmtbx.monomer_library.cif_types import link_link_id, chem_comp
    from mmtbx.monomer_library.cif_types import chem_link_bond, chem_link_angle
    from mmtbx.monomer_library.bondlength_defaults import get_default_bondlength
    link = link_link_id("custom",
                        chem_comp(),
                        )
    # bond
    bond = chem_link_bond()
    bond.atom_1_comp_id = 1
    bond.atom_id_1 = apply.atom1.name.strip()
    bond.atom_2_comp_id = 2
    bond.atom_id_2 = apply.atom2.name.strip()
    if order==1:
      bond.type = "single"
      bond.value_dist = get_default_bondlength(apply.atom1.element,
                                               apply.atom2.element,
                                               )
      if bond.value_dist is None:
        bond.value_dist = sqrt(linking_utils.get_distance2(apply.atom1,
                                                           apply.atom2,
                                                           ))
        if verbose: print("bond will be maintained")
    bond.value_dist_esd = 0.02
    if verbose:
      print('Link created')
      bond.show()
    i_seqs = [apply.atom1.i_seq, apply.atom2.i_seq]
    if self.geometry_proxy_registries.bond_simple.is_proxy_set(
      i_seqs=i_seqs,
      ):
      return False
    link_resolution = add_bond_proxies(
      counters=counters(label="apply_cif_link_bond"),
      m_i=m_i,
      m_j=m_j,
      bond_list=[bond],
      bond_simple_proxy_registry=self.geometry_proxy_registries
      .bond_simple,
      sites_cart=self.sites_cart,
      distance_cutoff=self.params.link_distance_cutoff,
      )
    link.bond_list.append(bond)
    # angles
    bonds1=[]
    bonds2=[]
    for bond in self.geometry_proxy_registries.bond_simple.proxies:
      if(apply.atom1.i_seq in bond.i_seqs and
         apply.atom2.i_seq in bond.i_seqs): continue
      if apply.atom1.i_seq in bond.i_seqs:
        bonds1.append(bond)
      elif apply.atom2.i_seq in bond.i_seqs:
        bonds2.append(bond)
    # remove excess hydrogens
    def _get_other(b,a):
      if a.i_seq==b.i_seqs[0]:
        return b.i_seqs[1]
      elif a.i_seq==b.i_seqs[1]:
        return b.i_seqs[0]
      return -1
    #
    def get_other(b,a):
      i_seq = _get_other(b,a)
      for atom in a.parent().parent().atoms():
        if atom.i_seq==i_seq: return atom
      return None
    #
    def remove_atom(m_x, atom_name):
      from mmtbx.monomer_library.cif_types import mod_mod_id, chem_mod, chem_mod_atom
      cm = chem_mod()
      cm.id = "DEL_%s_%s" % (m_x.residue_name, atom_name)
      cm.show()
      mmi = mod_mod_id("custom", cm)
      mod_atom = chem_mod_atom()
      mod_atom.function = "delete"
      mod_atom.atom_id = "%s" % atom_name
      #mod_atom.new_atom_id =
      #mod_atom.new_type_symbol =
      #mod_atom.new_type_energy =
      #mod_atom.new_partial_charge = 0.0
      mod_atom.show()
      mmi.atom_list.append(mod_atom)
      mmi.show()
      mmi = m_x.apply_mod(mmi)
      return mmi

    if order==1:
      if len(bonds1)>3:
        pass
      if len(bonds2)>3:
        pass
      if False:
        atoms = self.pdb_hierarchy.atoms()
        atom_names = []
        for bond in bonds2:
          print(bond.i_seqs, atoms[bond.i_seqs[0]].quote(), atoms[bond.i_seqs[1]].quote(),bond.distance_ideal)
          other = get_other(bond, apply.atom2)
          if other is None: continue
          if other.element.strip() not in ["H", "D"]: continue
          atom_names.append(other.name.strip())
        atom_names.sort()
        m_j = remove_atom(m_j, atom_names[0])
    #
    if not angles: return True
    angle_list = []
    for bond in bonds1:
      other = get_other(bond, apply.atom1)
      if other is None: continue
      angle = chem_link_angle()
      angle.atom_1_comp_id = 1
      angle.atom_id_1 = other.name.strip()
      angle.atom_2_comp_id = 1
      angle.atom_id_2 = apply.atom1.name.strip()
      angle.atom_3_comp_id = 2
      angle.atom_id_3 = apply.atom2.name.strip()
      if order==1:
        angle.value_angle = 109.6
      angle.value_angle_esd = 2.
      angle_list.append(angle)
    for bond in bonds2:
      other = get_other(bond, apply.atom2)
      if other is None: continue
      angle = chem_link_angle()
      angle.atom_1_comp_id = 2
      angle.atom_id_1 = other.name.strip()
      angle.atom_2_comp_id = 2
      angle.atom_id_2 = apply.atom2.name.strip()
      angle.atom_3_comp_id = 1
      angle.atom_id_3 = apply.atom1.name.strip()
      if order==1:
        angle.value_angle = 109.6
      angle.value_angle_esd = 2.
      angle_list.append(angle)
    #
    link_resolution = add_angle_proxies(
        counters=counters(label="apply_cif_link_angle"),
        m_i=m_i,
        m_j=m_j,
        angle_list=angle_list,
        angle_proxy_registry=self.geometry_proxy_registries.angle,
        special_position_dict=self.special_position_dict,
        )
    link.angle_list = angle_list
    return link

  def process_custom_links(self,
                           mon_lib_srv,
                           pdbres_pair,
                           data_link,
                           atoms,
                           indent=10,
                           verbose=False,
                           ):
    from mmtbx.monomer_library import linking_utils
    outl = ""
    classes1 = linking_utils.get_classes(atoms[0])
    classes2 = linking_utils.get_classes(atoms[1])
    # peptide links are auto-created
    possible_peptide_link = False
    if classes1.common_amino_acid or classes2.common_amino_acid:
      if(atoms[0].name.strip() in ["C"] and
         atoms[1].name.strip() in ["N"]
         ):
        possible_peptide_link=True
    # so are nucleotide links
    possible_rna_dna_link = False
    if ((classes1.common_rna_dna or classes1.ccp4_mon_lib_rna_dna)
          or (classes2.common_rna_dna or classes2.ccp4_mon_lib_rna_dna)):
      if(atoms[0].name.strip() in ["O3'", "O3*"] and
         atoms[1].name.strip() in ["P"]
         ):
        possible_rna_dna_link = True
    # add them
    ga = group_args(
      pdbres_pair=pdbres_pair,
      data_link=data_link,
      was_used=False,
      automatic=True,
      possible_peptide_link=possible_peptide_link,
      possible_rna_dna_link=possible_rna_dna_link,
      atom1=atoms[0],
      atom2=atoms[1],
      )
    # check if the link as already been added to apply list
    # could move this to class self.apply_cif_links
    count = 0
    matches = 2
    for i, apply in enumerate(self.apply_cif_links):
      count = 0
      if apply.pdbres_pair==pdbres_pair: count+=1
      if apply.data_link==data_link: count+=1
      if count==matches: break
    if count==matches:
      ga.was_used=True
    else:
      for pdbres in pdbres_pair:
        self.empty_apply_cif_links_mm_pdbres_dict[pdbres] = {}
    self.apply_cif_links.append(ga)
    # mods
    link = mon_lib_srv.link_link_id_dict.get(ga.data_link)
    if link is None: return outl
    #################################################
    # saccahrides that have non-standard atom names #
    #  in the names in the mods                     #
    #################################################
    ignore_mods = [False, False]
    if ga.atom1.parent().resname in ["BGC"]:
      ignore_mods[0]=True
    if ga.atom1.parent().resname in ["BGC"]:
      ignore_mods[1]=True
    #####################
    mod_ids = []
    for mod_i, mod_attr in enumerate(["mod_id_1", "mod_id_2"]):
      if ignore_mods[mod_i]: continue
      mod_id = getattr(link.chem_link, mod_attr)
      if (mod_id == ""): mod_id = None
      mod_ids.append(mod_id)
      if (mod_id is not None):
        outl += "%s%s: %s\n" % (" "*indent, mod_attr, mod_id)
        mod = mon_lib_srv.mod_mod_id_dict.get(mod_id)
        if (mod is None):
          print(outl)
          raise Sorry(
            "Missing CIF modification: data_mod_%s\n" % mod_id
            + "  Please check for spelling errors or specify the file name\n"
            + "  with the modification as an additional argument.")
    for pdbres,mod_id in zip(pdbres_pair, mod_ids):
      if (mod_id is not None):
        self.apply_cif_modifications.setdefault(pdbres, []).append(mod_id)
    return outl

  def create_disulfides(self, disulfide_distance_cutoff, log=None):
    if (self.model_indices is not None):
      model_indices = self.model_indices.select(self.cystein_sulphur_i_seqs)
    conformer_indices = self.conformer_indices.select(
      self.cystein_sulphur_i_seqs)
    sym_excl_indices = self.sym_excl_indices.select(
      self.cystein_sulphur_i_seqs)
    donor_acceptor_excl_groups = self.donor_acceptor_excl_groups.select(
      self.cystein_sulphur_i_seqs)
    asu_mappings = self.special_position_settings.asu_mappings(
      buffer_thickness=disulfide_distance_cutoff)
    sulphur_sites_cart = self.sites_cart.select(self.cystein_sulphur_i_seqs)
    asu_mappings.process_sites_cart(
      original_sites=sulphur_sites_cart,
      site_symmetry_table=self.site_symmetry_table().select(
        self.cystein_sulphur_i_seqs))
    pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    # keep this around ??????????
    nonbonded_proxies = geometry_restraints.nonbonded_sorted_asu_proxies(
      model_indices=model_indices,
      conformer_indices=conformer_indices,
      sym_excl_indices=sym_excl_indices,
      donor_acceptor_excl_groups=donor_acceptor_excl_groups,
      nonbonded_params=geometry_restraints.nonbonded_params(
        default_distance=1),
      nonbonded_types=flex.std_string(conformer_indices.size()),
      nonbonded_charges=flex.int(conformer_indices.size(), 0),
      nonbonded_distance_cutoff_plus_buffer=disulfide_distance_cutoff,
      min_cubicle_edge=5,
      shell_asu_tables=[pair_asu_table])
    labels = [self.pdb_atoms[i_seq].id_str()
      for i_seq in self.cystein_sulphur_i_seqs]
    for proxy in nonbonded_proxies.simple:
      pair_asu_table.add_pair(proxy.i_seqs)
    for proxy in nonbonded_proxies.asu:
      pair_asu_table.add_pair(proxy)
    pair_sym_table = pair_asu_table.extract_pair_sym_table()
    if (log is not None):
      n_simple, n_symmetry = 0, 0
      for sym_pair in pair_sym_table.iterator():
        if (sym_pair.rt_mx_ji.is_unit_mx()): n_simple += 1
        else:                                n_symmetry += 1
      print("  Number of disulfides: simple=%d, symmetry=%d" % (
        n_simple, n_symmetry), file=log)
      if (n_symmetry == 0):
        blanks = ""
      else:
        blanks = "  "
    max_distance_model = 0
    frac = asu_mappings.unit_cell().fractionalize
    orth = asu_mappings.unit_cell().orthogonalize
    for sym_pair in pair_sym_table.iterator():
      # Symmetry-aware distance calculation between 2 atoms
      distance_model = geometry_restraints.bond(
        [sulphur_sites_cart[sym_pair.i_seq],
         orth(sym_pair.rt_mx_ji*frac(sulphur_sites_cart[sym_pair.j_seq]))],
        distance_ideal=0,
        weight=1).distance_model
      max_distance_model = max(max_distance_model, distance_model)
      if (log is not None):
        if (sym_pair.rt_mx_ji.is_unit_mx()):
          disulfide_type = "Simple disulfide:%s" % blanks
        else:
          disulfide_type = "Symmetry disulfide:"
        print("    %s %s - %s distance=%.2f" % tuple(
          [disulfide_type]
          + [labels[i_seq] for i_seq in sym_pair.i_seqs()]
          + [distance_model]), end='', file=log)
        if (not sym_pair.rt_mx_ji.is_unit_mx()):
          print(" " + str(sym_pair.rt_mx_ji), end='', file=log)
        print(file=log)
    return pair_sym_table, max_distance_model

  def atom_selection(self, parameter_name, string, cache=None):
    try:
      return self.selection(string=string, cache=cache)
    except KeyboardInterrupt: raise
    except Exception as e: # keep e alive to avoid traceback
      fe = format_exception()
      raise Sorry('Invalid atom selection:\n  %s="%s"\n  (%s)' % (
        parameter_name, string, fe))

  def phil_atom_selection(self,
        cache,
        scope_extract,
        attr,
        string=None,
        allow_none=False,
        allow_auto=False,
        raise_if_empty_selection=True):
    def parameter_name():
      return scope_extract.__phil_path__(object_name=attr)
    if (string is None):
      string = getattr(scope_extract, attr)
    if (string is None):
      if (allow_none): return None
      raise Sorry('Atom selection cannot be None:\n  %s=None' % (
        parameter_name()))
    elif (string is Auto):
      if (allow_auto): return Auto
      raise Sorry('Atom selection cannot be Auto:\n  %s=Auto' % (
        parameter_name()))
    try: result = self.selection(string=string, cache=cache)
    except KeyboardInterrupt: raise
    except Exception as e: # keep e alive to avoid traceback
      fe = format_exception()
      raise Sorry('Invalid atom selection:\n  %s=%s\n  (%s)' % (
        parameter_name(), show_string(string), fe))
    if (raise_if_empty_selection and result.count(True) == 0):
      raise Sorry('Empty atom selection:\n  %s=%s' % (
        parameter_name(), show_string(string)))
    return result

  def phil_atom_selection_multiple(self,
        cache,
        scope_extract,
        attr,
        allow_none=False,
        allow_auto=False,
        raise_if_empty_selection=True):
    result = []
    def parameter_name():
      return scope_extract.__phil_path__(object_name=attr)
    string_list = getattr(scope_extract, attr)
    for string in string_list:
      if (string is None):
        if (allow_none): return None
        raise Sorry('Atom selection cannot be None:\n  %s=None' % (
          parameter_name()))
      elif (string is Auto):
        if (allow_auto): return Auto
        raise Sorry('Atom selection cannot be Auto:\n  %s=Auto' % (
          parameter_name()))
      try:
        result.append(self.selection(string=string, cache=cache).iselection())
      except KeyboardInterrupt: raise
      except Exception as e: # keep e alive to avoid traceback
        fe = format_exception()
        raise Sorry('Invalid atom selection:\n  %s=%s\n  (%s)' % (
          parameter_name(), show_string(string), fe))
      if (raise_if_empty_selection and result.count(True) == 0):
        raise Sorry('Empty atom selection:\n  %s=%s' % (
          parameter_name(), show_string(string)))
    return result

  def phil_atom_selections_as_i_seqs(self, cache, scope_extract, sel_attrs, n_atoms_needed="one"):
    """ This function should be used when selection parameter is .multiple=False
    On top of that it can demand exactly one atom to be selected or more than 2
    atoms shoud be selected. Other options are not required at the moment.
    n_atoms_needed:
      one - exactly one
      many - 3 and more for parallelity
    """
    assert n_atoms_needed in ["one", "many"]
    result = []
    for attr in sel_attrs:
      iselection = self.phil_atom_selection(
        cache=cache,
        scope_extract=scope_extract,
        attr=attr,
        raise_if_empty_selection=False).iselection()
      if (iselection.size() != 1) and n_atoms_needed=="one":
        atom_sel = getattr(scope_extract, attr)
        if (iselection.size() == 0):
          raise Sorry("No atom selected: %s" % show_string(atom_sel))
        else:
          raise Sorry(
            "More than one atom selected: %s\n"
            "  Number of selected atoms: %d" % (
              show_string(atom_sel), iselection.size()))
      if (iselection.size() < 3) and n_atoms_needed=="many":
        atom_sel = getattr(scope_extract, attr)
        raise Sorry(
          "Slelected less than 3 atoms: %s\n"
          "  Number of selected atoms: %d" % (
            show_string(atom_sel), iselection.size()))
      if n_atoms_needed == "one":
        result.append(iselection[0])
      elif n_atoms_needed == "many":
        result.append(iselection)
    return result

  def phil_atom_selections_as_i_seqs_multiple(self,
                                              cache,
                                              scope_extract,
                                              sel_attrs):
    """ This function should be used when selection parameter is .multiple=True
    """
    result = []
    for attr in sel_attrs:
        iselection = self.phil_atom_selection_multiple(
          cache=cache,
          scope_extract=scope_extract,
          attr=attr,
          raise_if_empty_selection=False)
        atom_sel = getattr(scope_extract, attr)
        for i, i_sel in enumerate(iselection):
          if (i_sel.size() == 0):
            raise Sorry("No atom selected: %s" % show_string(atom_sel[i]))
          for atom in i_sel:
            result.append(atom)
    return result

  def process_geometry_restraints_remove(self,
        params,
        geometry_restraints_manager):
    path = params.__phil_path__() + "."
    grm = geometry_restraints_manager
    cache = self.pdb_hierarchy.atom_selection_cache()
    for atom_selection in params.angles:
      grm.remove_angles_in_place(selection=self.atom_selection(
        parameter_name=path+"angles", string=atom_selection, cache=cache))
    for atom_selection in params.dihedrals:
      grm.remove_dihedrals_in_place(selection=self.atom_selection(
        parameter_name=path+"dihedrals", string=atom_selection, cache=cache))
    for atom_selection in params.chiralities:
      grm.remove_chiralities_in_place(selection=self.atom_selection(
        parameter_name=path+"chiralities", string=atom_selection, cache=cache))
    for atom_selection in params.planarities:
      grm.remove_planarities_in_place(selection=self.atom_selection(
        parameter_name=path+"planarities", string=atom_selection, cache=cache))
    for atom_selection in params.parallelities:
      grm.remove_parallelities_in_place(selection=self.atom_selection(
        parameter_name=path+"parallelities", string=atom_selection, cache=cache))

  def process_geometry_restraints_edits_bond(self, sel_cache, params, log):
    bond_sym_proxies = []
    bond_distance_model_max = 0
    if (len(params.bond) == 0):
      return group_args(
        bond_sym_proxies=bond_sym_proxies,
        bond_distance_model_max=bond_distance_model_max)
    print("  Custom bonds:", file=log)
    atoms = self.pdb_atoms
    unit_cell = self.special_position_settings.unit_cell()
    space_group = self.special_position_settings.space_group()
    uc_shortest_vector = unit_cell.shortest_vector_sq()**0.5
    max_bond_length = uc_shortest_vector
    ebdl = params.excessive_bond_distance_limit
    if (ebdl not in [None, Auto] and ebdl > 0 and ebdl < max_bond_length):
      max_bond_length = ebdl
    n_excessive = 0
    sel_attrs = ["atom_selection_"+n for n in ["1", "2"]]
    self.pdb_link_records.setdefault("LINK", [])
    for bond in params.bond:
      def show_atom_selections():
        print(get_atom_selections_text(), end='', file=log)
      def get_atom_selections_text():
        txt = ""
        for attr in sel_attrs:
          txt += "      %s = %s\n" % (attr, show_string(getattr(bond, attr, None)))
        return txt
      slack = bond.slack
      if (slack is None or slack < 0):
        slack = 0
      if (bond.distance_ideal is None):
        print("    Warning: Ignoring bond with distance_ideal = None:", file=log)
        show_atom_selections()
      elif (bond.distance_ideal < 0):
        print("    Warning: Ignoring bond with distance_ideal < 0:", file=log)
        show_atom_selections()
        print("      distance_ideal = %.6g" % bond.distance_ideal, file=log)
      elif (bond.sigma is None):
        print("    Warning: Ignoring bond with sigma = None:", file=log)
        show_atom_selections()
        print("      distance_ideal = %.6g" % bond.distance_ideal, file=log)
      elif (bond.sigma <= 0):
        print("    Warning: Ignoring bond with sigma <= 0:", file=log)
        show_atom_selections()
        print("      distance_ideal = %.6g" % bond.distance_ideal, file=log)
        print("      sigma = %.6g" % bond.sigma, file=log)
        print("      slack = %.6g" % slack, file=log)
      elif (bond.action == "delete"):
        raise Sorry("%s = %s not implemented." %
          bond.__phil_path_and_value__(object_name="action"))
      else:
        i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache, scope_extract=bond, sel_attrs=sel_attrs)
        bond_exist = self.geometry_proxy_registries.bond_simple.is_proxy_set(i_seqs)
        if bond_exist and bond.action == 'add':
          txt = get_atom_selections_text()
          raise Sorry("Bond below exists, use action=change instead.\n" + txt)
        if not bond_exist and bond.action == 'change':
          txt = get_atom_selections_text()
          raise Sorry("Bond below does not exists, use action=add instead.\n" + txt)
        if (bond.symmetry_operation is None):
          s = "x,y,z"
        else:
          s = bond.symmetry_operation
        rt_mx_ji = sgtbx.rt_mx(symbol=s, t_den=space_group.t_den())
        self.pdb_link_records["LINK"].append([self.pdb_atoms[i_seqs[0]],
                                              self.pdb_atoms[i_seqs[1]],
                                              rt_mx_ji,
                                            ])
        p = geometry_restraints.bond_sym_proxy(
          i_seqs=i_seqs,
          distance_ideal=bond.distance_ideal,
          weight=geometry_restraints.sigma_as_weight(sigma=bond.sigma),
          slack=slack,
          limit=bond.limit,
          top_out=bond.top_out,
          origin_id=origin_ids.get_origin_id('edits'),
          rt_mx_ji=rt_mx_ji)
        bond_sym_proxies.append(p)
        b = geometry_restraints.bond(
          unit_cell=unit_cell,
          sites_cart=self.sites_cart,
          proxy=p)
        print("    bond:", file=log)
        for i in [0,1]:
          print("      atom %d:" % (i+1), atoms[p.i_seqs[i]].quote(), file=log)
        print("      symmetry operation:", str(p.rt_mx_ji), file=log)
        if (not space_group.contains(smx=p.rt_mx_ji)):
          raise Sorry(
            'The bond symmetry operation "%s" is not compatible'
            ' with space group %s.' % (
              str(p.rt_mx_ji),
              self.special_position_settings.space_group_info()
                .symbol_and_number()))
        print("      distance_model: %7.3f" % b.distance_model, file=log)
        print("      distance_ideal: %7.3f" % b.distance_ideal, file=log)
        print("      ideal - model:  %7.3f" % b.delta, file=log)
        print("      slack:          %7.3f" % b.slack, file=log)
        print("      delta_slack:    %7.3f" % b.delta_slack, file=log)
        print("      sigma:          %8.4f" % \
          geometry_restraints.weight_as_sigma(weight=b.weight), file=log)
        if (bond_distance_model_max < b.distance_model):
          bond_distance_model_max = b.distance_model
        if (b.distance_model > max_bond_length):
          print("      *** WARNING: EXCESSIVE BOND LENGTH. ***", file=log)
          n_excessive += 1
    if (n_excessive != 0):
      if (max_bond_length == uc_shortest_vector):
        print("  Excessive bond length limit at hard upper bound:" \
          " length of shortest vector between unit cell lattice points: %.6g" \
            % uc_shortest_vector, file=log)
      else:
        print("  %s = %.6g" % \
          params.__phil_path_and_value__("excessive_bond_distance_limit"), file=log)
        print("    Please assign a larger value to this parameter if necessary.", file=log)
      raise Sorry(
        "Custom bonds with excessive length: %d\n"
        "  Please check the log file for details." % n_excessive)
    print("    Total number of added/changed bonds:", len(bond_sym_proxies), file=log)
    return group_args(
      bond_sym_proxies=bond_sym_proxies,
      bond_distance_model_max=bond_distance_model_max)

  def process_geometry_restraints_edits_angle(self,
                                              sel_cache,
                                              params,
                                              log,
                                              second_pass=False):
    result = []
    if (len(params.angle) == 0): return result
    if (self.special_position_indices is None):
      special_position_indices = []
    else:
      special_position_indices = self.special_position_indices
    print("  Custom angles:", file=log)
    atoms = self.pdb_atoms
    n_changed_angles = 0
    sel_attrs = ["atom_selection_"+n for n in ["1", "2", "3"]]
    for angle in params.angle:
      def show_atom_selections():
        print(get_atom_selections_text(), end='', file=log)
      def get_atom_selections_text():
        txt = ""
        for attr in sel_attrs:
          txt += "      %s = %s\n" % (attr, show_string(getattr(angle, attr, None)))
        return txt
      if (angle.angle_ideal is None):
        print("    Warning: Ignoring angle with angle_ideal = None:", file=log)
        show_atom_selections()
      elif (angle.sigma is None):
        print("    Warning: Ignoring angle with sigma = None:", file=log)
        show_atom_selections()
        print("      angle_ideal = %.6g" % angle.angle_ideal, file=log)
      elif (angle.sigma is None or angle.sigma <= 0):
        print("    Warning: Ignoring angle with sigma <= 0:", file=log)
        show_atom_selections()
        print("      angle_ideal = %.6g" % angle.angle_ideal, file=log)
        print("      sigma = %.6g" % angle.sigma, file=log)
      elif (angle.action == "change"):
        if not second_pass: pass
        i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache, scope_extract=angle, sel_attrs=sel_attrs)
        i_proxy = self.geometry_proxy_registries.angle.lookup_i_proxy(i_seqs)
        if i_proxy is None:
          txt = get_atom_selections_text()
          raise Sorry("Angle below is not restrained, nothing to change.\n" + txt)
        a_proxy = self.geometry_proxy_registries.angle.proxies[i_proxy]
        a_proxy.angle_ideal=angle.angle_ideal
        a_proxy.weight = geometry_restraints.sigma_as_weight(sigma=angle.sigma)
        a_proxy.origin_id=origin_ids.get_origin_id('edits')
        n_changed_angles += 1
      elif (angle.action != "add"):
        raise Sorry("%s = %s not implemented." %
          angle.__phil_path_and_value__("action"))
      else:
        i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache, scope_extract=angle, sel_attrs=sel_attrs)
        p = geometry_restraints.angle_proxy(
          i_seqs=i_seqs,
          angle_ideal=angle.angle_ideal,
          origin_id=origin_ids.get_origin_id('edits'),
          weight=geometry_restraints.sigma_as_weight(sigma=angle.sigma))
        a = geometry_restraints.angle(
          sites_cart=self.sites_cart,
          proxy=p)
        print("    angle:", file=log)
        n_special = 0
        for i,i_seq in enumerate(p.i_seqs):
          print("      atom %d:" % (i+1), atoms[i_seq].quote(), end=' ', file=log)
          if (i_seq in special_position_indices):
            n_special += 1
            print("# SPECIAL POSITION", end=' ', file=log)
          print(file=log)
        print("      angle_model: %7.2f" % a.angle_model, file=log)
        print("      angle_ideal: %7.2f" % a.angle_ideal, file=log)
        print("      ideal - model:  %7.2f" % a.delta, file=log)
        print("      sigma: %.6g" % \
          geometry_restraints.weight_as_sigma(weight=a.weight), file=log)
        if (n_special != 0):
          raise Sorry(
            "Custom angle involves %d special position%s:\n"
            "  Please inspect the output for details."
              % plural_s(n_special))
        result.append(p)
    print("    Total number of new custom angles:", len(result), file=log)
    print("    Total number of changed angles:", n_changed_angles, file=log)
    return result

  def process_geometry_restraints_edits_dihedral(self, sel_cache, params, log):
    result = []
    if (len(params.dihedral) == 0): return result
    if (self.special_position_indices is None):
      special_position_indices = []
    else:
      special_position_indices = self.special_position_indices
    print("  Custom dihedrals:", file=log)
    atoms = self.pdb_atoms
    sel_attrs = ["atom_selection_"+n for n in ["1", "2", "3", "4"]]
    for dihedral in params.dihedral:
      def show_atom_selections():
        print(get_atom_selections_text(), end='', file=log)
      def get_atom_selections_text():
        txt = ""
        for attr in sel_attrs:
          txt += "      %s = %s\n" % (attr, show_string(getattr(dihedral, attr, None)))
        return txt
      if (dihedral.angle_ideal is None):
        print("    Warning: Ignoring dihedral with angle_ideal = None:", file=log)
        show_atom_selections()
      elif (dihedral.sigma is None):
        print("    Warning: Ignoring dihedral with sigma = None:", file=log)
        show_atom_selections()
        print("      angle_ideal = %.6g" % dihedral.angle_ideal, file=log)
      elif (dihedral.sigma is None or dihedral.sigma <= 0):
        print("    Warning: Ignoring dihedral with sigma <= 0:", file=log)
        show_atom_selections()
        print("      angle_ideal = %.6g" % dihedral.angle_ideal, file=log)
        print("      sigma = %.6g" % dihedral.sigma, file=log)
      elif dihedral.action == "change":
        i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache, scope_extract=dihedral, sel_attrs=sel_attrs)
        i_proxy = self.geometry_proxy_registries.dihedral.lookup_i_proxy(i_seqs)[0]
        if i_proxy is None:
          txt = get_atom_selections_text()
          raise Sorry("Angle below is not restrained, nothing to change.\n" + txt)
        a_proxy = self.geometry_proxy_registries.dihedral.proxies[i_proxy]
        a_proxy.angle_ideal=dihedral.angle_ideal
        a_proxy.weight = geometry_restraints.sigma_as_weight(sigma=dihedral.sigma)
        a_proxy.periodicity = dihedral.periodicity
        a_proxy.alt_angle_ideals = dihedral.alt_angle_ideals
        a_proxy.origin_id=origin_ids.get_origin_id('edits')
      elif (dihedral.action != "add"):
        raise Sorry("%s = %s not implemented." %
          dihedral.__phil_path_and_value__("action"))
      else:
        i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache, scope_extract=dihedral, sel_attrs=sel_attrs)
        p = geometry_restraints.dihedral_proxy(
          i_seqs=i_seqs,
          angle_ideal=dihedral.angle_ideal,
          alt_angle_ideals=dihedral.alt_angle_ideals,
          weight=geometry_restraints.sigma_as_weight(sigma=dihedral.sigma),
          origin_id=origin_ids.get_origin_id('edits'),
          periodicity=dihedral.periodicity,
        )
        a = geometry_restraints.dihedral(
          sites_cart=self.sites_cart,
          proxy=p)
        print("    dihedral:", file=log)
        n_special = 0
        for i,i_seq in enumerate(p.i_seqs):
          print("      atom %d:" % (i+1), atoms[i_seq].quote(), end=' ', file=log)
          if (i_seq in special_position_indices):
            n_special += 1
            print("# SPECIAL POSITION", end=' ', file=log)
          print(file=log)
        print("      angle_model: %7.2f" % a.angle_model, file=log)
        print("      angle_ideal: %7.2f" % a.angle_ideal, file=log)
        print("      ideal - model:  %7.2f" % a.delta, file=log)
        print("      sigma: %.6g" % \
          geometry_restraints.weight_as_sigma(weight=a.weight), file=log)
        if (n_special != 0):
          raise Sorry(
            "Custom dihedral involves %d special position%s:\n"
            "  Please inspect the output for details."
              % plural_s(n_special))
        result.append(p)
    print("    Total number of custom dihedrals:", len(result), file=log)
    return result

  def process_geometry_restraints_edits_planarity(self,
                                                  sel_cache,
                                                  params,
                                                  log):
    result = []
    if (len(params.planarity) == 0): return result
    if (self.special_position_indices is None):
      special_position_indices = []
    else:
      special_position_indices = self.special_position_indices
    sel_attrs = ["atom_selection"]
    print("  Custom planarities:", file=log)
    for planarity in params.planarity:
      def show_atom_selections():
        print("      %s = %s" % (
            "atom_selection", planarity.atom_selection), file=log)
      if (planarity.sigma is None) or (planarity.sigma <= 0):
        print("    Warning: Ignoring planarity with with sigma <= 0:", file=log)
        print(show_atom_selections(), file=log)
        continue
        # raise Sorry("Custom planarity sigma is undefined or zero/negative - "+
        #   "this must be a positive decimal number.")
      i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache,
          scope_extract=planarity,
          sel_attrs=sel_attrs,
          n_atoms_needed="many")[0]
      weights = []
      for i_seq in i_seqs:
        weights.append(geometry_restraints.sigma_as_weight(
          sigma=planarity.sigma))
      proxy = geometry_restraints.planarity_proxy(
        i_seqs=i_seqs,
        origin_id=origin_ids.get_origin_id('edits'),
        weights=flex.double(weights))
      plane = geometry_restraints.planarity(
        sites_cart=self.sites_cart,
        proxy=proxy)
      result.append(proxy)
    print("    Total number of custom planarities:", len(result), file=log)
    return result

  def process_geometry_restraints_edits_parallelity(self,
                                                  sel_cache,
                                                  params,
                                                  log):
    result = []
    if len(params.parallelity) == 0:
      return result
    print("  Custom parallelities:", file=log)
    for parallelity in params.parallelity:
      if (parallelity.atom_selection_1 is None or
          parallelity.atom_selection_2 is None):
        print("Warning: Ignoring parallelity with empty atom selection.", file=log)
        continue
      if (parallelity.sigma is None) or (parallelity.sigma <= 0):
        raise Sorry("Custom parallelity sigma is undefined or zero/negative - "+
          "this must be a positive decimal number.")
      elif (parallelity.action != "add"):
        raise Sorry("%s = %s not implemented." %
          parallelity.__phil_path_and_value__("action"))
      i_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache,
          scope_extract=parallelity,
          sel_attrs=["atom_selection_1"],
          n_atoms_needed="many")[0]
      j_seqs = self.phil_atom_selections_as_i_seqs(
          cache=sel_cache,
          scope_extract=parallelity,
          sel_attrs=["atom_selection_2"],
          n_atoms_needed="many")[0]
      weight = 1./parallelity.sigma**2
      target_angle_deg = parallelity.target_angle_deg
      # print i_seqs, j_seqs, weight
      proxy = geometry_restraints.parallelity_proxy(
        i_seqs=flex.size_t(i_seqs),
        j_seqs=flex.size_t(j_seqs),
        weight=weight,
        origin_id=origin_ids.get_origin_id('edits'),
        target_angle_deg=target_angle_deg)
      result.append(proxy)
    print("    Total number of custom parallelities:", len(result), file=log)
    return result


  def process_geometry_restraints_scale(self, params, log):
    """
    Scale the weights for selected basic geometry restraints for given
    atom selections.  This allows geometry to be manually tightened or
    loosened for problematic atoms, without making global changes to the
    geometry target weight or editing CIF files.  Since this is still somewhat
    dangerous, atom selections are checked for overlap, and each restraint
    proxy will be modified no more than once.

    Note that planarity proxies are currently left alone, since they have
    an array of weights instead of a scalar value.
    """
    return # obsoleted
    sel_cache = self.pdb_hierarchy.atom_selection_cache()
    other_selections = []
    other_selection_strs = []
    if (len(params.scale_restraints) > 0):
      print("Scaling restraint weights for %d selections" % \
        len(params.scale_restraints), file=log)
    for scale_params in params.scale_restraints :
      if (scale_params.scale < 0):
        raise Sorry("scale_restraints.scale must be at least zero.")
      selection = self.phil_atom_selection(
        cache=sel_cache,
        scope_extract=scale_params,
        attr="atom_selection",
        raise_if_empty_selection=True)
      for other, other_str in zip(other_selections, other_selection_strs):
        if (not (selection & other).all_eq(False)):
          raise Sorry(("Error scaling selected restraint weights: the atom "+
            "selection \"%s\" overlaps with at least one other previously "+
            "defined atom selection (\"%s\").") % (scale_params.atom_selection,
              other_str))
      other_selections.append(selection)
      other_selection_strs.append(scale_params.atom_selection)
      proxy_lists = [ self.geometry_proxy_registries.angle.proxies,
                      self.geometry_proxy_registries.dihedral.proxies,
                      self.geometry_proxy_registries.chirality.proxies,
                      self.geometry_proxy_registries.bond_simple.proxies ]
      proxy_types = ["angle", "dihedral", "chirality", "bond"]
      modified_proxies = []
      for k, proxies in enumerate(proxy_lists):
        if (not proxy_types[k] in scale_params.apply_to):
          continue
        for j, proxy in enumerate(proxies):
          for i_seq in proxy.i_seqs :
            if (selection[i_seq]):
              if ((k,j) in modified_proxies):
                print("  skipping %s restraint proxy #%d - already modified" % (
                    proxy_types[k], j), file=log)
                print("  atoms involved:", file=log)
                for i_seq_2 in proxy.i_seqs :
                  print("    %s" % self.pdb_atoms[i_seq_2].id_str(), file=log)
                continue
              proxy.weight *= scale_params.scale
              modified_proxies.append((k,j))
              break
      # TODO: planarity?

  def process_geometry_restraints_edits(self, params, log, second_pass=False):
    sel_cache = self.pdb_hierarchy.atom_selection_cache()
    result = self.process_geometry_restraints_edits_bond(
        sel_cache=sel_cache, params=params, log=log)
    result.angle_proxies=self.process_geometry_restraints_edits_angle(
        sel_cache=sel_cache, params=params, log=log, second_pass=second_pass)
    result.dihedral_proxies=self.process_geometry_restraints_edits_dihedral(
        sel_cache=sel_cache, params=params, log=log)
    result.planarity_proxies=self.process_geometry_restraints_edits_planarity(
        sel_cache=sel_cache, params=params, log=log)
    result.parallelity_proxies=self.process_geometry_restraints_edits_parallelity(
        sel_cache=sel_cache, params=params, log=log)
    return result

  def process_hydrogen_bonds(self, bonds_table, log, verbose=False):
    atoms = self.pdb_atoms
    def show_atoms(i_seqs, log):
      for i_seq in i_seqs :
        print("     %s" % atoms[i_seq].fetch_labels().quote(), file=log)
    unit_cell = self.special_position_settings.unit_cell()
    space_group = self.special_position_settings.space_group()
    uc_shortest_vector = unit_cell.shortest_vector_sq()**0.5
    max_bond_length = uc_shortest_vector
    n_excessive = 0
    bond_distance_model_max = 0
    bond_sym_proxies = []
    for bond in bonds_table.get_bond_restraint_data():
      i_seqs = [bond.donor_i_seq, bond.acceptor_i_seq]
      slack = bond.slack
      if (slack is None or slack < 0):
        slack = 0
      if (bond.distance_ideal is None):
        print("    Warning: Ignoring bond with distance_ideal = None:", file=log)
        show_atoms(i_seqs, log)
      elif (bond.distance_ideal <= 0):
        print("    Warning: Ignoring bond with distance_ideal <= 0:", file=log)
        show_atoms(i_seqs, log)
        print("      distance_ideal = %.6g" % bond.distance_ideal, file=log)
      elif (bond.sigma is None):
        print("    Warning: Ignoring bond with sigma = None:", file=log)
        show_atoms(i_seqs, log)
        print("      distance_ideal = %.6g" % bond.distance_ideal, file=log)
      elif (bond.sigma <= 0):
        print("    Warning: Ignoring bond with sigma <= 0:", file=log)
        show_atoms(i_seqs, log)
        print("      distance_ideal = %.6g" % bond.distance_ideal, file=log)
        print("      sigma = %.6g" % bond.sigma, file=log)
        print("      slack = %.6g" % slack, file=log)
      else:
        rt_mx_ji = sgtbx.rt_mx(symbol="x,y,z", t_den=space_group.t_den())
        p = geometry_restraints.bond_sym_proxy(
          i_seqs=i_seqs,
          distance_ideal=bond.distance_ideal,
          weight=geometry_restraints.sigma_as_weight(sigma=bond.sigma),
          slack=slack,
          rt_mx_ji=rt_mx_ji)
        bond_sym_proxies.append(p)
        b = geometry_restraints.bond(
          unit_cell=unit_cell,
          sites_cart=self.sites_cart,
          proxy=p)
        if (b.distance_model > max_bond_length):
          print("      *** WARNING: EXCESSIVE BOND LENGTH. ***", file=log)
          n_excessive += 1
        if verbose :
          print("    hydrogen bond:", file=log)
          for i in [0,1]:
            print("      atom %d:" % (i+1), atoms[p.i_seqs[i]].quote(), file=log)
          print("      distance_model: %7.3f" % b.distance_model, file=log)
          print("      distance_ideal: %7.3f" % b.distance_ideal, file=log)
          print("      ideal - model:  %7.3f" % b.delta, file=log)
          print("      slack:          %7.3f" % b.slack, file=log)
          print("      delta_slack:    %7.3f" % b.delta_slack, file=log)
          print("      sigma:          %8.4f" % \
            geometry_restraints.weight_as_sigma(weight=b.weight), file=log)
        if (bond_distance_model_max < b.distance_model):
          bond_distance_model_max = b.distance_model
    if (n_excessive != 0):
      print("  Excessive bond length limit at hard upper bound:" \
        " length of shortest vector between unit cell lattice points: %.6g" \
          % uc_shortest_vector, file=log)
      raise Sorry(
        "Hydrogen bonds with excessive length: %d\n"
        "  Please check the log file for details." % n_excessive)
    print("  Total number of hydrogen bonds:", len(bond_sym_proxies), file=log)
    return group_args(
      bond_sym_proxies=bond_sym_proxies,
      bond_distance_model_max=bond_distance_model_max)

  def process_custom_nonbonded_exclusions(self, log, exclude_pair_indices,
      shell_asu_tables, verbose=True):
    space_group = self.special_position_settings.space_group()
    rt_mx_ji = sgtbx.rt_mx(symbol="x,y,z", t_den=space_group.t_den())
    have_header = False
    for (i_seq, j_seq) in exclude_pair_indices :
      if (verbose) and (not have_header):
        print(file=log)
        print("  Custom nonbonded exclusions (H-bonds, etc.):", file=log)
        have_header = True
      if (verbose):
        print("    %s  %s" % (self.pdb_atoms[i_seq].id_str(),
                                      self.pdb_atoms[j_seq].id_str()), file=log)
      try :
        shell_asu_tables[1].add_pair(i_seq, j_seq, rt_mx_ji)
      except RuntimeError as e :
        print("    WARNING: could not process nonbonded pair", file=log)
        print("    Original error:", file=log)
        print("      %s" % str(e), file=log)
      #shell_sym_tables[1][i_seq][j_seq].add(rt_mx_ji)

  def process_custom_nonbonded_symmetry_exclusions(self,
        log, curr_sym_excl_index):
    have_header = False
    sel_cache = None
    sel_full_occ = None
    for sel_string in self.params.custom_nonbonded_symmetry_exclusions:
      if (sel_string is None):
        continue
      if (sel_cache is None):
        sel_cache = self.pdb_hierarchy.atom_selection_cache()
        sel_full_occ = (self.pdb_atoms.extract_occ() == 1)
      sel = self.phil_atom_selection(
        cache=sel_cache,
        scope_extract=self.params,
        attr="custom_nonbonded_symmetry_exclusions",
        string=sel_string)
      if (sel is not None):
        isel = sel.iselection()
        if (not have_header):
          print(file=log)
          print("  Custom nonbonded symmetry exclusions:", file=log)
          have_header = True
        print("    Atom selection:", sel_string, file=log)
        print("      Number of atoms selected:", isel.size(), file=log)
        prev_sym_excl_indices = self.sym_excl_indices.select(isel)
        n_prev = isel.size() - prev_sym_excl_indices.count(0)
        if (n_prev != 0):
          if (n_prev == 1): s = ""
          else:             s = "s"
          print("      WARNING: %d atom%s in previous symmetry exclusion group%s" \
              % (n_prev, s, s), file=log)
          i_seq = isel.select(prev_sym_excl_indices != 0)[0]
          print("        Example:", self.pdb_atoms[i_seq].id_str(), file=log)
        sel_sel_full_occ = sel_full_occ.select(isel)
        n_full_occ = sel_sel_full_occ.count(True)
        if (n_full_occ != 0):
          print("      WARNING: %d atom%s with full occupancy" \
            % plural_s(n_full_occ), file=log)
          i_seq = isel.select(sel_sel_full_occ)[0]
          print("        Example:", self.pdb_atoms[i_seq].id_str(), file=log)
        sp = sorted(set(self.special_position_indices).intersection(set(isel)))
        if (len(sp) != 0):
          if (len(sp) == 1):
            print("      WARNING: one atom at a special position", file=log)
          else:
            print("      WARNING: %d atoms at special positions" \
              % len(sp), file=log)
          i_seq = sp[0]
          print("        Example:", self.pdb_atoms[i_seq].id_str(), file=log)
        curr_sym_excl_index += 1
        self.sym_excl_indices.set_selected(isel, curr_sym_excl_index)
    if (have_header):
      print(file=log)

  def construct_geometry_restraints_manager(self,
        ener_lib,
        disulfide_link,
        plain_pairs_radius=None,
        params_edits=None,
        params_remove=None,
        custom_nonbonded_exclusions=None,
        assume_hydrogens_all_missing=True,
        external_energy_function=None,
        log=None):
    assert self.special_position_settings is not None
    timer = user_plus_sys_time()
    #if (params_edits is not None):
    #  self.process_geometry_restraints_scale(params_edits, log)
    bond_params_table = geometry_restraints.extract_bond_params(
      n_seq=self.sites_cart.size(),
      bond_simple_proxies=self.geometry_proxy_registries.bond_simple.proxies)
    bond_distances_model = geometry_restraints.bond_distances_model(
      sites_cart=self.sites_cart,
      proxies=self.geometry_proxy_registries.bond_simple.proxies)
    if (bond_distances_model.size() > 0):
      excessive_bonds = (
        bond_distances_model > self.special_position_settings.unit_cell()
          .shortest_vector_sq()**.5).iselection()
      if(excessive_bonds.size() > 0):
        if(not self.params.proceed_with_excessive_length_bonds):
          atoms = self.pdb_atoms
          proxies = self.geometry_proxy_registries.bond_simple.proxies
          print("  Bonds with excessive lengths:", file=log)
          for i_proxy in excessive_bonds:
            proxy = proxies[i_proxy]
            bond = geometry_restraints.bond(
              sites_cart=self.sites_cart, proxy=proxy)
            print("    Distance model: %.6g (ideal: %.6g)" % (
              bond.distance_model, bond.distance_ideal), file=log)
            for i_seq in proxy.i_seqs:
              print("      %s" % atoms[i_seq].format_atom_record(), file=log)
          raise Sorry("Number of bonds with excessive lengths: %d" %
            excessive_bonds.size())
        else:
          self.params.max_reasonable_bond_distance = \
            flex.max(bond_distances_model)*2
    # disulfides
    disulfide_sym_table, max_disulfide_bond_distance = \
      self.create_disulfides(
        disulfide_distance_cutoff=self.params.disulfide_distance_cutoff,
        log=log)
    disulfide_cif_block = None
    disulfide_cif_loop = None
    max_bond_distance = max_disulfide_bond_distance
    if (bond_distances_model.size() > 0):
      max_bond_distance = max(max_bond_distance,
        flex.max(bond_distances_model))
    if (params_edits is None):
      processed_edits = None
    else:
      processed_edits = self.process_geometry_restraints_edits(
        params=params_edits, log=log)
      max_bond_distance = max(max_bond_distance,
        processed_edits.bond_distance_model_max)
    asu_mappings = self.special_position_settings.asu_mappings(
      buffer_thickness=max_bond_distance*3,
      )
    # factor 3 is to reach 1-4 interactions
    asu_mappings.process_sites_cart(
      original_sites=self.sites_cart,
      site_symmetry_table=self.site_symmetry_table())
    bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)

    geometry_restraints.add_pairs(
      bond_asu_table, self.geometry_proxy_registries.bond_simple.proxies)
    #
    disulfide_bond = disulfide_link.bond_list[0]
    disulfide_angle = disulfide_link.angle_list[0]
    disulfide_torsions = disulfide_link.tor_list

    assert disulfide_bond.value_dist is not None
    assert disulfide_bond.value_dist > 0
    assert disulfide_bond.value_dist_esd is not None
    assert disulfide_bond.value_dist_esd > 1e-5

    assert disulfide_angle.value_angle is not None
    assert disulfide_angle.value_angle > 0
    assert disulfide_angle.value_angle_esd is not None
    assert disulfide_angle.value_angle_esd > 1e-5

    added = False
    # FIXME this does not work in some situations
    for sym_pair in disulfide_sym_table.iterator():
      i_seq = self.cystein_sulphur_i_seqs[sym_pair.i_seq]
      j_seq = self.cystein_sulphur_i_seqs[sym_pair.j_seq]
      # add atoms to PDB link object
      # need to include sym. op.
      specific_origin_id = origin_ids.get_origin_id('SS BOND')
      self.pdb_link_records.setdefault("SSBOND", [])
      self.pdb_link_records["SSBOND"].append([self.pdb_atoms[i_seq],
                                              self.pdb_atoms[j_seq],
                                              sym_pair.rt_mx_ji,
                                            ]
                                            )
      bond_params_table.update(
        i_seq=i_seq,
        j_seq=j_seq,
        params=geometry_restraints.bond_params(
          distance_ideal=disulfide_bond.value_dist,
          weight=1/disulfide_bond.value_dist_esd**2,
          origin_id=specific_origin_id,
          ))
      bond_asu_table.add_pair(
        i_seq=i_seq,
        j_seq=j_seq,
        rt_mx_ji=sym_pair.rt_mx_ji)
      sym_str = str(sym_pair.rt_mx_ji)
      if not sym_str:
        sym_str = "."
      if (self.params.add_angle_and_dihedral_restraints_for_disulfides
          and sym_str == "x,y,z"):
        # Here we will add angles and torsions for disulfides because here
        # we do know i_seq, j_seq of linked atoms, and we do know that there
        # is no symmetry operation for this link
        def _get_ss_atom_pairs(atom_name):
          ss_atoms = [{}, {}]
          for i, seq in enumerate([i_seq, j_seq]):
            a = self.pdb_atoms[seq]
            ag = a.parent()
            rg = ag.parent()
            atom = ag.get_atom(atom_name)
            if atom: ss_atoms[i][ag.altloc]=atom
            else:
              for ag_t in rg.atom_groups():
                if ag_t.get_atom(atom_name):
                  ss_atoms[i][ag_t.altloc] = ag_t.get_atom(atom_name)
          return ss_atoms
        def genereate_ss_atom_lookup():
          cb_atoms = _get_ss_atom_pairs("CB")
          ca_atoms = _get_ss_atom_pairs("CA")
          altlocs = []
          for atoms in [ca_atoms, cb_atoms]:
            for residue in atoms:
              for altloc in residue:
                if altloc not in altlocs: altlocs.append(altloc)
          for altloc in altlocs:
            lookup = {
              "1SG" : i_seq,
              "2SG" : j_seq,
              }
            if altloc in ["", " "]:
              if altloc in ca_atoms[0]: lookup["1CA"]=ca_atoms[0][altloc].i_seq
              if altloc in ca_atoms[1]: lookup["2CA"]=ca_atoms[1][altloc].i_seq
              if altloc in cb_atoms[0]: lookup["1CB"]=cb_atoms[0][altloc].i_seq
              if altloc in cb_atoms[1]: lookup["2CB"]=cb_atoms[1][altloc].i_seq
            else:
              if altloc in ca_atoms[0]: lookup["1CA"]=ca_atoms[0][altloc].i_seq
              elif ""   in ca_atoms[0]: lookup["1CA"]=ca_atoms[0][""].i_seq
              if altloc in ca_atoms[1]: lookup["2CA"]=ca_atoms[1][altloc].i_seq
              elif ""   in ca_atoms[1]: lookup["2CA"]=ca_atoms[1][""].i_seq
              if altloc in cb_atoms[0]: lookup["1CB"]=cb_atoms[0][altloc].i_seq
              elif ""   in cb_atoms[0]: lookup["1CB"]=cb_atoms[0][""].i_seq
              if altloc in cb_atoms[1]: lookup["2CB"]=cb_atoms[1][altloc].i_seq
              elif ""   in cb_atoms[1]: lookup["2CB"]=cb_atoms[1][""].i_seq
            if len(lookup)==6:
              yield lookup

        for lookup in genereate_ss_atom_lookup():
          # move checking of having all atoms to generator
          angle_weight = 1/disulfide_angle.value_angle_esd**2
          proxy = geometry_restraints.angle_proxy(
            i_seqs=[lookup["1CB"],i_seq,j_seq],
            angle_ideal=disulfide_angle.value_angle,
            weight=angle_weight,
            origin_id=specific_origin_id)
          self.geometry_proxy_registries.angle.add_if_not_duplicated(proxy=proxy)
          proxy = geometry_restraints.angle_proxy(
            i_seqs=[i_seq,j_seq,lookup["2CB"]],
            angle_ideal=disulfide_angle.value_angle,
            weight=angle_weight,
            origin_id=specific_origin_id)
          self.geometry_proxy_registries.angle.add_if_not_duplicated(proxy=proxy)
          if 0:
            indent=14
            print("      Atoms : %s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s" % (
              self.pdb_atoms[lookup["1CA"]].quote(),
              ' '*indent,
              self.pdb_atoms[lookup["1CB"]].quote(),
              ' '*indent,
              self.pdb_atoms[lookup["1SG"]].quote(),
              ' '*indent,
              self.pdb_atoms[lookup["2SG"]].quote(),
              ' '*indent,
              self.pdb_atoms[lookup["2CB"]].quote(),
              ' '*indent,
              self.pdb_atoms[lookup["2CA"]].quote(),
              ), file=log)
          for disulfide_torsion in disulfide_torsions:
            assert disulfide_torsion.value_angle is not None
            assert disulfide_torsion.value_angle_esd is not None
            assert disulfide_torsion.value_angle_esd > 1e-5
            assert disulfide_torsion.period is not None
            assert disulfide_torsion.period >= 0
            alt_value_angle = None
            if (disulfide_torsion.alt_value_angle is not None and
                disulfide_torsion.alt_value_angle != ''):
              try:
                alt_value_angle = [float(_d) for _d in disulfide_torsion.alt_value_angle.split(",")]
              except ValueError as AttributeError:
                raise Sorry("Wrong format of alt_value_angle in SS bond in cif file")
            i_seqs = []
            for nwm in range(1,5):
              key = "%s%s" % (getattr(disulfide_torsion, "atom_%d_comp_id" % nwm),
                             getattr(disulfide_torsion, "atom_id_%s" % nwm),
                )
              i_seqs.append(lookup[key])
            proxy = geometry_restraints.dihedral_proxy(
              i_seqs=i_seqs,
              angle_ideal=disulfide_torsion.value_angle,
              weight=1/disulfide_torsion.value_angle_esd**2,
              periodicity=disulfide_torsion.period,
              alt_angle_ideals=alt_value_angle,
              origin_id=specific_origin_id)
            self.geometry_proxy_registries.dihedral.add_if_not_duplicated(proxy=proxy)
      if disulfide_cif_loop is not None:
        disulfide_cif_loop.add_row(("SS",
                                    self.pdb_atoms[i_seq].pdb_label_columns(),
                                    self.pdb_atoms[j_seq].pdb_label_columns(),
                                    sym_str,
                                    ))
        # added = True
    # if added:
    #   self._cif.cif["link_SS"] = disulfide_link.as_cif_block()
    #
    # ====================== End of disulfides ========================
    #
    if (processed_edits is not None):
      for proxy in processed_edits.bond_sym_proxies:
        if (proxy.weight <= 0): continue
        i_seq, j_seq = proxy.i_seqs
        # print (dir(bond_params_table))
        # STOP()
        bond_params_table.update(i_seq=i_seq, j_seq=j_seq, params=proxy)
        bond_asu_table.add_pair(
          i_seq=i_seq,
          j_seq=j_seq,
          rt_mx_ji=proxy.rt_mx_ji)
      not_added_proxies = []
      for proxy in processed_edits.angle_proxies:
        added = self.geometry_proxy_registries.angle.add_if_not_duplicated(proxy=proxy)
        if not added:
          not_added_proxies.append(proxy)
      for proxy in processed_edits.dihedral_proxies:
        added = self.geometry_proxy_registries.dihedral.add_if_not_duplicated(proxy=proxy)
        if not added:
          not_added_proxies.append(proxy)
      for proxy in processed_edits.planarity_proxies:
        added = self.geometry_proxy_registries.planarity.add_if_not_duplicated(proxy=proxy)
        if not added:
          not_added_proxies.append(proxy)
      for proxy in processed_edits.parallelity_proxies:
        added = self.geometry_proxy_registries.parallelity.add_if_not_duplicated(proxy=proxy)
        if not added:
          not_added_proxies.append(proxy)
      if len(not_added_proxies) > 0:
        raise Sorry("Some restraints were not added because they are already present.")

    if params_edits and params_edits.angle:
      processed_edits = self.process_geometry_restraints_edits(
        params=params_edits, log=log, second_pass=True)
      max_bond_distance = max(max_bond_distance,
        processed_edits.bond_distance_model_max)
    #
    al_params = self.params.automatic_linking
    any_links = False
    link_attrs = ["link_metals",
                  "link_residues",
                  "link_carbohydrates",
                  "link_amino_acid_rna_dna",
                  "link_ligands",
                  'link_small_molecules',
                  ]
    if al_params.link_all:
      for attr in link_attrs:
        setattr(al_params, attr, True)
        any_links = True
    if al_params.link_none:
      for attr in link_attrs:
        setattr(al_params, attr, False)
        any_links = False
    for attr in link_attrs:
      if getattr(al_params, attr, False):
        any_links = True
        break
    #
    selection_cache = self.pdb_hierarchy.atom_selection_cache()
    exclude_selections = []
    for efal in self.params.exclude_from_automatic_linking:
      if efal.selection_1 is None: continue
      t_selection = flex.size_t()
      t_selection = self.pdb_hierarchy.atom_selection_cache().selection(efal.selection_1).iselection()
      exclude_selections.append([t_selection])
      if efal.selection_2 is not None:
        t_selection = flex.size_t()
        t_selection = self.pdb_hierarchy.atom_selection_cache().selection(efal.selection_2).iselection()
        exclude_selections[-1].append(t_selection)
      else:
        exclude_selections[-1].append(None)
    include_selections = []
    for iial in self.params.include_in_automatic_linking:
      if iial.selection_1 is None: continue
      t_selection = flex.size_t()
      t_selection = self.pdb_hierarchy.atom_selection_cache().selection(iial.selection_1).iselection()
      include_selections.append([t_selection])
      if iial.selection_2 is not None:
        t_selection = flex.size_t()
        t_selection = self.pdb_hierarchy.atom_selection_cache().selection(iial.selection_2).iselection()
        include_selections[-1].append(t_selection)
      else:
        include_selections[-1].append(None)
      include_selections[-1].append(iial.bond_cutoff)
    #
    if any_links:
      self.process_nonbonded_for_links(
        bond_params_table,
        bond_asu_table,
        self.geometry_proxy_registries,
        link_metals               = al_params.link_metals,
        link_residues             = al_params.link_residues,
        link_carbohydrates        = al_params.link_carbohydrates,
        link_amino_acid_rna_dna   = al_params.link_amino_acid_rna_dna,
        link_ligands              = al_params.link_ligands,
        link_small_molecules      = al_params.link_small_molecules,
        amino_acid_bond_cutoff    = al_params.amino_acid_bond_cutoff,
        metal_coordination_cutoff = al_params.metal_coordination_cutoff,
        carbohydrate_bond_cutoff  = al_params.carbohydrate_bond_cutoff,
        inter_residue_bond_cutoff = al_params.inter_residue_bond_cutoff,
        ligand_bond_cutoff        = al_params.ligand_bond_cutoff,
        small_molecule_bond_cutoff= al_params.small_molecule_bond_cutoff,
        second_row_buffer         = al_params.buffer_for_second_row_elements,
        exclude_selections        = exclude_selections,
        include_selections        = include_selections,
        log=log,
        )
    self.geometry_proxy_registries.discard_tables()
    self.scattering_type_registry.discard_tables()
    self.nonbonded_energy_type_registry.discard_tables()
    self.geometry_proxy_registries.bond_simple.proxies = None # free memory
    #
    shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
      pair_asu_table=bond_asu_table,
      max_shell=3)
    if (custom_nonbonded_exclusions is not None):
      self.process_custom_nonbonded_exclusions(
        log=log,
        exclude_pair_indices=custom_nonbonded_exclusions,
        shell_asu_tables=shell_asu_tables)
    shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
      for shell_asu_table in shell_asu_tables]
    nonbonded_params = ener_lib_as_nonbonded_params(
      ener_lib=ener_lib,
      assume_hydrogens_all_missing=assume_hydrogens_all_missing,
      factor_1_4_interactions=self.params.vdw_1_4_factor,
      default_distance=self.params.default_vdw_distance,
      minimum_distance=self.params.min_vdw_distance,
      const_shrink_donor_acceptor=self.params.const_shrink_donor_acceptor)
    if(self.params.nonbonded_weight is None):
      nonbonded_weight = 100 # c_rep in prolsq repulsion function
      if(assume_hydrogens_all_missing):
        nonbonded_weight = 100
    else:
      nonbonded_weight = self.params.nonbonded_weight
    #
    result = geometry_restraints.manager.manager(
      crystal_symmetry=self.special_position_settings,
      model_indices=self.model_indices,
      conformer_indices=self.conformer_indices,
      sym_excl_indices=self.sym_excl_indices,
      donor_acceptor_excl_groups=self.donor_acceptor_excl_groups,
      site_symmetry_table=self.site_symmetry_table(),
      bond_params_table=bond_params_table,
      shell_sym_tables=shell_sym_tables,
      nonbonded_params=nonbonded_params,
      nonbonded_types=self.nonbonded_energy_type_registry.symbols,
      nonbonded_charges=self.nonbonded_energy_type_registry.charges,
      nonbonded_function=geometry_restraints.prolsq_repulsion_function(
        c_rep=nonbonded_weight),
      nonbonded_distance_cutoff=self.params.nonbonded_distance_cutoff,
      nonbonded_buffer=self.params.nonbonded_buffer,
      angle_proxies=self.geometry_proxy_registries.angle.proxies,
      dihedral_proxies=self.geometry_proxy_registries.dihedral.proxies,
      chirality_proxies=self.geometry_proxy_registries.chirality.proxies,
      planarity_proxies=self.geometry_proxy_registries.planarity.proxies,
      parallelity_proxies=self.geometry_proxy_registries.parallelity.proxies,
      ramachandran_manager=None,
      external_energy_function=external_energy_function,
      max_reasonable_bond_distance=self.params.max_reasonable_bond_distance,
      plain_pairs_radius=plain_pairs_radius,
      log=log)
    if (params_remove is not None):
      self.process_geometry_restraints_remove(
        params=params_remove, geometry_restraints_manager=result)
    self.time_building_geometry_restraints_manager = timer.elapsed()
    restraints_source = "GeoStd + Monomer Library"
    use_cdl = self.params.restraints_library.cdl
    cis_pro_eh99 = self.params.restraints_library.cis_pro_eh99
    if (use_cdl is Auto):
      use_cdl = self.pdb_inp.used_cdl_restraints()
      if (use_cdl):
        print("  Switching to conformation-dependent library", file=log)
    if use_cdl:
      restraints_source += ' + %s' % mmtbx.conformation_dependent_library.cdl_database.version
      from mmtbx.conformation_dependent_library.cdl_setup import setup_restraints
      from mmtbx.conformation_dependent_library import update_restraints
      from libtbx import utils
      t0=time.time()
      cdl_proxies=setup_restraints(result)
      update_restraints(
        self.pdb_hierarchy,
        result,
        cdl_proxies=cdl_proxies,
        cis_pro_eh99=cis_pro_eh99,
        cdl_svl=self.params.restraints_library.cdl_svl,
        log=log,
        verbose=True,
        )
      self.use_cdl = True # what is this used for???
      cdl_time = time.time()-t0
      print("""\
  Conformation dependent library (CDL) restraints added in %0.1f %sseconds
  """ % utils.greek_time(cdl_time), file=log)
    #
    # need autodetect code
    #
    use_omega_cdl = self.params.restraints_library.omega_cdl
    if (use_omega_cdl is Auto):
      use_omega_cdl = self.pdb_inp.used_omega_cdl_restraints()
      if (use_omega_cdl):
        print("  Switching to omega-CDL", file=log)
    if use_omega_cdl:
      restraints_source += ' + omega-cdl'
      from mmtbx.conformation_dependent_library.omega import setup_restraints
      from mmtbx.conformation_dependent_library.omega import update_restraints
      from libtbx import utils
      t0=time.time()
      cdl_proxies=setup_restraints(result)
      update_restraints(
        self.pdb_hierarchy,
        result,
        cdl_proxies=cdl_proxies,
        log=log,
        verbose=True,
        )
      self.use_omega_cdl = True
      cdl_time = time.time()-t0
      print("""\
  omega-Conformation dependent library (o-CDL) restraints added in %0.1f %sseconds
  """ % utils.greek_time(cdl_time), file=log)
    #
    if getattr(self.params.restraints_library, "rdl", False):
      from mmtbx.conformation_dependent_library import rotamers
      from libtbx import utils
      t0=time.time()
      rdl_proxies=None #rotamers.setup_restraints(result)
      rotamers.update_restraints(
        self.pdb_hierarchy,
        result,
        #current_geometry=model.xray_structure,
        rdl_proxies=rdl_proxies,
        data_version="8000",
        log=log,
        verbose=False,
        )
      rdl_time = time.time()-t0
      print("""\
  Rotamer dependent library (RDL) restraints added in %0.1f %sseconds
  """ % utils.greek_time(rdl_time), file=log)
    if getattr(self.params.restraints_library, "hpdl", False):
      from mmtbx.conformation_dependent_library import histidines
      from libtbx import utils
      t0=time.time()
      histidines.update_restraints(
        self.pdb_hierarchy,
        result,
        log=log,
        verbose=True,
        )
      hpr_time = time.time()-t0
      print("""\
  Histidine protonation dependent restraints added in %0.1f %sseconds
  """ % utils.greek_time(hpr_time), file=log)
    #
    if self.pdb_inp and self.pdb_inp.used_amber_restraints():
      restraints_source = 'Amber'
      self.params.use_neutron_distances = True
    #
    result.set_source(source = restraints_source)
    return result

  def extract_xray_structure(self, unknown_scattering_type_substitute = "?"):
    from cctbx import xray
    from cctbx import adptbx
    from itertools import count
    assert self.special_position_settings is not None
    result = xray.structure(
      special_position_settings=self.special_position_settings)
    sites_frac = result.unit_cell().fractionalize(
      sites_cart=self.sites_cart_exact())
    site_symmetry_table = self.site_symmetry_table()
    if (site_symmetry_table is not None):
      assert site_symmetry_table.indices().size() == sites_frac.size()
    scattering_types = self.scattering_type_registry.symbols
    if (scattering_types is None):
      scattering_types = self.get_element_symbols(strip_symbols=True)
    site_symmetry_ops = None
    for i_seq,atom,site_frac,scattering_type in zip(
          count(),
          self.pdb_atoms,
          sites_frac,
          scattering_types):
      assert atom.i_seq == i_seq
      from cctbx.eltbx.xray_scattering import get_standard_label
      scattering_type = get_standard_label(
        label=scattering_type, exact=True, optional=True)
      if (scattering_type is not None):
        charge = atom.charge_tidy(strip=True)
        if (charge is not None and len(charge) != 0):
          scattering_type_with_charge = get_standard_label(
            label=scattering_type+charge, exact=False, optional=True)
          if (scattering_type_with_charge is not None):
            scattering_type = scattering_type_with_charge
      if (scattering_type is None):
        if (unknown_scattering_type_substitute is None):
          raise RuntimeError("Unknown scattering type: %s" %
            atom.format_atom_record(cut_after_label_columns=True))
        scattering_type = unknown_scattering_type_substitute
      if (not atom.uij_is_defined()):
        u = adptbx.b_as_u(atom.b)
      else:
        u = adptbx.u_cart_as_u_star(result.unit_cell(), atom.uij)
      if (site_symmetry_table is not None):
        site_symmetry_ops = site_symmetry_table.get(i_seq=i_seq)
      result.add_scatterer(
        scatterer=xray.scatterer(
          label=atom.id_str(),
          site=site_frac,
          u=u,
          occupancy=atom.occ,
          scattering_type=scattering_type,
          fp=atom.fp,
          fdp=atom.fdp),
        site_symmetry_ops=site_symmetry_ops)
    return result

def show_residue_groups(residue_groups, log, prefix, max_items):
  if (len(residue_groups) == max_items + 1):
    max_items += 1
  for rg in residue_groups[:max_items]:
    if (rg.unique_resnames().size() == 1):
      print(prefix+"residue:", file=log)
    else:
      print(prefix+"residue group:", file=log)
    rg_atoms = rg.atoms()
    def show_atom(i):
      a = rg_atoms[i]
      print(prefix+"  %s occ=%.2f" % (a.id_str(), a.occ), file=log)
    show_atom(0)
    n = rg_atoms.size()
    if (n > 3):
      print(prefix+"  ... (%d atoms not shown)" % (n-2), file=log)
    elif (n == 3):
      show_atom(1)
    if (n > 1):
      show_atom(-1)
  n = len(residue_groups) - max_items
  if (n > 0):
    print(prefix+"... (remaining %d not shown)" % n, file=log)

class process(object):

  def __init__(self,
        mon_lib_srv,
        ener_lib,
        params=None,
        file_name=None,
        raw_records=None,
        pdb_inp=None,
        pdb_hierarchy=None,
        atom_selection_string=None,
        strict_conflict_handling=False,
        special_position_settings=None,
        crystal_symmetry=None,
        force_symmetry=False,
        substitute_non_crystallographic_unit_cell_if_necessary=False,
        keep_monomer_mappings=False,
        max_atoms=None,
        log=None,
        carbohydrate_callback=None,
        restraints_loading_flags=None):
    self.mon_lib_srv = mon_lib_srv
    self.ener_lib = ener_lib
    self.log = log
    nul = StringIO()
    # allow pdb_inp and pdb_hierarchy present simultaneously
    assert [file_name, raw_records, pdb_inp, pdb_hierarchy].count(None) >= 2
    if log is None:
      self.log = nul
    if file_name is not None:
      assert pdb_inp is None
      pdb_inp = iotbx.pdb.input(file_name=file_name)
    if raw_records is not None:
      assert pdb_inp is None
      assert pdb_hierarchy is None
      assert file_name is None
      if (isinstance(raw_records, str)):
        raw_records = flex.split_lines(raw_records)
      elif (not isinstance(raw_records, flex.std_string)):
        raw_records = flex.std_string(raw_records)
      pdb_inp = pdb.input(source_info=None, lines=raw_records)

    self.ss_manager = None
    self.ss_torsion_restraints = None
    self.all_chain_proxies = build_all_chain_proxies(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=params,
      pdb_inp=pdb_inp,
      pdb_hierarchy=pdb_hierarchy,
      atom_selection_string=atom_selection_string,
      special_position_settings=special_position_settings,
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      substitute_non_crystallographic_unit_cell_if_necessary
        =substitute_non_crystallographic_unit_cell_if_necessary,
      strict_conflict_handling=strict_conflict_handling,
      keep_monomer_mappings=keep_monomer_mappings,
      max_atoms=max_atoms,
      log=log,
      carbohydrate_callback=carbohydrate_callback,
      restraints_loading_flags=restraints_loading_flags)
    if (log is not None
        and self.all_chain_proxies.time_building_chain_proxies is not None):
      print("  Time building chain proxies: %.2f, per 1000 atoms: %.2f" % (
          self.all_chain_proxies.time_building_chain_proxies,
          self.all_chain_proxies.time_building_chain_proxies * 1000
            / max(1,self.all_chain_proxies.pdb_atoms.size())), file=log)

    self._geometry_restraints_manager = None
    self._xray_structure = None
    # Find NCS
    self.ncs_obj = None
    if(len(list(self.all_chain_proxies.pdb_hierarchy.models())) == 1 and
       self.all_chain_proxies.params.ncs_search.enabled):
      self.ncs_obj = self.search_for_ncs(
        hierarchy = self.all_chain_proxies.pdb_hierarchy)
    # attempt to fix pH problems
    if self.all_chain_proxies.nonbonded_energy_type_registry.n_unknown_type_symbols():
      from mmtbx.conformation_dependent_library import pH_dependent_restraints
      unknown_atoms = \
        self.all_chain_proxies.nonbonded_energy_type_registry.get_unknown_atoms(
          self.all_chain_proxies.pdb_hierarchy.atoms(),
          return_iseqs=True,
        )
      rc = pH_dependent_restraints.adjust_geometry_proxies_registeries(
        self.all_chain_proxies.pdb_hierarchy,
        self.all_chain_proxies.geometry_proxy_registries,
        unknown_atoms,
        )
      # update the nonbonded energy of "new" atoms - should do all
      # needed for book keeping of successful restraints matching
      for atom_i_seq, item in rc.items():
        if atom_i_seq in unknown_atoms:
          self.all_chain_proxies.nonbonded_energy_type_registry.symbols[atom_i_seq] = \
              item.type_energy

  def geometry_restraints_manager(self,
        plain_pairs_radius=None,
        params_edits=None,
        params_remove=None,
        custom_nonbonded_exclusions=None,
        assume_hydrogens_all_missing=True,
        show_energies=True,
        hard_minimum_bond_distance_model=0.001,
        external_energy_function=None,
        den_manager=None,
        ):
    if (    self.all_chain_proxies.sites_cart is not None
        and self.all_chain_proxies.special_position_settings is not None
        and self._geometry_restraints_manager is None):
      self._geometry_restraints_manager \
        = self.all_chain_proxies.construct_geometry_restraints_manager(
            ener_lib=self.ener_lib,
            disulfide_link=self.mon_lib_srv.link_link_id_dict["SS"],
            plain_pairs_radius=plain_pairs_radius,
            params_edits=params_edits,
            params_remove=params_remove,
            custom_nonbonded_exclusions=custom_nonbonded_exclusions,
            assume_hydrogens_all_missing=assume_hydrogens_all_missing,
            external_energy_function=external_energy_function,
            log=self.log)

      #for key, item in self.all_chain_proxies.pdb_link_records.items():
      #  print key
      #  for atom1, atom2, sym_op in item:
      #    print atom1.quote(), atom2.quote(), sym_op

      # initializing grm
      self._geometry_restraints_manager.pair_proxies(
          sites_cart=self.all_chain_proxies.sites_cart)

      # improved metal coordination
      automatic_linking = self.all_chain_proxies.params.automatic_linking
      if self.all_chain_proxies.params.restraints_library.mcl:
        from mmtbx.conformation_dependent_library import mcl
        mcl.update(self._geometry_restraints_manager,
                   self.all_chain_proxies.pdb_hierarchy,
                   log=self.log,
                  )

      # Here we are going to add another needed restraints.
      # Ramachandran restraints
      ramachandran_manager = None
      pep_link_params = self.all_chain_proxies.params.peptide_link
      rama_params = None
      if hasattr(self.all_chain_proxies.params, 'ramachandran_plot_restraints'):
        rama_params = self.all_chain_proxies.params.ramachandran_plot_restraints
      if (pep_link_params.ramachandran_restraints or
          (rama_params and rama_params.enabled)):
        if (not pep_link_params.discard_psi_phi):
          # Not sure anymore why this is necessary
          raise Sorry("You may not use Ramachandran restraints when "+
            "discard_psi_phi=False.")
        params_to_pass = rama_params
        if pep_link_params.ramachandran_restraints:
          params_to_pass = pep_link_params
        # print ('params_to_pass', type(params_to_pass))
        ramachandran_manager = ramachandran.ramachandran_manager(
            pdb_hierarchy=self.all_chain_proxies.pdb_hierarchy,
            params=params_to_pass,
            log=self.log)
        self._geometry_restraints_manager.set_ramachandran_restraints(
            ramachandran_manager)

      # C-beta restraints
      if self.all_chain_proxies.params.c_beta_restraints:
        print("  Adding C-beta torsion restraints...", file=self.log)
        from mmtbx.geometry_restraints import c_beta
        c_beta_torsion_proxies, c_beta_skipped = \
            c_beta.get_c_beta_torsion_proxies(self.all_chain_proxies.pdb_hierarchy)
        n_c_beta_restraints = len(c_beta_torsion_proxies)
        if n_c_beta_restraints > 0:
          self._geometry_restraints_manager.add_dihedrals_in_place(
              c_beta_torsion_proxies,
              check_for_duplicates=False)
        outl = ""
        for key, item in c_beta_skipped.items():
          if key=="-ve":
            outl += "      Input volumes are d-peptide like\n"
          elif key=="d-peptide":
            outl += "      Input residue name is d-peptide\n"
          for cb_atom in item:
            outl += "        %s\n"% cb_atom.id_str()
        if outl:
          print("    Skipped\n%s" % outl[:-1], file=self.log)
        print("  Number of C-beta restraints generated: ",\
            n_c_beta_restraints, file=self.log)
        print(file=self.log)

      # Reference coordinate restraints
      if self.all_chain_proxies.params.reference_coordinate_restraints.enabled:
        rcr = self.all_chain_proxies.params.reference_coordinate_restraints
        self._geometry_restraints_manager.\
            add_reference_coordinate_restraints_in_place(
                all_chain_proxies=self.all_chain_proxies,
                # pdb_hierarchy=self.all_chain_proxies.pdb_hierarchy,
                selection=rcr.selection,
                exclude_outliers=rcr.exclude_outliers,
                sigma=rcr.sigma,
                limit=rcr.limit,
                top_out=rcr.top_out)
        n_rcr = self._geometry_restraints_manager.\
            get_n_reference_coordinate_proxies()
        print("  Number of reference coordinate restraints generated:",\
           n_rcr, file=self.log)

      # DEN manager
      self._geometry_restraints_manager.adopt_den_manager(den_manager)

      # Secondary structure restraints:
      # Proteins first
      # then nucleic acids restraints.
      # from mmtbx.secondary_structure import sec_str_master_phil
      # def_ss_params = sec_str_master_phil.fetch().extract()
      # def_ss_params.secondary_structure = self.all_chain_proxies.params.secondary_structure
      ss_params = self.all_chain_proxies.params.secondary_structure
      if ss_params.enabled:
        from mmtbx.secondary_structure import manager
        t0=time.time()
        print("  Finding SS restraints...", file=self.log)
        self.ss_manager = manager(
            pdb_hierarchy=self.all_chain_proxies.pdb_hierarchy,
            geometry_restraints_manager=self._geometry_restraints_manager,
            sec_str_from_pdb_file=self.all_chain_proxies.extract_secondary_structure(),
            params=ss_params,
            mon_lib_srv=self.mon_lib_srv,
            verbose=-1,
            log=self.log)
        t1=time.time()
        print("    Time for finding SS restraints: %.2f" % (t1-t0), file=self.log)
        print("  Creating SS restraints...", file=self.log)

        self._geometry_restraints_manager.set_secondary_structure_restraints(
            ss_manager=self.ss_manager,
            hierarchy=self.all_chain_proxies.pdb_hierarchy,
            log=self.log)

        t3=time.time()
        print("  Total time for adding SS restraints: %.2f" % (t3-t1), file=self.log)
        print(file=self.log)
      if (self.log is not None):
        print("  Time building geometry restraints manager: %.2f seconds" % (
            self.all_chain_proxies.time_building_geometry_restraints_manager), file=self.log)
        print(file=self.log)
        def note_geo():
          print("""\
  NOTE: a complete listing of the restraints can be obtained by requesting
        output of .geo file.""", file=self.log)
        note_geo()
        print(file=self.log)
        flush_log(self.log)
        site_labels = [atom.id_str()
                       for atom in self.all_chain_proxies.pdb_atoms]
        pair_proxies = self._geometry_restraints_manager.pair_proxies(
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          site_labels=site_labels)
        params = self.all_chain_proxies.params
        pair_proxies.bond_proxies.show_histogram_of_model_distances(
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          n_slots=params.show_histogram_slots.bond_lengths,
          f=self.log,
          prefix="  ",
          origin_id=origin_ids.get_origin_id('covalent geometry'))
        smallest_distance_model = \
          pair_proxies.bond_proxies.show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items.bond_restraints_sorted_by_residual,
            origin_id=origin_ids.get_origin_id('covalent geometry'))
        if (    smallest_distance_model is not None
            and hard_minimum_bond_distance_model is not None
            and smallest_distance_model < hard_minimum_bond_distance_model):
          raise Sorry("""Bond restraint model distance < %.6g:
  Please inspect the output above and correct the input model file.""" % (
            hard_minimum_bond_distance_model))
        print(file=self.log)
        self._geometry_restraints_manager.angle_proxies \
          .show_histogram_of_deltas(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            n_slots=params.show_histogram_slots
              .bond_angle_deviations_from_ideal,
            f=self.log,
            prefix="  ",
            origin_id=origin_ids.get_origin_id('covalent geometry'))
        self._geometry_restraints_manager.angle_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .bond_angle_restraints_sorted_by_residual,
            origin_id=origin_ids.get_origin_id('covalent geometry'))
        print(file=self.log)
        self._geometry_restraints_manager.dihedral_proxies \
          .show_histogram_of_deltas(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            n_slots=params.show_histogram_slots
              .dihedral_angle_deviations_from_ideal,
            f=self.log,
            prefix="  ")
        self._geometry_restraints_manager.dihedral_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .dihedral_angle_restraints_sorted_by_residual)
        print(file=self.log)
        self._geometry_restraints_manager.chirality_proxies \
          .show_histogram_of_deltas(
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            n_slots=params.show_histogram_slots
              .chiral_volume_deviations_from_ideal,
            f=self.log,
            prefix="  ")
        self._geometry_restraints_manager.chirality_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .chirality_restraints_sorted_by_residual)
        print(file=self.log)
        self._geometry_restraints_manager.planarity_proxies \
          .show_sorted(
            by_value="residual",
            sites_cart=self.all_chain_proxies.sites_cart_exact(),
            site_labels=site_labels,
            f=self.log,
            prefix="  ",
            max_items=params.show_max_items
              .planarity_restraints_sorted_by_residual)
        print(file=self.log)
        pair_proxies.nonbonded_proxies.show_histogram_of_model_distances(
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          n_slots=params.show_histogram_slots.nonbonded_interaction_distances,
          f=self.log,
          prefix="  ")
        pair_proxies.nonbonded_proxies.show_sorted(
          by_value="delta",
          sites_cart=self.all_chain_proxies.sites_cart_exact(),
          site_labels=site_labels,
          f=self.log,
          prefix="  ",
          max_items=params.show_max_items
            .nonbonded_interactions_sorted_by_model_distance)
        print(file=self.log)
        note_geo()
        flush_log(self.log)
        if (show_energies):
          print(file=self.log)
          timer = user_plus_sys_time()
          energies = self._geometry_restraints_manager.energies_sites(
            sites_cart=self.all_chain_proxies.sites_cart_exact())
          energies.show(f=self.log, prefix="  ")
          print("  Time first energy calculation" \
                             " (mainly nonbonded setup): %.2f" % (
            timer.elapsed()), file=self.log)
          flush_log(self.log)
    return self._geometry_restraints_manager


  def clash_guard(self,
                  hard_minimum_nonbonded_distance=0.001,
                  nonbonded_distance_threshold=0.5,
                  new_sites_cart=None):
    params = self.all_chain_proxies.params.clash_guard
    if nonbonded_distance_threshold != 0.5:
      # WHY is this here???
      # only if nonbonded_distance_threshold is not the default
      # use it, otherwise use the parameters
      # passed by self.all_chain_proxies.params.clash_guard
      params.nonbonded_distance_threshold = nonbonded_distance_threshold
    if (params.nonbonded_distance_threshold is None): return
    geo = self._geometry_restraints_manager
    if geo is None:
      return None
    # This is done for phenix.refine when run with shaking coordinates
    if new_sites_cart is not None:
      geo.pair_proxies(sites_cart=new_sites_cart)
    n_below_threshold = (
      geo.nonbonded_model_distances() < params.nonbonded_distance_threshold) \
        .count(True)
    if ((params.max_number_of_distances_below_threshold is not None
         and n_below_threshold
               > params.max_number_of_distances_below_threshold)
        or
        (params.max_fraction_of_distances_below_threshold is not None
         and n_below_threshold
               / max(1,geo.sites_cart_used_for_pair_proxies().size())
                 > params.max_fraction_of_distances_below_threshold)):
      phil_path = params.__phil_path__()
      return """%s failure:
  Number of nonbonded interaction distances < %.6g: %d
    Please inspect the histogram of nonbonded interaction distances above
    (or in the log window) for atoms placed too close to each other.
    To disable this error, run the same command again with the following
    additional argument:
      %s.nonbonded_distance_threshold=None""" % (
        phil_path, params.nonbonded_distance_threshold,
        n_below_threshold, phil_path)
    #
    n_below_hard_minimum_nonbonded_distance = (
      geo.nonbonded_model_distances() < hard_minimum_nonbonded_distance) \
        .count(True)
    if (n_below_hard_minimum_nonbonded_distance != 0 and
      params.nonbonded_distance_threshold >=0):
      return """Number of nonbonded interaction distances < %.6g: %d
  Please inspect the output above (or in the log window) and correct the input
  model file.""" % (
        hard_minimum_nonbonded_distance,
        n_below_hard_minimum_nonbonded_distance)
    return None

  def xray_structure(self, show_summary = True):
    log = self.log
    if (    self.all_chain_proxies.sites_cart is not None
        and self.all_chain_proxies.special_position_settings is not None
        and self._xray_structure is None):
      self._xray_structure = self.all_chain_proxies.extract_xray_structure()
      self._xray_structure.scattering_type_registry(
        types_without_a_scattering_contribution=["?"])
      # Stop if polymer crosses symmetry element
      get_class = iotbx.pdb.common_residue_names_get_class
      lines = []
      for i_sp in self._xray_structure.special_position_indices():
        atom = self.all_chain_proxies.pdb_hierarchy.atoms()[i_sp]
        resname = atom.parent().resname
        if(get_class(name=resname) == "common_amino_acid" or
           get_class(name=resname) == "common_rna_dna"):
          lines.append(atom.format_atom_record())
      if(len(lines)>0 and not
         self.all_chain_proxies.params.allow_polymer_cross_special_position):
        msg="Polymer crosses special position element:\n%s"%("\n".join(lines))
        msg+="\nUse 'allow_polymer_cross_special_position=True' to keep going."
        raise Sorry(msg)
      #
      if (log is not None and show_summary):
        self._xray_structure.show_summary(f = log, prefix="  ")
        self._xray_structure.show_special_position_shifts(
          sites_cart_original=self.all_chain_proxies.sites_cart,
          out= log, prefix="  ")
        self._xray_structure.scattering_type_registry().show(
          show_gaussians=False, out = log, prefix="  ")
        flush_log(log)
    return self._xray_structure

  def search_for_ncs(self, hierarchy):
    ncs_phil_groups = self.all_chain_proxies.params.ncs_group
    # This function may alter pdb_hierarchy, e.g. when chain_id is blank
    # it substitutes them with "A". Therefore deep_copy() is necessary.
    # ! Should no longer be the case.
    new_h = hierarchy.deep_copy()
    new_h.atoms().reset_i_seq()
    ncs_obj = iotbx.ncs.input(
      ncs_phil_groups             = ncs_phil_groups,
      hierarchy                   = new_h,
      params                      = self.all_chain_proxies.params.ncs_search,
      log                         = self.log)
    print("Found NCS groups:", file=self.log)
    print(ncs_obj.print_ncs_phil_param(), file=self.log)
    return ncs_obj

def run(
      args,
      params=None,
      strict_conflict_handling=True,
      substitute_non_crystallographic_unit_cell_if_necessary=False,
      return_all_processed_pdb_files=False,
      max_atoms=None,
      assume_hydrogens_all_missing=True,
      hard_minimum_nonbonded_distance=0.001,
      nonbonded_distance_threshold=0.5,
      log=None):
  if (log is None): log = sys.stdout
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  pdb_file_names = []
  cif_file_names = []
  for arg in args:
    if (not os.path.isfile(arg)):
      raise Sorry("No such file: %s" % show_string(arg))
    if (iotbx.pdb.is_pdb_file(file_name=arg)):
      pdb_file_names.append(arg)
    else:
      try:
        cif_object = server.read_cif(file_name=arg)
      except KeyboardInterrupt: raise
      except Exception:
        raise Sorry("Unknown file format: %s" % show_string(arg))
      else:
        print("Processing CIF file: %s" % show_string(arg), file=log)
        for srv in [mon_lib_srv, ener_lib]:
          srv.process_cif_object(cif_object=cif_object, file_name=arg)
  all_processed_pdb_files = []
  for file_name in pdb_file_names:
    processed_pdb_file = process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=params,
      file_name=file_name,
      strict_conflict_handling=strict_conflict_handling,
      substitute_non_crystallographic_unit_cell_if_necessary
        =substitute_non_crystallographic_unit_cell_if_necessary,
      max_atoms=max_atoms,
      log=log)
    processed_pdb_file.geometry_restraints_manager(
      assume_hydrogens_all_missing=assume_hydrogens_all_missing)
    err_msg = processed_pdb_file.clash_guard(
        hard_minimum_nonbonded_distance=hard_minimum_nonbonded_distance,
        nonbonded_distance_threshold=nonbonded_distance_threshold)
    if err_msg is not None:
      raise Sorry(err_msg)
    processed_pdb_file.xray_structure()
    if (return_all_processed_pdb_files):
      all_processed_pdb_files.append(processed_pdb_file)
  if (return_all_processed_pdb_files):
    return all_processed_pdb_files
  if (len(pdb_file_names) > 0):
    return processed_pdb_file
  return None

if (__name__ == "__main__"):
  run(sys.argv[1:])
