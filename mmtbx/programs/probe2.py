from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import null_out
from libtbx import group_args, phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry
import mmtbx
import mmtbx_probe_ext as probeExt
from mmtbx.probe import Helpers, AtomTypes
from scitbx.array_family import flex
from iotbx import pdb
from iotbx.pdb import common_residue_names_get_class

master_phil_str = '''
source_selection = None
  .type = str
  .help = Source selection description

target_selection = None
  .type = str
  .help = Target selection description ('=' means same as source)

use_neutron_distances = False
  .type = bool
  .help = Use neutron distances (-nuclear in probe)

approach = *self both once surface count_atoms
  .type = choice
  .help = self (src -> src) both (src <=> targ) once (src -> targ) surface (VdW surface) count_atoms (count atoms)

include_mainchain_mainchain = True
  .type = bool
  .help = Include mainchain -> mainchain interactions (-mc in probe)

include_het = True
  .type = bool
  .help = Include non-water HET group interactions (-het in probe)

include_water = True
  .type = bool
  .help = Include water interactions (-waters/-nowaters in probe)

excluded_bond_chain_length = 4
  .type = int
  .help = Exclude chain of atoms bonded to source for this many hops (-4H, -3, -2 , -1 in probe)

drop_non_selected_atoms = False
  .type = bool
  .help = Drop non selected atoms (-drop in probe)

use_polar_hydrogens = True
  .type = bool
  .help = Use polar hydrogens (-usepolarh in probe)

polar_hydrogen_radius = 1.05
  .type = float
  .help = Radius to set for polar hydrogens (1.05 in probe)

minimum_polar_hydrogen_occupancy = 0.25
  .type = float
  .help = Minimum occupancy for polar hydrogens (0.25 in probe)

maximum_polar_hydrogen_b = 80.0
  .type = float
  .help = Minimum b-factor for polar hydrogens (80.0 in probe)

do_water_water = False
  .type = bool
  .help = Include water-to-water interactions (-wat2wat in probe)

do_water = True
  .type = bool
  .help = Include water-to-non-water interactions (-wat2wat, -waters, -nowaters in probe)

do_het = True
  .type = bool
  .help = Include het-atom interactions (-hets, -nohets in probe)

atom_radius_scale = 1.0
  .type = float
  .help = Atom radius = (r*atom_radius_scale)+atom_radius_offset (-scalevds, -vswscale in probe)

atom_radius_offset = 0.0
  .type = float
  .help = Atom radius = (r*atom_radius_scale)+atom_radius_offset (-addvdw in probe)

probe
  .style = menu_item auto_align
{
  radius = 0.25
    .type = float
    .help = Probe radius (half distance between touched atoms) (-radius in probe)

  density = 16.0
    .type = float
    .help = Probe dots per square angstrom on atom surface (-density in probe)

  worse_clash_cutoff = 0.5
    .type = float
    .help = Cutoff for worse clashes, a positive (-divworse in probe)

  clash_cutoff = 0.4
    .type = float
    .help = Cutoff for the clashes, a positive number (-divlow in probe)

  contact_cutoff = 0.25
    .type = float
    .help = Cutoff for the contact (-divhigh in probe)

  uncharged_hydrogen_cutoff = 0.6
    .type = float
    .help = Cutoff for uncharged hydrogen overlap (-hbregular in probe)

  charged_hydrogen_cutoff = 0.8
    .type = float
    .help = Cutoff for charged hydrogen overlap (-hbcharged in probe)

  bump_weight = 10.0
    .type = float
    .help = Weight applied to bump score (-bumpweight in probe)

  hydrogen_bond_weight = 4.0
    .type = float
    .help = Weight applied to hydrogen bond score (-hbweight in probe)

  gap_weight = 0.25
    .type = float
    .help = Weight applied to gap score (-gapweight in probe)

  separate_weak_hydrogen_bonds = False
    .type = bool
    .help = Separately account for weak hydrogen bonds (-LweakHbonds in probe)

  implicit_hydrogens = False
    .type = bool
    .help = Use implicit hydrogens, no water proxies (-implicit in probe)
}

output
  .style = menu_item auto_align
{
  file_name_prefix = None
    .type = path
    .short_caption = Prefix for file name
    .help = Prefix for file name
    .input_size = 400

  count_dots = False
    .type = bool
    .help = Count dots rather than listing all contacts (-countdots in probe)

  one_line = False
    .type = bool
    .help = Output one line :contacts:by:severity:type: (-oneline in probe)

  hydrogen_bond_output = True
    .type = bool
    .help = Output hydrogen-bond contacts (-nohbout in probe)

  record_added_hydrogens = False
    .type = bool
    .help = Output hydrogen-bond contacts (-dumph2o in probe)

  clash_output = True
    .type = bool
    .help = Output clash contacts (-noclashout in probe)

  vdw_output = True
    .type = bool
    .help = Output van der Waals contacts (-novdwout in probe)

  separate_worse_clashes = False
    .type = bool
    .help = Separately report worse clashes (-sepworse in probe)

  set_group_name = dots
    .type = str
    .help = Specify the group name (-name in probe)

  add_group_name_master_line = False
    .type = bool
    .help = Add a master=name line on lists (-dotmaster in probe)

  add_kinemage_keyword = False
    .type = bool
    .help = Add kinemage 1 to beginning of kin file (-kinemage in probe)

  add_lens_keyword = False
    .type = bool
    .help = Add lens keywoard to kin file (-lens in probe)

  add_group_statement = True
    .type = bool
    .help = Add lens keywoard to kin file (-nogroup in probe)

  color_by_dna_base = False
    .type = bool
    .help = Color by DNA base (-basecolor, -colorbase in probe)

  spike_scale_factor = 0.5
    .type = float
    .help = Fraction of overlap assigned to each atom (-spike in probe)

  group_label = ""
    .type = str
    .help = Label for the surface-dots group (-name, -scsurface, -exposed, -asurface, -access in probe)
}
'''

program_citations = phil.parse('''
citation {
  authors = Word, et. al.
  journal = J. Mol. Biol.
  volume = 285
  pages = 1711-1733
  year = 1999
  external = True
}
''')

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
Compute the MolProbity Probe score for a file, or a subset of the file.

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
Output:
  Kinemage file describing the score and other information, depending on the parameters.
Note:
  The source_selection phil parameter must always be filled in, and some
  approaches require the target_selection parameter as well.  Setting the
  target_selection to "=" will re-use the source for the target.  In all
  other cases, the string passed in will be used as a CCTBX selection on
  the model to select a subset of its atoms.
'''
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  citations = program_citations
  epilog = '''
  For additional information and help, see http://kinemage.biochem.duke.edu/software/probe
  and http://molprobity.biochem.duke.edu
  '''

# ------------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)
    if self.params.source_selection is None:
      raise Sorry("Must specify a source parameter for approach "+self.params.approach)
    if self.params.approach in ['once','both'] and self.params.target_selection is None:
      raise Sorry("Must specify a target parameter for approach "+self.params.approach)
    if self.params.output.file_name_prefix is None:
      raise Sorry("Supply the prefix for an output file name using output.file_name_prefix=")
    aScale = self.params.atom_radius_scale
    if aScale < 0.0001 or aScale > 1000:
      raise Sorry("Invalid atom_radius_scale value: {}".format(aScale))
    ao = self.params.atom_radius_offset
    if ao < -10 or ao > 1000:
      raise Sorry("Invalid atom_radius_offset value: {}".format(ao))

    # Ensure consistency among parameters
    if self.params.probe.contact_cutoff < self.params.probe.radius:
      self.params.probe.contact_cutoff = self.params.probe.radius

# ------------------------------------------------------------------------------

  def _scaled_atom_radius(self, a):
    '''
      Find the scaled and offset radius for the specified atom.  This will be called on each
      atom after their extra information has been loaded to determine the scaled and offset
      value to use for the remainder of the program.
      @todo How to handle asking for the implicit vs. explicit radius, which will presumably
      be done not here but in the helper function that gets extra atom info.
      :param a: Atom whose radius is to be scaled
      :return: Scaled and offset radius of the atom.
    '''
    rad = self._extraAtomInfo.getMappingFor(a).vdwRadius
    if rad <= 0:
      alt = a.parent().altloc
      resName = a.parent().resname.strip().upper()
      resID = str(a.parent().parent().resseq_as_int())
      chainID = a.parent().parent().parent().id
      myFullName = "chain "+str(chainID)+" "+resName+" "+resID+" "+a.name+" "+alt
      raise Sorry("Invalid radius for atom look-up: "+myFullName+"; rad = "+str(rad))
    return self.params.atom_radius_offset + (rad * self.params.atom_radius_scale)

# ------------------------------------------------------------------------------

  def _atom_class_for(self, a):
    '''
      Assign the atom class for a specified atom.
      :param a: Atom whose class is to be specified
      :return: If our parameters have been set to color and sort by NA base,
      then it returns the appropriate base name.  Otherwise, it returns the
      element of the atom.
    '''
    if not self.params.output.color_by_dna_base:
      return a.element
    else:
      resName = a.parent().resname
      cl = common_residue_names_get_class(name = resName)
      if cl == "common_rna_dna" or cl == "modified_rna_dna":
        cleanName = resName.upper().strip()
        if cleanName in ['U','URA','UTP','UDP','UMP','UR',
                         'T','THY','TTP','TDP','TMP','5MU','DT','TR']:
          return 't/u'
        elif cleanName in ['A','ADE','ATP','ADP','AMP','1MA','RIA','T6A','DA','AR']:
          return 'a'
        elif cleanName in ['C','CYT','CTP','CDP','CMP','5MC','OMC','DC','CR']:
          return 'c'
        elif cleanName in ['G','GUA','GTP','GDP','GMP','GSP','1MG','2MG','M2G','7MG','OMG','DG','GR']:
          return 'g'
        return 'other na'
      else:
        return "nonbase"

# ------------------------------------------------------------------------------

  def _color_for_atom_class(self, c):
    '''
      Report the color associated with an atom class.
      Based on atomprops.h:INIT_ATOM_TABLE from original probe.
      :param c: Class of the atom.
      :return: Kinemage name of the color associated with the class.
    '''

    # Make sure the atom class is one that we know about
    if not c in self._allAtomClasses:
      return 'magenta'

    # Check to see if this atom belongs to one of the special colors.
    if c in ['C','Ag','other na']:
      return 'white'
    elif c in ['N','He','t/u']:
      return 'sky'
    elif c in ['O']:
      return 'red'
    elif c in ['P','Ne','a']:
      return 'pink'
    elif c in ['S','c']:
      return 'yellow'
    elif c in ['Se','F','Cl']:
      return 'green'
    elif c in ['Br','I']:
      return 'brown'
    elif c in ['Co']:
      return 'blue'
    elif c in ['Cu','Ar']:
      return 'orange'
    elif c in ['Au']:
      return 'gold'
    elif c in ['Kr']:
      return 'greentint'
    elif c in ['Xe']:
      return 'magenta'
    elif c in ['Rn']:
      return 'pinktint'
    elif c in ['g']:
      return 'sea'

    # Most atom types, the default.
    return 'grey'

# ------------------------------------------------------------------------------

  class DotInfo:
    # Dot class storing information about an individual dot.
    def __init__(self, src, target, loc, spike, overlapType, gap, ptmaster, angle):
      self.src = src                  # Source atom for the interaction
      self.target = target            # Target atom for the interactions
      self.loc = loc                  # Location of the dot start
      self.spike = spike              # Location of the dot end
      self.overlapType = overlapType  # Type of overlap the interaction represents
      self.gap = gap                  # Gap between the atoms
      self.ptmaster = ptmaster        # Main/side chain interaction type
      self.angle = angle              # Angle associated with the bump

# ------------------------------------------------------------------------------

  def _save_dot(self, src, target, atomClass, loc, spike, overlapType, gap, ptmaster, angle):
    '''
      Generate and store a DotInfo entry with the specified parameters.  It will be stored
      into the self._results data structure.
      :param src: Source atom for the dot.
      :param target: Target atom for the dot, if any.
      :param atomClass: Atom class of this dot, indicates where to store.
      :param loc: Location of the dot start.
      :param spike: Location of the dot end.
      :param overlapType: Type of overlap for the dot.
      :param gap: Gap spacing for the dot.
      :param ptmaster: ptmaster entry for the dot.
      :param angle: angle for the dot.
      :return: As a side effect, this will add a new entry into one of the lists in the
      self._results data structure.
    '''
    self._results[atomClass][self._dotScorer.interaction_type(overlapType,gap,
        self.params.probe.separate_weak_hydrogen_bonds)].append(
      self.DotInfo(src, target, loc, spike, overlapType, gap, ptmaster, angle)
    )

# ------------------------------------------------------------------------------

  def _generate_surface_dots_for(self, src, bonded):
    '''
      Find all surface dots for the specified atom.
      This does not include locations where the probe is interacting with
      a bonded atom, so it is a subset of the skin dots where only the
      dots themselves are outside of the bonded atoms.
      :param src: Atom whose surface dots are to be found.
      :param bonded: Atoms that are bonded to src by one or more hops.
      :return: Side effect: Add dots on the surface of the atom to the
              self._results data structure by atomclass and dot type.
    '''

    # Generate no dots for ignored atoms.
    if self._atomClasses[src] == 'ignore':
      return

    # Check all of the dots for the atom and see if they should be
    # added to the list.
    srcInWater = self._inWater[src]
    r = self._extraAtomInfo.getMappingFor(src).vdwRadius
    pr = self.params.probe.radius
    srcDots = self._dots[src]
    for dotvect in srcDots:
      # Dot on the surface of the atom, at its radius; both dotloc and spikeloc from original code.
      # This is where the probe touches the surface.
      dotloc = Helpers.rvec3(src.xyz) + Helpers.rvec3(dotvect)
      # Dot that is one probe radius past the surface of the atom, exploring for contact with bonded
      # atoms.  This is the location of the center of the probe.
      exploc = Helpers.rvec3(src.xyz) + Helpers.rvec3(dotvect).normalize() * (r + pr)

      # If the exploring dot is within a probe radius + vdW radius of a bonded atom,
      # we don't add a dot.
      okay = True
      for b in bonded:
        bInWater = self._inWater[b]
        # If we should ignore the bonded element, we don't check it.
        if self._atomClasses[b] == 'ignore':
          continue
        # If we're ignoring water-water interactions and both src and
        # bonded are in a water, we should ignore this as well (unless
        # both are hydrogens from the same water, in which case we
        # continue on to check.)
        elif (not self.params.do_water_water and srcInWater and bInWater
              and src.parent() != b.parent() ):
          continue
        # If we're ignoring non-water to water interactions, skip check
        # if b is a water.
        elif not self.params.do_water and bInWater:
          continue
        # If we're ignoring hetatom interactions, skip check
        # if b is in a non-standard residue.
        elif not self.params.do_het and self._inHet[b]:
          continue

        # The bonded neighbor is one that we should check interaction with, see if
        # we're in range.  If so, mark this dot as not okay because it is inside a
        # bonded atom.
        if ( (Helpers.rvec3(b.xyz) - exploc).length() <=
            pr + self._extraAtomInfo.getMappingFor(b).vdwRadius ):
          okay = False

      # If this dot is okay, add it to the internal data structure based on its
      # atom class and overlap type.
      if okay:
        self._save_dot(src, None, self._atomClasses[src], dotloc, dotloc,
                       probeExt.OverlapType.NoOverlap, 0.0, ' ', 0.0)

# ------------------------------------------------------------------------------

  def _count_skin_dots_for(self, src, bonded):
    '''
      Count all skin dots for the specified atom.
      :param src: Atom whose surface dots are to be found.
      :param bonded: Atoms that are bonded to src by one or more hops.
      :return: Side effect: Add dots on the surface of the atom to the
              self._results data structure by atomclass and dot type.
    '''

    # No dots yet...
    ret = 0

    # Generate no dots for ignored atoms or for phantom Hydrogens.
    if self._atomClasses[src] == 'ignore' or self._extraAtomInfo.getMappingFor(src).isDummyHydrogen:
      return 0

    # Check all of the dots for the atom and see if they should be
    # added to the list.
    srcDots = self._dots[src]
    for dotvect in srcDots:
      # Dot on the surface of the atom, at its radius.
      # This is where the probe touches the surface.
      dotloc = Helpers.rvec3(src.xyz) + Helpers.rvec3(dotvect)

      # If the exploring dot is within a probe radius + vdW radius of a bonded atom,
      # we don't add a dot.
      okay = True
      for b in bonded:
        # If we should ignore the bonded element, we don't check it.
        if self._atomClasses[b] == 'ignore' or self._extraAtomInfo.getMappingFor(b).isDummyHydrogen:
          continue

        # The bonded neighbor is one that we should check interaction with, see if
        # we're in range.  If so, mark this dot as not okay because it is inside a
        # bonded atom.
        if ( (Helpers.rvec3(b.xyz) - dotloc).length() <=
             self._extraAtomInfo.getMappingFor(b).vdwRadius ):
          okay = False

      # If this dot is okay, add it to the internal data structure based on its
      # atom class and overlap type.
      if okay:
        ret += 1

    # Return the number of dots not inside a bonded atom.
    return ret

# ------------------------------------------------------------------------------

  def _count_skin_dots(self, atoms, bondedNeighborLists):
    '''
      Count all skin dots for the atoms passed in.
      :param atoms: Atoms to check.
      :param bondedNeighborLists: Neighbor list including these atoms.
      This is used to normalize output scores.
      :return: Number of skin dots on any of the atoms in the source selection.
    '''

    ret = 0

    maxVDWRadius = AtomTypes.AtomTypes().MaximumVDWRadius()
    for src in atoms:
      # Find nearby atoms that might come into contact.  This greatly speeds up the
      # search for touching atoms.
      maxRadius = (self._extraAtomInfo.getMappingFor(src).vdwRadius + maxVDWRadius +
        2 * self.params.probe.radius)
      nearby = self._spatialQuery.neighbors(src.xyz, 0.001, maxRadius)

      # Select those that are actually within the contact distance based on their
      # particular radius.
      atomList = []
      for n in nearby:
        d = (Helpers.rvec3(n.xyz) - Helpers.rvec3(src.xyz)).length()
        if (d <= self._extraAtomInfo.getMappingFor(n).vdwRadius +
            self._extraAtomInfo.getMappingFor(src).vdwRadius + 2*self.params.probe.radius):
          atomList.append(n)

      # Find the atoms that are bonded to the source atom within the specified hop
      # count.  Limit the length of the chain to 3 if neither the source nor the final
      # atom is a Hydrogen.
      neighbors = Helpers.getAtomsWithinNBonds(src, bondedNeighborLists,
        self.params.excluded_bond_chain_length, 3)

      # Count the skin dots for this atom.
      ret += self._count_skin_dots_for(src, neighbors)

    # Return the total count
    return ret

# ------------------------------------------------------------------------------

  def run(self):
    # String that will be output to the specified file.
    outString = ''

    if self.params.output.add_kinemage_keyword and not self.params.output.count_dots:
      outString += '@kinemage 1\n'

    make_sub_header('Interpret Model', out=self.logger)

    # Get our model.
    self.model = self.data_manager.get_model()

    # Fix up bogus unit cell when it occurs by checking crystal symmetry.
    cs = self.model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      self.model = shift_and_box_model(model = self.model)

    ################################################################################
    # Get the list of all atoms in the model
    atoms = self.model.get_atoms()

    ################################################################################
    # Get the bonding information we'll need to exclude our bonded neighbors.
    make_sub_header('Compute neighbor lists', out=self.logger)
    try:
      p = mmtbx.model.manager.get_default_pdb_interpretation_params()
      p.pdb_interpretation.use_neutron_distances = self.params.use_neutron_distances
      self.model.set_pdb_interpretation_params(params = p)
      self.model.process_input_model(make_restraints=True) # make restraints
      geometry = self.model.get_restraints_manager().geometry
      sites_cart = self.model.get_sites_cart() # cartesian coordinates
      bondProxies, asu = \
          geometry.get_all_bond_proxies(sites_cart = sites_cart)
    except Exception as e:
      raise Sorry("Could not get bonding information for input file: " + str(e))
    allBondedNeighborLists = Helpers.getBondedNeighborLists(atoms, bondProxies)

    ################################################################################
    # Get the extra atom information needed to score all of the atoms in the model.
    # @todo Handle getting implicit radii here or in the _scaled_atom_radius() function
    make_sub_header('Compute extra atom information', out=self.logger)
    ret = Helpers.getExtraAtomInfo(self.model)
    self._extraAtomInfo = ret.extraAtomInfo
    if len(ret.warnings) > 0:
      print('Warnings returned by getExtraAtomInfo():\n'+ret.warnings, file=self.logger)

    # Scale and offset the radius values for all atoms based on our command-line arguments.
    for a in atoms:
      ei = self._extraAtomInfo.getMappingFor(a)
      ei.vdwRadius = self._scaled_atom_radius(a)
      self._extraAtomInfo.setMappingFor(a, ei)

    ################################################################################
    # Get the extra atom information needed to sort all of the atoms in the model
    # into proper classes for reporting.  These classes may be atom names, when we're
    # sorting by atoms and it can be nucleic acid base names when we're sorting by that.
    # Comes from newAtom() and dotType() functions in probe.c.
    # Rather than a table indexed by type, we directly write the result.
    # Handle all atoms, not only selected atoms.
    self._atomClasses = {}
    for a in atoms:
      if a.element != 'H':
        # All elements except hydrogen use their own names.
        self._atomClasses[a] = self._atom_class_for(a)
      else:
        # For hydrogen, assign based on what it is bonded to.
        self._atomClasses[a] = self._atom_class_for(allBondedNeighborLists[a][0])

    ################################################################################
    # Get the dot sets we will need for each atom.  This is the set of offsets from the
    # atom center where dots should be placed.  We use a cache to reduce the calculation
    # time by returning the same answer for atoms that have the same radius.
    dotCache = probeExt.DotSphereCache(self.params.probe.density)
    self._dots = {}
    for a in atoms:
      self._dots[a] = dotCache.get_sphere(self._extraAtomInfo.getMappingFor(a).vdwRadius).dots()

    ################################################################################
    # Get the other characteristics we need to know about each atom to do our work.
    self._inWater = {}
    self._inHet = {}
    hetatm_sel = self.model.selection("hetatm")
    for a in atoms:
      self._inWater[a] = common_residue_names_get_class(name=a.parent().resname) == "common_water"
      self._inHet[a] = hetatm_sel[a.i_seq]

    ################################################################################
    # Get the source selection (and target selection if there is one).  These will be
    # lists of atoms that are in each selection, a subset of the atoms in the model.
    source_sel = self.model.selection(self.params.source_selection)
    source_atoms = set()
    for a in atoms:
      if source_sel[a.i_seq]:
        source_atoms.add(a)

    target_atoms = set()
    if self.params.target_selection is not None:
      # If the target selection is "=", that means that it should be the same as the source selection.
      if self.params.target_selection == "=":
        target_atoms = source_atoms
      else:
        target_sel = self.model.selection(self.params.target_selection)
        for a in atoms:
          if target_sel[a.i_seq]:
            target_atoms.add(a)

    ################################################################################
    # Find a list of all of the selected atoms with no duplicates
    # Get the bonded neighbor lists for the atoms that are in this selection.
    all_selected_atoms = set()
    for a in source_atoms:
      all_selected_atoms.add(a)
    for a in target_atoms:
      all_selected_atoms.add(a)
    bondedNeighborLists = Helpers.getBondedNeighborLists(all_selected_atoms, bondProxies)

    ################################################################################
    # Build a spatial-query structure that tells which atoms are nearby.
    # Include all atoms in the structure, not just the ones that have been selected.
    make_sub_header('Make spatial-queary accelerator', out=self.logger)
    self._spatialQuery = probeExt.SpatialQuery(atoms)

    ################################################################################
    # Construct a DotScorer object.
    make_sub_header('Make dot scorer', out=self.logger)
    self._dotScorer = probeExt.DotScorer(self._extraAtomInfo, self.params.probe.gap_weight,
      self.params.probe.bump_weight, self.params.probe.hydrogen_bond_weight,
      self.params.probe.uncharged_hydrogen_cutoff, self.params.probe.charged_hydrogen_cutoff,
      self.params.probe.clash_cutoff, self.params.probe.worse_clash_cutoff,
      self.params.probe.contact_cutoff)

    ################################################################################
    # If we're not doing implicit hydrogens, add Phantom hydrogens to waters and mark
    # the water oxygens as not being donors in atoms that are in the source or target selection.
    # Also clear the donor status of all N, O, S atoms because we have explicit hydrogen donors.
    if not self.params.probe.implicit_hydrogens:
      if self.params.output.record_added_hydrogens:
        outString += '@vectorlist {water H?} color= gray\n'

      # Check all atoms
      for a in all_selected_atoms:

        # If we are a hydrogen that is bonded to a nitrogen, oxygen, or sulfur then we're a donor.
        if a.element == 'H':
          # If we are in a water, make sure our occupancy and temperature (b) factor are acceptable.
          # If they are not, set the class for the atom to 'ignore'.
          if self._inWater[a] and (a.occ < self.params.minimum_polar_hydrogen_occupancy or
              a.b > self.params.maximum_polar_hydrogen_b):
            self._atomClasses[a] = 'ignore'
          else:
            for n in bondedNeighborLists[a]:
              if n.element in ['N','O','S']:
                # Copy the value, set the new values, then copy the new one back in.
                # We are a donor and may have our radius adjusted
                ei = self._extraAtomInfo.getMappingFor(a)
                ei.isDonor = True
                if self.params.use_polar_hydrogens:
                  ei.vdwRadius = self.params.atom_radius_offset + (
                    self.params.polar_hydrogen_radius * self.params.atom_radius_scale)
                self._extraAtomInfo.setMappingFor(a, ei)

                # Set our neigbor to not be a donor, since we are the donor
                ei = self._extraAtomInfo.getMappingFor(n)
                ei.isDonor = False
                self._extraAtomInfo.setMappingFor(n, ei)

        # If we are the Oxygen in a water, then add phantom hydrogens to nearby acceptors
        elif self._inWater[a] and a.element == 'O':
          # We're an acceptor and not a donor.
          ei = self._extraAtomInfo.getMappingFor(a)
          ei.isDonor = False
          ei.isAcceptor = True
          self._extraAtomInfo.setMappingFor(a, ei)

          # If we don't yet have Hydrogens attached, add phantom hydrogen(s)
          if len(bondedNeighborLists[a]) == 0:
            newPhantoms = Helpers.getPhantomHydrogensFor(a, self._spatialQuery, self._extraAtomInfo, 0.0, True)
            for p in newPhantoms:
              # Set all of the information other than the name and element and xyz of the atom based
              # on our parent; then overwrite the name and element and xyz.
              newAtom = pdb.hierarchy.atom(a.parent(),a)
              newAtom.name = " H?"
              newAtom.element = p.element
              newAtom.xyz = p.xyz

              # Add the atom to the spatial-query data structure
              self._spatialQuery.add(newAtom)

              # Set the extra atom information for this atom
              rad = self.params.atom_radius_offset + (
                    self.params.polar_hydrogen_radius * self.params.atom_radius_scale)
              ei = probeExt.ExtraAtomInfo(rad, False, True, True)
              self._extraAtomInfo.setMappingFor(newAtom, ei)

              # Set the atomClass and other data based on the parent Oxygen.
              self._atomClasses[newAtom] = self._atom_class_for(a)
              self._inWater[newAtom] = self._inWater[a]
              self._inHet[newAtom] = self._inHet[a]

              # Mark the new atom as being bonded to the parent atom, and the
              # parent atom as having to children.  Do this in both sets of bonded
              # neighbor lists.
              bondedNeighborLists[newAtom] = [a]
              bondedNeighborLists[a] = []
              allBondedNeighborLists[newAtom] = [a]
              allBondedNeighborLists[a] = []

              # Generate source dots for the new atom
              self._dots[newAtom] = dotCache.get_sphere(self._extraAtomInfo.getMappingFor(newAtom).vdwRadius).dots()

              # Add the new atom to any selections that the old atom was in.
              if a in source_atoms:
                source_atoms.add(newAtom)
              if a in target_atoms:
                target_atoms.add(newAtom)

              # Report on the creation if we've been asked to
              if self.params.output.record_added_hydrogens:

                resName = a.parent().resname.strip().upper()
                resID = str(a.parent().parent().resseq_as_int())
                chainID = a.parent().parent().parent().id
                iCode = a.parent().parent().icode
                alt = a.parent().altloc
                outString += '{{{:4s}{:1s}{:>3s}{:>2s}{:>4s}{}}}P {:8.3f}{:8.3f}{:8.3f}\n'.format(
                  a.name, alt, resName, chainID, resID, iCode,
                  a.xyz[0], a.xyz[1], a.xyz[2])

                resName = newAtom.parent().resname.strip().upper()
                resID = str(newAtom.parent().parent().resseq_as_int())
                chainID = newAtom.parent().parent().parent().id
                iCode = newAtom.parent().parent().icode
                alt = newAtom.parent().altloc
                outString += '{{{:4s}{:1s}{:>3s}{:>2s}{:>4s}{}}}L {:8.3f}{:8.3f}{:8.3f}\n'.format(
                  newAtom.name, alt, resName, chainID, resID, iCode,
                  newAtom.xyz[0], newAtom.xyz[1], newAtom.xyz[2])

        # Otherwise, if we're an N, O, or S then remove our donor status because
        # the hydrogens will be the donors
        elif a.element in ['N','O','S']:
          ei = self._extraAtomInfo.getMappingFor(a)
          ei.isDonor = False
          self._extraAtomInfo.setMappingFor(a, ei)

    ################################################################################
    # Re-fill all_selected_atoms
    all_selected_atoms = source_atoms.union(target_atoms)

    ################################################################################
    # List of all of the keys for atom classes, including all elements and all
    # nucleic acid types.  These are in the order that the original Probe reported
    # them.  Based on atomprops.h:INIT_ATOM_TABLE from original probe.
    self._allAtomClasses = ['ignore',
                         'H','C','N','O','P','S','As','Se','F','Cl','Br','I',
                         'Li','Na','Al','K','Mg','Ca','Mn','Fe','Co','Ni','Cu','Zn',
                         'Rb','Sr','Mo','Ag','Cd','In','Cs','Ba','Au','Hg','Tl','Pb',
                         'V','Cr','Te','Sm','Gd','Yb','W','Pt','U',
                         'He','Be','B','Ne','Se','Ar','Sc','Ti','Ga','Ge','Kr','Y','Zr',
                         'Sn','Sb','Xe','La','Ce','Fr','Ra','Th',
                         'Nb','Tc','Ru','Rh','Pd','Pr','Nd','Pm','Eu','Tb','Dy','Ho','Er',
                         'Tm','Lu','Hf','Ta','Re','Os','Ir','Bi','Po','At','Rn','Ac','Pa',
                         'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
                         'a','c','t/u','g','other na','nonbase']

    ################################################################################
    # Dictionary of dictionaries of lists structure holding lists of DotInfo class objects,
    # indexed by atom class and then by interaction type.  Fill in empty lists for all of
    # the possible classes and types.
    self._interactionTypes = [
        probeExt.InteractionType.WideContact,
        probeExt.InteractionType.CloseContact,
        probeExt.InteractionType.WeakHydrogenBond,
        probeExt.InteractionType.SmallOverlap,
        probeExt.InteractionType.Bump,
        probeExt.InteractionType.BadBump,
        probeExt.InteractionType.HydrogenBond
      ]
    self._results = {}
    for c in self._allAtomClasses:
      interactionTypeDicts = {}
      for i in self._interactionTypes:
        interactionTypeDicts[i] = []
      self._results[c] = interactionTypeDicts

    ################################################################################
    # Do the calculations; which one depends on the approach and other phil parameters.
    # Append the information to the string that will be written to file.
    if self.params.approach == 'count_atoms':
      make_sub_header('Counting atoms', out=self.logger)
      # Report the number of atoms in the source selection
      outString += 'atoms selected: '+str(len(source_atoms))+'\n'

    elif self.params.approach == 'surface':
      make_sub_header('Making surface dots', out=self.logger)

      # Produce dots on the surfaces of the selected atoms.
      maxVDWRadius = AtomTypes.AtomTypes().MaximumVDWRadius()
      for src in source_atoms:
        # Find nearby atoms that might come into contact.  This greatly speeds up the
        # search for touching atoms.
        maxRadius = (self._extraAtomInfo.getMappingFor(src).vdwRadius + maxVDWRadius +
          2 * self.params.probe.radius)
        nearby = self._spatialQuery.neighbors(src.xyz, 0.001, maxRadius)

        # Select those that are actually within the contact distance based on their
        # particular radius.
        atomList = []
        for n in nearby:
          d = (Helpers.rvec3(n.xyz) - Helpers.rvec3(src.xyz)).length()
          if (d <= self._extraAtomInfo.getMappingFor(n).vdwRadius +
              self._extraAtomInfo.getMappingFor(src).vdwRadius + 2*self.params.probe.radius):
            atomList.append(n)

        # Find the atoms that are bonded directly to the source atom.
        neighbors = Helpers.getAtomsWithinNBonds(src, allBondedNeighborLists, 1)

        # Find out what class of dot we should place for this atom.
        atomClass = self._atomClasses[src]

        # Generate all of the dots for this atom.
        self._generate_surface_dots_for(src, neighbors)

      # Count the dots if we've been asked to do so.
      if self.params.output.count_dots:
        numSkinDots = self._count_skin_dots(source_atoms, allBondedNeighborLists)
        name = "dots"
        if len(self.params.output.group_label) > 0:
          name = self.params.output.group_label
        outString += "selection: external\nname: {}\n".format(name)
        outString += "density: {:.1f} dots per A^2\nprobeRad: {:.3f} A\nVDWrad: (r * {:.3f}) + {:.3f} A\n".format(
          self.params.probe.density, self.params.probe.radius, self.params.atom_radius_scale,
          self.params.atom_radius_offset)

        nsel = len(source_atoms)
        # @todo

      # Otherwise, produce the dots as output
      else:
        # @todo
        pass

    # @todo
    else:

      # @todo Add Phantom Hydrogens to the ones we consider in cases where we should.

      make_sub_header('Compute Probe Score', out=self.logger)

      # @todo
      # Find the radius of each atom in the structure and construct dot spheres for
      # them. Find the atoms that are bonded to them and add them to an excluded list.
      # Then compute the score for each of them and report the summed score over the
      # whole molecule the way that Reduce will.
      total = 0
      badBumpTotal = 0
      for a in atoms:
        sphere = dotCache.get_sphere(self._extraAtomInfo.getMappingFor(a).vdwRadius)

        # Excluded atoms that are bonded to me or to one of my neightbors.
        # It has the side effect of excluding myself if I have any neighbors.
        # Construct as a set to avoid duplicates.
        exclude = set()
        for n in allBondedNeighborLists[a]:
          exclude.add(n)
          for n2 in allBondedNeighborLists[n]:
            exclude.add(n2)
        exclude = list(exclude)

        dots = sphere.dots()
        res = self._dotScorer.score_dots(a, 1.0, self._spatialQuery, rad*3, self.params.probe.radius,
              exclude, sphere.dots(), sphere.density(), False)
        total += res.totalScore()
        if res.hasBadBump:
          badBumpTotal += 1
      outString += 'Summed probe score for molecule = {:.2f} with {} bad bumps'.format(total, badBumpTotal)

    base = self.params.output.file_name_prefix
    of = open("%s.kin"%base,"w")
    of.write(outString)
    of.close()

# ------------------------------------------------------------------------------

  #def get_results(self):
  #  return group_args(model = self.model)
