from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import null_out
from libtbx import group_args, phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry
import mmtbx
import mmtbx_probe_ext as probeExt
from mmtbx.probe import Helpers

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

approach = *self both once out count
  .type = choice
  .help = self (src -> src) both (src <=> targ) once (src -> targ) surface (VdW surface) count (count atoms)

quiet = False
  .type = bool
  .help = Quiet mode, reduced output

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

count_dots = False
  .type = bool
  .help = Count dots rather than listing all contacts (-countdots in probe)

raw_output = False
  .type = bool
  .help = Unformatted output (-unformated in probe)

one_line = False
  .type = bool
  .help = Output one line :contacts:by:severity:type: (-oneline in probe)

hydrogen_bond_output = True
  .type = bool
  .help = Output hydrogen-bond contacts (-nohbout in probe)

clash_output = True
  .type = bool
  .help = Output clash contacts (-noclashout in probe)

vdw_output = True
  .type = bool
  .help = Output van der Waals contacts (-novdwout in probe)

separate_worst_clashes = False
  .type = bool
  .help = Separately report the worst clashes, below -0.5 (-sepworse in probe)

set_group_name = dots
  .type = str
  .help = Specify the group name (-name in probe)

add_group_name_master_line = False
  .type = bool
  .help = Add a master=name line on lists (-dotmaster in probe)

kinemage
  .style = menu_item auto_align
{
  add_lens_keyword = False
    .type = bool
    .help = Add lens keywoard to kin file (-lens in probe)

  add_group_statement = True
    .type = bool
    .help = Add lens keywoard to kin file (-nogroup in probe)
}

probe
  .style = menu_item auto_align
{
  radius = 0.25
    .type = float
    .help = Probe radius (half distance between touched atoms) (-radius in probe)

  density = 16.0
    .type = float
    .help = Probe dots per square angstrom on atom surface (-density in probe)

  worst_clash_cutoff = -0.5
    .type = float
    .help = Cutoff for the worst clashes (-divworse in probe)

  clash_cutoff = -0.4
    .type = float
    .help = Cutoff for the clashes (-divlow in probe)

  contact_cutoff = 0.25
    .type = float
    .help = Cutoff for the contact (-divhigh in probe)
}

output
  .style = menu_item auto_align
{
  file_name_prefix = None
    .type = path
    .short_caption = Prefix for file name
    .help = Prefix for file name
    .input_size = 400
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
'''
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  citations = program_citations
  epilog = '''
  For additional information and help, see http://molprobity.biochem.duke.edu
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

# ------------------------------------------------------------------------------

  def run(self):
    # If the target selection is "=", that means that it should be the same as the source selection.
    if self.params.target_selection == "=":
      self.params.target_selection = self.params.source_selection

    make_sub_header('Interpret Model', out=self.logger)

    # Get our model.
    self.model = self.data_manager.get_model()

    # Fix up bogus unit cell when it occurs by checking crystal symmetry.
    cs = self.model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      self.model = shift_and_box_model(model = self.model)

    # Get the list of all atoms in the model
    atoms = self.model.get_atoms()

    # Get the source selection (and target selection if there is one)
    # @todo
    source_atoms = atoms

    if self.params.approach == 'count':
      ret = 'atoms selected: '+str(len(source_atoms))

    else:

      # Get the bonding information we'll need to exclude our bonded neighbors.
      try:
        p = mmtbx.model.manager.get_default_pdb_interpretation_params()
        self.model.set_pdb_interpretation_params(params = p)
        self.model.process_input_model(make_restraints=True) # make restraints
        geometry = self.model.get_restraints_manager().geometry
        sites_cart = self.model.get_sites_cart() # cartesian coordinates
        bond_proxies_simple, asu = \
            geometry.get_all_bond_proxies(sites_cart = sites_cart)
      except Exception as e:
        raise Sorry("Could not get bonding information for input file: " + str(e))
      bondedNeighbors = Helpers.getBondedNeighborLists(atoms, bond_proxies_simple)

      # Construct a SpatialQuery and fill in the atoms.  Ensure that we can make a
      # query within 1000 Angstroms of the origin.
      sq = probeExt.SpatialQuery(atoms)

      make_sub_header('Fill in chemical information', out=self.logger)
      ret = Helpers.getExtraAtomInfo(self.model)
      extra = ret.extraAtomInfo
      if len(ret.warnings) > 0:
        print('Warnings returned by getExtraAtomInfo():\n'+ret.warnings, file=self.logger)

      make_sub_header('Compute Probe Score', out=self.logger)
      # Construct a DotScorer object.
      ds = probeExt.DotScorer(extra)

      # Construct a dot-sphere cache
      cache = probeExt.DotSphereCache(self.params.probe.density)

      # @todo
      # Find the radius of each atom in the structure and construct dot spheres for
      # them. Find the atoms that are bonded to them and add them to an excluded list.
      # Then compute the score for each of them and report the summed score over the
      # whole molecule the way that Reduce will.
      total = 0
      badBumpTotal = 0
      for a in atoms:
        rad = extra.getMappingFor(a).vdwRadius
        if rad <= 0:
          alt = a.parent().altloc
          resName = a.parent().resname.strip().upper()
          resID = str(a.parent().parent().resseq_as_int())
          chainID = a.parent().parent().parent().id
          myFullName = "chain "+str(chainID)+" "+resName+" "+resID+" "+a.name+" "+alt
          raise Sorry("Invalid radius for atom look-up: "+myFullName+"; rad = "+str(rad))
        sphere = cache.get_sphere(rad)

        # Excluded atoms that are bonded to me or to one of my neightbors.
        # It has the side effect of excluding myself if I have any neighbors.
        # Construct as a set to avoid duplicates.
        exclude = set()
        for n in bondedNeighbors[a]:
          exclude.add(n)
          for n2 in bondedNeighbors[n]:
            exclude.add(n2)
        exclude = list(exclude)

        dots = sphere.dots()
        res = ds.score_dots(a, 1.0, sq, rad*3, self.params.probe.radius, exclude, sphere.dots(), sphere.density(), False)
        total += res.totalScore()
        if res.hasBadBump:
          badBumpTotal += 1
      ret = 'Summed probe score for molecule = {:.2f} with {} bad bumps'.format(total, badBumpTotal)

    base = self.params.output.file_name_prefix
    of = open("%s.kin"%base,"w")
    of.write(ret)
    of.close()

# ------------------------------------------------------------------------------

  #def get_results(self):
  #  return group_args(model = self.model)
