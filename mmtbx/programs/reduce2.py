##################################################################################
# Copyright(c) 2021, Richardson Lab at Duke
# Licensed under the Apache 2 license
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissionsand
# limitations under the License.

from __future__ import absolute_import, division, print_function
import sys
import os
import time
from datetime import datetime
from libtbx.program_template import ProgramTemplate
from libtbx import group_args, phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry
import mmtbx
from mmtbx.probe import Helpers
from iotbx import pdb
# @todo See if we can remove the shift and box once reduce_hydrogen is complete
from cctbx.maptbx.box import shift_and_box_model
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.reduce import Optimizers
from libtbx.development.timers import work_clock
from scitbx.array_family import flex
from iotbx.pdb import common_residue_names_get_class

version = "0.6.0"

master_phil_str = '''
approach = *add remove
  .type = choice
  .short_caption = Add or remove Hydrogens
  .help = Determines whether Reduce will add (and optimize) or remove Hydrogens from the model
n_terminal_charge = *residue_one first_in_chain no_charge
  .type = choice(multi=False)
  .short_caption = N terminal charge approach
  .help = Mode for placing H3 at terminal nitrogen.
use_neutron_distances = False
  .type = bool
  .short_caption = Use neutron distances
  .help = Use neutron distances (-nuclear in reduce)
preference_magnitude = 1.0
  .type = float
  .short_caption = Rotational-preference magnitude
  .help = Multiplier on the rotational-preference energy for rotatable Movers (-penalty in reduce)
alt_id = None
  .type = str
  .short_caption = Alternate to optimize
  .help = Alternate to optimize.  The default is to optimize all of them.
model_id = None
  .type = int
  .short_caption = Model ID to optimize
  .help = Model ID to optimize.  The default is to optimize all of them.
add_flip_movers = False
  .type = bool
  .short_caption = Add flip movers
  .help = Insert flip movers (-flip, -build, -noflip, -demandflipallnhqs in reduce)
non_flip_preference = 0.5
  .type = float
  .short_caption = Preference to not flip
  .help = For flip movers, only do the flip if the score in the flipped orientation is this much better.
skip_bond_fix_up = False
  .type = bool
  .short_caption = Skip fixup step for Movers
  .help = For debugging purposes, it can be useful to only do flips with no bond fix-up to compare scores.
profile = False
  .type = bool
  .short_caption = Profile the entire run
  .help = Profile the performance of the entire run

output
  .style = menu_item auto_align
{
  description_file_name = None
    .type = str
    .short_caption = Description output file name
    .help = Description output file name
  flipkin_directory = None
    .type = str
    .short_caption = Where to place the Flipkin Kinemages
    .help = Where to place the Flipkin Kinemages. If None, no Flipkin files are made.
  print_atom_info = False
    .type = bool
    .short_caption = Print extra atom info
    .help = Print extra atom info
}
''' + Helpers.probe_phil_parameters

program_citations = phil.parse('''
citation {
  authors = Word, et. al.
  journal = J. Mol. Biol.
  volume = 285
  pages = 1735-1747
  year = 1999
  external = True
}
''')

# ------------------------------------------------------------------------------

class _MoverLocation(object):
  # Holds information needed to identify a Mover within a model file.
  def __init__(self, moverType, modelId, altId, chain, resName, resId):
    self.moverType = moverType  # String indicating Mover type: AmideFlip or HisFlip
    self.modelId = modelId      # Integer indicating the modelId that the entry corresponds to
    self.altId = altId          # String indicating the altId that the entry corresponds to
    self.chain = chain          # String Chain that the residue is in
    self.resName = resName      # String Name of the residue
    self.resId = resId          # Integer ID of the residue

  def __str__(self):
    return "{} {} '{}' {} {} {}".format(self.moverType, self.modelId, self.altId,
      self.chain, self.resName, self.resId)
  def __repr__(self):
      return "_MoverLocation({})".format(str(self))


def _FindMoversInOutputString(s, moverTypes = ['SingleHydrogenRotator',
      'NH3Rotator', 'AromaticMethylRotator', 'AmideFlip', 'HisFlip']):
  '''Return a list of _MoverLocation items that include all of them found in the
  output string from an Optimizer.
  :param s: String returned from the getInfo() method on an optimizer.
  :param moverTypes: List of names for the movertypes to find.
  :return: List of Optimizers.FlipMoverState objects indicating all flip Movers listed inside
  any BEGIN...END REPORT block in the string, including whether they were marked as
  flipped.
  '''
  ret = []
  modelId = None
  altId = None
  inBlock = False
  for line in s.splitlines():
    words = line.split()
    if inBlock:
      if words[0:2] == ['END','REPORT']:
        inBlock = False
      elif words[0] in moverTypes:
        ret.append(_MoverLocation(words[0], modelId, altId, words[3], words[4], int(words[5])))
    else:
      if words[0:2] == ['BEGIN','REPORT:']:
        modelId = int(words[3])
        # Remove the single-quote and colon characters from the AltId, possibly leaving it empty
        trim = words[5].replace("'", "")
        trim = trim.replace(":", "")
        altId = trim
        inBlock = True
  return ret


def _FindFlipsInOutputString(s, moverType):
  '''Return a list of _FlipMoverState items that include all of them found in the
  output string from an Optimizer.
  :param s: String returned from the getInfo() method on an optimizer.
  :param moverType: Either AmideFlip or HisFlip, selects the flip type to report.
  :return: List of Optimizers.FlipMoverState objects indicating all flip Movers listed inside
  any BEGIN...END REPORT block in the string, including whether they were marked as
  flipped.
  '''
  ret = []
  modelId = None
  altId = None
  inBlock = False
  for line in s.splitlines():
    words = line.split()
    if inBlock:
      if words[0:2] == ['END','REPORT']:
        inBlock = False
      elif words[0] == moverType:
        ret.append(Optimizers.FlipMoverState(moverType, modelId, altId, words[3], words[4],
                                   int(words[5]), words[14] == 'Flipped') )
    else:
      if words[0:2] == ['BEGIN','REPORT:']:
        modelId = int(words[3])
        # Remove the single-quote and colon characters from the AltId, possibly leaving it empty
        trim = words[5].replace("'", "")
        trim = trim.replace(":", "")
        altId = trim
        inBlock = True
  return ret


def _IsMover(rg, movers):
  '''Is a residue group in the list of residues that have Movers?
  :param rg: residue group
  :param movers: List of Optimizers.FlipMoverState objects
  :return: True if the residue is in the list, False if not
  '''
  for m in movers:
    chain = rg.parent()
    modelId = chain.parent().id
    # We must offset the index of the model by 1 to get to the 1-based model ID
    if ( (modelId == m.modelId + 1 or modelId == '') and
         (chain.id == m.chain) and
         (rg.resseq_as_int() == m.resId)
       ):
      return True
  return False

def _AltFromFlipOutput(fo):
  '''Return a string that describes the alternate from a record in _FindFlipsInOutputString()
  output.  Reports ' ' for an empty alternate. Returns lower-case representation.
  '''
  if fo.altId in ["", ""]:
    return ' '
  return fo.altId.lower()


def _AddPointOrLineTo(a, tag, group):
  '''Return a string that describes the point at or line to the specified atom.
  This is used when building Kinemages.  Reports the alternate only if it is not empty.
  :param a: Atom to describe.
  :param tag: 'P' for point, 'L' for line.
  :param group: The dominant group name the point or line is part of.
  '''
  if a.parent().altloc in ['', ' ']:
    altTag = ''
  else:
    altTag = " '{}'".format(a.parent().altloc)
  return '{{{:.4s} {} {} {:3d} B{:.2f} {}}} {}{} {:.3f}, {:.3f}, {:.3f}'.format(
    a.name.strip().lower(),               # Atom name
    a.parent().resname.strip().lower(),   # Residue name
    a.parent().parent().parent().id,      # chain
    a.parent().parent().resseq_as_int(),  # Residue number
    a.b,                                  # B factor
    group,                                # Dominant group name
    tag,                                  # Tag (P or L)
    altTag,                               # Alternate, if any
    a.xyz[0],                             # location
    a.xyz[1],
    a.xyz[2]
  )

def _DescribeMainchainResidue(r, group, prevC):
  '''Return a string that describes the mainchain for a specified residue.
  Add the point for the first mainchain atom in the previous residue
  (none for the first) and lines to the N, Ca, C, and O.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param prevC: Atom that is the mainchain C for the previous residue (or
  the N for this residue if we're the first residue in a chain).
  '''

  ret = ''
  # Find the atoms that we're going to use and insert their records into
  # the output. If we're missing an atom, skip its record and make sure that
  # the first one we use has a P and the others have L.
  try:
    aN = [a for a in r.atoms() if a.name.strip().upper() == 'N'][0]
    # If we're at the N terminus, we won't have a previous C, so just re-use our N
    if prevC is None:
      prevC = aN
    ret += _AddPointOrLineTo(prevC, 'P', group) + ' ' + _AddPointOrLineTo(aN, 'L', group) + '\n'
    tag = 'L'
  except Exception:
    tag = 'P'

  try:
    aCA = [a for a in r.atoms() if a.name.strip().upper() == 'CA'][0]
    ret += _AddPointOrLineTo(aCA, tag, group) + '\n'
    tag = 'L'
  except Exception:
    tag = 'P'

  try:
    aC = [a for a in r.atoms() if a.name.strip().upper() == 'C'][0]
    ret += _AddPointOrLineTo(aC, tag, group) + '\n'
    tag = 'L'
  except Exception:
    tag = 'P'

  try:
    aO = [a for a in r.atoms() if a.name.strip().upper() == 'O'][0]
    ret += _AddPointOrLineTo(aO, tag, group) + '\n'
    tag = 'L'
  except Exception:
    tag = 'P'

  # If we're the C terminus, we'll have an OXT atom.  In this case, add
  # a point at the C and a line to OXT.
  try:
    aOXT = [a for a in r.atoms() if a.name.strip().upper() == 'OXT'][0]
    ret += _AddPointOrLineTo(aC, 'P', group) + ' ' + _AddPointOrLineTo(aOXT, 'L', group) + '\n'
  except Exception:
    pass

  return ret


def _DescribeMainchainResidueHydrogens(r, group, bondedNeighborLists):
  '''Return a string that describes the mainchain hydrogens for a specified residue.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  # Find all of the Hydrogens in the residue
  Hs = [a for a in r.atoms() if a.element_is_hydrogen()]
  for h in Hs:
    try:
      n = bondedNeighborLists[h][0]
      # If the hydrogen is bonded to a mainchain atom, add it
      if n.name.strip().upper() in ['N', 'CA', 'C', 'O']:
        ret += _AddPointOrLineTo(n, 'P', group) + ' ' + _AddPointOrLineTo(h, 'L', group) + '\n'
    except Exception:
      pass

  return ret


def _DescribeSidechainResidue(r, group, bondedNeighborLists):
  '''Return a string that describes the sidechain non-hydrogen portions for a specified residue.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  # Start with the CA atom and mark as handled all links that go back to the main chain
  # or to Hydrogens.  Add the CA to the list of atoms queued to be handled.
  described = []   # A list of sets of two atoms that we have already described bonds between
  queued = []
  try:
    # Get all of the neighbors of CA that are not N or C.  Queue them for testing.
    aCA = [a for a in r.atoms() if a.name.strip().upper() == 'CA'][0]
    queued.append(aCA)
    known = [a for a in bondedNeighborLists[aCA] if a.name.strip().upper() in ['N','C']]
    for a in known:
      described.append({aCA, a})
  except Exception as e:
    pass

  # Cycle through the list of queued atoms until we run out of them.
  # For each, look for a non-hydrogen neighbor that we've not yet described a bond
  # with.  If none are found, remove this entry from the list and cycle again.
  # If one is found, add a point and then a line for the first neighbor found, adding
  # that neighbor and all others to the queued list and continuing to chase to (line to
  # the first, adding it and all others to the queued list) until we find no more links
  # to atoms in the same residue.
  while len(queued) > 0:
    last = queued[0]
    links = [a for a in bondedNeighborLists[last]
              if (not {last, a} in described) and not a.element_is_hydrogen()]
    if len(links) == 0:
      # First entry on the list yielded no useful neightbors; remove it and check the next
      queued = queued[1:]
      continue
    ret += _AddPointOrLineTo(last, 'P', group) + ' '
    while len(links) != 0:
      # Put all but the first link into the list to be checked later.
      for a in links[1:]:
        queued.append(a)
      # Add the description for our first one and keep chasing this path
      curr = links[0]
      described.append({last,curr})
      ret += _AddPointOrLineTo(curr, 'L', group) + '\n'
      links = [a for a in bondedNeighborLists[curr]
                if (not {curr, a} in described) and not a.element_is_hydrogen()
                and curr.parent() == a.parent()
              ]
      last = curr
  return ret


def _DescribeSidechainResidueHydrogens(r, group, bondedNeighborLists):
  '''Return a string that describes the sidechain hydrogens for a specified residue.
  :param r: Residue to describe.
  :param group: The dominant group name the point or line is part of.
  :param bondedNeighborLists: A dictionary that contains an entry for each atom in the
  structure that the atom from the first parameter interacts with that lists all of the
  bonded atoms.  Can be obtained by calling mmtbx.probe.Helpers.getBondedNeighborLists().
  '''
  ret = ''

  # Find all of the Hydrogens in the residue
  Hs = [a for a in r.atoms() if a.element_is_hydrogen()]
  for h in Hs:
    try:
      n = bondedNeighborLists[h][0]
      # If the hydrogen is bonded to a mainchain atom, add it
      if not n.name.strip().upper() in ['N', 'CA', 'C', 'O']:
        ret += _AddPointOrLineTo(n, 'P', group) + ' ' + _AddPointOrLineTo(h, 'L', group) + '\n'
    except Exception:
      pass

  return ret


def _AddFlipkinBase(states, views, fileName, fileBaseName, model, alts, bondedNeighborLists,
    moverList, inSideChain, inWater):
  '''Return a string that forms the basis for a Flipkin file without the optional positions
  for the specified movers.  This includes the views that will be used to look at them.
  :param states: Return value from _FindFlipsInOutputString() indicating behavior of
  each Mover.
  :param views: List of viewpoints, one per entry in states.
  :param fileName: Name of the optimized output file associated with this Flipkin.
  :param fileBaseName: The base name of the file, without path or extension.
  :param model: The model we're optimizing.
  :param alts: A list of alternates, empty if there are none. Sorted in increasing order.
  :param bondedNeighborLists: List of neighboring atoms bonded to each atom.
  :param moverList: List of Movers, which will be used to exclude residues.
  :param inSideChain: Dictionary looked up by atom telling whether it is in a side chain.
  :param inWater: Dictionary looked up by atom telling whether it is in water.
  '''
  ret = '@kinemage 1\n'
  ret += '@caption\nfrom file: {}\n'.format(fileName)

  # Compute the views for each Mover as the center of mass of all of the moving atoms and
  # record them, indicating which are flipped in Reduce.
  ret += ' views marked with * are for groups flipped by reduce\n'
  for i, s in enumerate(states):
    # See whether the state is flipped in Reduce and add a star if so
    star = ' '
    if s.flipped:
      star = '*'

    # Find out the type of the residue, used to determine the type of flip.
    type = '?'
    if s.resName == 'ASN':
      type = 'N'
    elif s.resName == 'GLN':
      type = 'Q'
    elif s.resName == 'HIS':
      type = 'H'

    if i > 0:
      indexString = str(i+1)
    else:
      indexString = ''
    ret += '@{}viewid {{{}{}{} {} {}}}\n'.format(indexString, star, type, s.resId, _AltFromFlipOutput(s), s.chain)
    ret += '@{}span 12\n'.format(indexString)
    ret += '@{}zslab 100\n'.format(indexString)
    ret += '@{}center{:9.3f}{:9.3f}{:9.3f}\n'.format(indexString, views[i][0], views[i][1], views[i][2])

  # Add the master descriptions
  ret += '@master {mainchain}\n'
  ret += '@master {sidechain}\n'
  ret += "@master {H's}\n"
  ret += '@master {hets}\n'
  ret += '@master {water}\n'

  # Add width and group description, which is the base name of the file without
  # its extension
  ret += '@onewidth\n'
  ret += '@group {{{}}} dominant\n'.format(fileBaseName)

  # Add the masters for the alternates if there are any. The sorted list always starts
  # with '' and then adds the others in increasing order.
  defaultAltSet = False
  for a in alts:
    # The first one is turned on and the others are turned off by default
    if defaultAltSet:
      state = 'off'
    else:
      defaultAltSet = True
      state = 'on'
    ret += "@pointmaster '{}' {{{}}} {}".format(a.lower(), 'alt'+a.lower(), state)

  # Add the mainchain (no hydrogens) for atoms all atoms (even Movers of the type
  # we're looking at right now). Add the point for the first mainchain
  # atom in the previous residue (None for the first).
  ret += '@subgroup {{mc {}}} dominant\n'.format(fileBaseName)
  ret += '@vectorlist {mc} color= white  master= {mainchain}\n'
  for c in model.chains():
    # @todo What to do about mainchains with alternate conformations?
    prevC = None
    for rg in c.residue_groups():
      ret += _DescribeMainchainResidue(rg, fileBaseName, prevC)
      try:
        prevC = [a for a in rg.atoms() if a.name.strip().upper() == 'C'][0]
      except Exception:
        pass

  # Add the Hydrogens on the mainchain
  ret += "@vectorlist {mc H} color= gray  master= {mainchain} master= {H's}\n"
  for c in model.chains():
    prevC = None
    for rg in c.residue_groups():
      ret += _DescribeMainchainResidueHydrogens(rg, fileBaseName, bondedNeighborLists)

  # Add the sidechain non-hydrogen atoms for residues that do not have Movers
  ret += '@subgroup {{sc {}}} dominant\n'.format(fileBaseName)
  ret += '@vectorlist {sc} color= cyan  master= {sidechain}\n'
  for c in model.chains():
    for rg in c.residue_groups():
      if not _IsMover(rg, moverList):
        ret += _DescribeSidechainResidue(rg, fileBaseName, bondedNeighborLists)

  # Add the Hydrogens on the sidechains for residues that do not have Movers
  ret += "@vectorlist {mc H} color= gray  master= {sidechain} master= {H's}\n"
  for c in model.chains():
    for rg in c.residue_groups():
      if not _IsMover(rg, moverList):
        ret += _DescribeSidechainResidueHydrogens(rg, fileBaseName, bondedNeighborLists)

  # Describe links between atoms in two different sidechains where neither of the
  # involved residues include Movers.  Don't repeat bonds that have already been
  # described.
  ret += '@vectorlist {SS} color= yellow  master= {sidechain}\n'
  described = []
  for a in model.get_atoms():
    for n in bondedNeighborLists[a]:
      if a.parent().parent() != n.parent().parent() and inSideChain[a] and inSideChain[n]:
        if {a,n} not in described:
          ret += _AddPointOrLineTo(a, 'P', fileBaseName) + ' ' + _AddPointOrLineTo(n, 'L', fileBaseName) + '\n'
          described.append({a,n})

  # Add het groups (ions and bonded structures)
  # @todo

  # Add waters
  ret += '@subgroup waters dominant\n'
  ret += '@balllist {water O} color= pink  radius= 0.15  master= {water}\n'
  for a in model.get_atoms():
    if inWater[a]:
      ret += _AddPointOrLineTo(a, 'P', fileBaseName)

  # @todo

  return ret

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
Reduce2 version {}
Add Hydrogens to a model and optimize their placement by adjusting movable groups and
flippable groups of atoms.

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
Output:
  PDB or mmCIF file with added hydrogens.  If output.filename is specified, then the
  type of file to write will be determined by its suffix (.pdb or .cif).
  If output.filename is not specified, the output file will be
  written into the current working directory with the same base name and type as the
  original file and with FH or H added to the base name (FH when flips are requested);
  1xs0.pdb would be written to ./1xsoH.pdb and 1xso.cif to ./1xsoH.cif by default.

NOTES:
  If multiple alternates are present in the file and a specific one is not specified on the
  command line, they will all be processed in reverse order, such that the lowest-named
  one (perhaps A) will be processed last.  The hydrogen addition and arrangements for
  residues that are not part of any alternate will be left in the configuration that is best
  for the final alternate tested.  This may leave other alternates in sub-optimal configurations.
  When a single alternate is selected using alt_id= on the command line, everything is
  optimized a single time and for that configuration.

  Note that the program also takes probe Phil arguments; run with --show_defaults to see
  all Phil arguments.

  Equivalent PHIL arguments for original Reduce command-line options:
    -quiet: No equivalent; metadata is never written to the model file, it is always
            written to the description file, and progress information is always written
            to standard output.
    -trim: approach=remove
    -build: approach=add add_flip_movers=True
    -flip: approach=add add_flip_movers=True
    -allalt: This is the default.
    -penalty200: preference_magnitude=200
    -nobuild9999: approach=add preference_magnitude=9999
    -noflip: approach=add add_flip_movers=True preference_magnitude=9999
    -onlya: alt_id=A
    -nuclear: use_neutron_distances=True
    -demandflipallnhqs: add_flip_movers=True
'''.format(version)
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']
  citations = program_citations
  epilog = '''
  For additional information and help, see http://kinemage.biochem.duke.edu/software/reduce
  and http://molprobity.biochem.duke.edu
  '''

# ------------------------------------------------------------------------------

  def validate(self):
    # Set the default output file name if one has not been given.
    if self.params.output.filename is None:
      inName = self.data_manager.get_default_model_name()
      suffix = os.path.splitext(os.path.basename(inName))[1]
      if self.params.add_flip_movers:
        pad = 'FH'
      else:
        pad = 'H'
      base = os.path.splitext(os.path.basename(inName))[0] + pad
      self.params.output.filename = base + suffix
      print('Writing model output to', self.params.output.filename)

    self.data_manager.has_models(raise_sorry=True)
    if self.params.output.description_file_name is None:
      raise Sorry("Must specify output.description_file_name")

    # Check the model ID to make sure they didn't set it to 0
    if self.params.model_id == 0:
      raise Sorry("Model ID must be >=1 if specified (None means all models)")

    # Turn on profiling if we've been asked to in the Phil parameters
    if self.params.profile:
      import cProfile
      self._pr = cProfile.Profile()
      self._pr.enable()

# ------------------------------------------------------------------------------

  def run(self):

    # String describing the run that will be output to the specified file.
    outString = 'reduce2 v.{}, run {}\n'.format(version, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for a in sys.argv:
      outString += ' {}'.format(a)
    outString += '\n'

    make_sub_header('Loading Model', out=self.logger)

    # Get our model.
    self.model = self.data_manager.get_model()

    # Stores the initial coordinates for all of the atoms, including added hydrogens,
    # for use in the Flipkins. Holds a list of tuples where the first is the atom
    # record and the second is the 3-vector position of the atom.
    initialAtomPositions = []

    # Fix up bogus unit cell when it occurs by checking crystal symmetry.
    cs = self.model.crystal_symmetry()
    if (cs is None) or (cs.unit_cell() is None):
      self.model = shift_and_box_model(model = self.model)

    if self.params.approach == 'add':
      # Add Hydrogens to the model
      make_sub_header('Adding Hydrogens', out=self.logger)
      startAdd = work_clock()
      reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
        model = self.model,
        use_neutron_distances=self.params.use_neutron_distances,
        n_terminal_charge=self.params.n_terminal_charge,
        exclude_water=True,
        stop_for_unknowns=True,
        keep_existing_H=False
      )
      reduce_add_h_obj.run()
      reduce_add_h_obj.show(None)
      missed_residues = set(reduce_add_h_obj.no_H_placed_mlq)
      if len(missed_residues) > 0:
        bad = ""
        for res in missed_residues:
          bad += " " + res
        raise Sorry("Restraints were not found for the following residues:"+bad)
      self.model = reduce_add_h_obj.get_model()
      doneAdd = work_clock()

      # Interpret the model after shifting and adding Hydrogens to it so that
      # all of the needed fields are filled in when we use them below.
      # @todo Remove this once place_hydrogens() does all the interpretation we need.
      make_sub_header('Interpreting Hydrogenated Model', out=self.logger)
      startInt = work_clock()
      self.model.get_hierarchy().sort_atoms_in_place()
      self.model.get_hierarchy().atoms().reset_serial()
      p = mmtbx.model.manager.get_default_pdb_interpretation_params()
      p.pdb_interpretation.allow_polymer_cross_special_position=True
      p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
      p.pdb_interpretation.use_neutron_distances = self.params.use_neutron_distances
      p.pdb_interpretation.proceed_with_excessive_length_bonds=True
      #p.pdb_interpretation.sort_atoms=True
      self.model.process(make_restraints=True, pdb_interpretation_params=p) # make restraints
      doneInt = work_clock()

      # Keep track of the initial atom positions after hydrogen placement but before
      # optimization so that we can put things back before each Flipkin run.
      for a in self.model.get_hierarchy().atoms():
        initialAtomPositions.append( (a, a.xyz) )

      make_sub_header('Optimizing', out=self.logger)
      startOpt = work_clock()
      opt = Optimizers.FastOptimizer(self.params.probe, self.params.add_flip_movers,
        self.model, probeRadius=0.25, altID=self.params.alt_id, modelIndex=self.params.model_id,
        preferenceMagnitude=self.params.preference_magnitude,
        nonFlipPreference=self.params.non_flip_preference,
        skipBondFixup=self.params.non_flip_preference,
        verbosity=3)
      doneOpt = work_clock()
      outString += opt.getInfo()
      outString += 'Time to Add Hydrogen = '+str(doneAdd-startAdd)+'\n'
      outString += 'Time to Interpret = '+str(doneInt-startInt)+'\n'
      outString += 'Time to Optimize = '+str(doneOpt-startOpt)+'\n'
      if self.params.output.print_atom_info:
        print('Atom information used during calculations:')
        print(opt.getAtomDump())

    else: # Removing Hydrogens from the model rather than adding them.
      make_sub_header('Removing Hydrogens', out=self.logger)
      sel = self.model.selection("element H")
      for a in self.model.get_atoms():
        if sel[a.i_seq]:
          a.parent().remove_atom(a)

    # Re-process the model because we have removed some atoms that were previously
    # bonded.  Don't make restraints during the reprocessing.
    # We had to do this to keep from crashing on a call to pair_proxies when generating
    # mmCIF files, so we always do it for safety.
    self.model.process(make_restraints=False, pdb_interpretation_params=p)

    make_sub_header('Writing output', out=self.logger)

    # Write the description output to the specified file.
    self.data_manager._write_text("description", outString,
      self.params.output.description_file_name)

    # Determine whether to write a PDB or CIF file and write the appropriate text output.
    suffix = os.path.splitext(self.params.output.filename)[1]
    if suffix.lower() == ".pdb":
      txt = self.model.model_as_pdb()
    else:
      txt = self.model.model_as_mmcif()
    self.data_manager._write_text("model", txt, self.params.output.filename)

    print('Wrote', self.params.output.filename,'and',
      self.params.output.description_file_name, file = self.logger)

    # If we've been asked to make Flipkins, then make each of them.
    if self.params.add_flip_movers and self.params.output.flipkin_directory is not None:

      # Find the base name of the two output files we will produce.
      inName = self.data_manager.get_default_model_name()
      suffix = os.path.splitext(os.path.basename(inName))[1]
      pad = 'FH'
      base = os.path.splitext(os.path.basename(inName))[0] + pad
      flipkinBase = self.params.output.flipkin_directory + "/" + base

      # Find the list of all Movers in the model, which will be used to segment
      # it into parts for the Flipkin.
      moverLocations = _FindMoversInOutputString(outString)
      print('XXX Movers:', moverLocations)

      # Find the list of all alternates in the model, ignoring empty ones.
      # Sort them in increasing alphabetical order.
      alts = Optimizers.AlternatesInModel(self.model)
      alts.discard('')
      alts.discard(' ')
      alts = sorted(list(alts))

      make_sub_header('Generating Amide Flipkin', out=self.logger)

      # Make list of Movers to lock in one flip orientation and then the other,
      # keeping track of which state they are in when Reduce was choosing.
      # We need a different list for the Amide Movers and the Histidine Movers
      # because we generate two different Flipkin files, one for each.
      amides = _FindFlipsInOutputString(outString, 'AmideFlip')

      # Find the viewpoint locations for each Mover we're going to
      # look at as the center of all atoms in the sidechain of the residue.
      views = []
      for a in amides:
        # Fill in information needed to construct the view.
        # As of 2/5/2023, the CCTBX selection returns no atoms on a file when the model
        # clause is used unless there is a MODEL statement in the file.  The get_number_of_models()
        # function returns 1 if there are 0 or 1 MODEL statements, so we check to see if there
        # are 2 or more (indicating the need to select) before adding the clause.
        # The model ID that the selection is looking for is 1-based, so we must add 1 to the
        # model index.
        if self.model.get_number_of_models() >= 2:
          modelClause = 'model {} and '.format(a.modelId + 1)
        else:
          modelClause = ''
        x = 0.0
        y = 0.0
        z = 0.0
        if a.altId in ["", " "]:
          selString = modelClause + "chain {} and resseq {} and sidechain".format(
                a.chain, a.resId)
        else:
          selString = modelClause + "chain {} and altid '{}' and resseq {} and sidechain".format(
                a.chain, a.altId, a.resId)
        sel = self.model.selection(selString)
        count = 0;
        for a in self.model.get_hierarchy().atoms():
          if sel[a.i_seq]:
            x += a.xyz[0]
            y += a.xyz[1]
            z += a.xyz[2]
            count += 1
        if count > 0:
          x /= count
          y /= count
          z /= count
        views.append( (x, y, z) )

      # Find out which atoms are bonded.
      p = mmtbx.model.manager.get_default_pdb_interpretation_params()
      p.pdb_interpretation.allow_polymer_cross_special_position=True
      p.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
      p.pdb_interpretation.proceed_with_excessive_length_bonds=True
      self.model.process(make_restraints=True,pdb_interpretation_params=p) # make restraints
      carts = flex.vec3_double()
      for a in self.model.get_atoms():
        carts.append(a.xyz)
      bondProxies = self.model.get_restraints_manager().geometry.get_all_bond_proxies(sites_cart = carts)[0]
      bondedNeighborLists = Helpers.getBondedNeighborLists(self.model.get_atoms(), bondProxies)

      ################################################################################
      # Get the other characteristics we need to know about each atom to do our work.
      inWater = {}
      inHet = {}
      inMainChain = {}
      inSideChain = {}
      hetatm_sel = self.model.selection("hetatm")
      mainchain_sel = self.model.selection("backbone")  # Will NOT include Hydrogen atoms on the main chain
      sidechain_sel = self.model.selection("sidechain") # Will include Hydrogen atoms on the side chain
      for a in self.model.get_atoms():
        inWater[a] = common_residue_names_get_class(name=a.parent().resname) == "common_water"
        inHet[a] = hetatm_sel[a.i_seq]
        if not a.element_is_hydrogen():
          inMainChain[a] = mainchain_sel[a.i_seq]
        else:
          # Check our bonded neighbor to see if it is on the mainchain if we are a Hydrogen
          if len(bondedNeighborLists[a]) != 1:
            raise Sorry("Found Hydrogen with number of neigbors other than 1: "+
                        str(len(bondedNeighborLists[a])))
          else:
            inMainChain[a] = mainchain_sel[bondedNeighborLists[a][0].i_seq]
        inSideChain[a] = sidechain_sel[a.i_seq]

      # Write the base information in the Flipkin, not including the moving atoms in
      # the Movers that will be placed, or atoms bonded to the moving atoms.
      flipkinText = _AddFlipkinBase(amides, views, self.params.output.filename, base, self.model,
        alts, bondedNeighborLists, moverLocations, inSideChain, inWater)

      # @todo Make two configurations, the one that Reduce picked and the one
      # that it did not.
      configurations = []

      # @todo Figure out how to run the optimization without fixup on the
      # original orientation of all flippers and again on the flipped orientation
      # for each and combine the info from both of them into the same Flipkin. This
      # is made complicated by the addition and deletion of atoms.
      for i, c in enumerate(configurations):
        # Put the atoms back to where they started before optimization.
        for p in initialAtomPositions:
          p[0].xyz = p[1]

        # Run optimization, locking the specified Amides into each configuration.
        # Don't do fixup.
        # @todo

        # Write the updates to the Flipkin for this configuration, showing the
        # atoms for the amide in the Reduce configuration (i=0) or the other
        # configuration (i=1).
        # @todo

      # Write the accumulated Flipkin string to the output file.
      with open(flipkinBase+"-flipnq.kin", "w") as f:
        f.write(flipkinText)

      make_sub_header('Generating Histidine Flipkin', out=self.logger)
      hists = _FindFlipsInOutputString(outString, 'HisFlip')
      print('XXX His:', hists)

      # @todo Enable subsetting the states of a Histidine so that it will try
      # all four protenations at a specific flip orientation, so that we can
      # ask it to do each set separately.
      # @todo Enabling setting the state (flipped or not) of each of the Movers
      # and passing a list of locked ones to the Optimizer in two passes for
      # each Flipkin file.
      # @todo Figure out how to add the Hydrogens back that have been removed
      # during Histidine placement.  The description will tell us which ones have
      # not been placed (so were removed).
      # @todo

    # Report profiling info if we've been asked to in the Phil parameters
    if self.params.profile:
      print('Profile results:')
      import pstats
      profile_params = {'sort_by': 'time', 'num_entries': 20}
      self._pr.disable()
      ps = pstats.Stats(self._pr).sort_stats(profile_params['sort_by'])
      ps.print_stats(profile_params['num_entries'])

# ------------------------------------------------------------------------------

  def get_results(self):
    return group_args(model = self.model)

# ------------------------------------------------------------------------------

  def Test(self):
    '''
      Run tests on the methods of the class.  Throw an assertion error if there is a problem with
      one of them and return normally if there is not a problem.
    '''

    #=====================================================================================
    # @todo Unit tests for other methods
