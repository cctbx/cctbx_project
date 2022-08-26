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

version = "0.3.0"

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
add_flip_movers = False
  .type = bool
  .short_caption = Add flip movers
  .help = Insert flip movers (-flip, -build, -noflip, -demandflipallnhqs in reduce)
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

class Program(ProgramTemplate):
  description = '''
Reduce2 version {}
Add Hydrogens to a model and optimize their placement by adjusting movable groups and
flippable groups of atoms.

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
Output:
  PDB or mmCIF file with added hydrogens.  If output.file_name is specified, then the
  type of file to write will be determined by its suffix (.pdb or .cif).
  If output.file_name is not specified, the output file will be
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
    -nuclear: use_neutron_distances
    -demandflipallnhqs: add_flip_movers = True
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
    if self.params.output.file_name is None:
      inName = self.data_manager.get_default_model_name()
      suffix = os.path.splitext(os.path.basename(inName))[1]
      if self.params.add_flip_movers:
        pad = 'FH'
      else:
        pad = 'H'
      base = os.path.splitext(os.path.basename(inName))[0] + pad
      self.params.output.file_name = base + suffix
      print('Writing model output to', self.params.output.file_name)

    self.data_manager.has_models(raise_sorry=True)
    if self.params.output.description_file_name is None:
      raise Sorry("Must specify output.description_file_name")

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

      make_sub_header('Optimizing', out=self.logger)
      startOpt = work_clock()
      Optimizers.probePhil = self.params.probe
      opt = Optimizers.FastOptimizer(self.params.add_flip_movers, self.model, probeRadius=0.25,
        altID=self.params.alt_id, preferenceMagnitude=self.params.preference_magnitude)
      doneOpt = work_clock()
      outString += opt.getInfo()
      outString += 'Time to Add Hydrogen = '+str(doneAdd-startAdd)+'\n'
      outString += 'Time to Interpret = '+str(doneInt-startInt)+'\n'
      outString += 'Time to Optimize = '+str(doneOpt-startOpt)+'\n'

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
    suffix = os.path.splitext(self.params.output.file_name)[1]
    if suffix.lower() == ".pdb":
      txt = self.model.model_as_pdb()
    else:
      txt = self.model.model_as_mmcif()
    self.data_manager._write_text("model", txt, self.params.output.file_name)

    print('Wrote', self.params.output.file_name,'and',
      self.params.output.description_file_name, file = self.logger)

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
