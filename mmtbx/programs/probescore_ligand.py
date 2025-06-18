"""Probescore analysis of ligands
"""
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
from cctbx.array_family import flex
import mmtbx.model
import iotbx.pdb
from libtbx.utils import Sorry

from mmtbx.validation.molprobity import probescore_ligand
import os

class Program(ProgramTemplate):

  description = '''
molprobity.probescore_ligand: wrapper for probescore analysis of ligands
  Allows residue selections using Phenix atom selection syntax
    (see https://www.phenix-online.org/documentation/reference/atom_selections.html)
  Returns human-readable digest, raw probe results, or compact oneline

  For best results, use an input file that already contains correct hydrogens
    and the has_h=True flag.  Otherwise, probescore will attempt to add
    hydrogens for you, and this is not always reliable for ligands.

  Currently only whole-residue selections are supported. If you wish to select
    more or less, like a sidechain or a domain interface, use molprobity.probe
    directly.

  Multiple selections can be made in a single command and each will be reported
    separately.

  A single selection can contain multiple residues, suitable for polysaccharides
    See selection syntax for details.
    Selections of > about 30 residues may not work.

Usage examples:
  Select single ligand:
    molprobity.probescore_ligand model_H.pdb has_h=True "chain A resseq 350"
  Select multiple separate ligands:
    molprobity.probescore_ligand model.pdb "chain A resseq 350" "chain A resseq 352" output=oneline
  Select single multi-residue ligand, like carbohydrate:
    molprobity.probescore_ligand model.pdb "chain A resseq 350:355"
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """
  probescore {
    output = *digest raw oneline
      .type = choice
      .help = '''choose output type'''
    has_h = False
      .type = bool
      .help = '''does input file have hydrogens?'''
    nuclear=False
      .type = bool
      .help = '''hydrogen distance for Reduce and Probe, ecloud by default'''
  }
  atom_selection_program {
    inselection = None
      .type = atom_selection
      .help = what to select
      .multiple = True
  }
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(raise_sorry=True)
    if (self.params.atom_selection_program.inselection is None or
        len(self.params.atom_selection_program.inselection) == 0):
      raise Sorry("Need selections")

  # ---------------------------------------------------------------------------
  def run(self):
    filename = os.path.basename(self.data_manager.get_model_names()[0])
    model = self.data_manager.get_model()
    atoms = model.get_atoms()
    all_bsel = flex.bool(atoms.size(), False)
    selection_string_list = []
    for selection_string in self.params.atom_selection_program.inselection:
      selection_string_list.append(selection_string)
    probescore = probescore_ligand.probescore(model, selection_string_list, self.params.probescore.has_h, nuclear=False, out=self.logger)
    if self.params.probescore.output == "digest":
      probescore.print_as_digest()
    elif self.params.probescore.output == "raw":
      probescore.print_as_raw()
    elif self.params.probescore.output == "oneline":
      probescore.print_as_oneline(filename=filename)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return None
