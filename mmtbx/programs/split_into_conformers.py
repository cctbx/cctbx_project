"""Split a model with alternative conformations into complete conformers"""
from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate
from mmtbx.building import split_into_conformers

class Program(ProgramTemplate):

  description = """
phenix.split_into_conformers: split a model with alternative conformations
into complete conformers.

Every polymer atom of the output model has an altloc label, and each label
makes a complete copy of the polymer. H atoms and waters are removed. Atoms
of the input model do not move: copies are added at the same positions, with
occupancies divided so that structure factors change as little as possible.
Positions are exchanged between conformers so that each conformer has good
geometry and no clashes, and so that occupancy rises from the first conformer
to the last where this is possible. Ligands and ions are not split: they stay
as they are, and only their altloc labels are changed to follow the polymer.

The program stops if an alternative conformation of a residue is incomplete
and too far from the atoms it would have to share with the other
conformations (see bond_stretch_threshold).

Usage:
  phenix.split_into_conformers model.pdb
  phenix.split_into_conformers model.cif ligand.cif clash_distance=2.0

Ligands that are not in the monomer library need a restraints (CIF) file.
"""

  datatypes = ['model', 'phil', 'restraint']

  master_phil_str = """
split_into_conformers {
  clash_distance = 2.1
    .type = float
    .help = Two atoms of one conformer that are not bonded and are closer \\
            than this distance (A) clash
  bond_stretch_threshold = 0.5
    .type = float
    .help = Bonds longer than the ideal value by more than this (A) are \\
            listed in the report. An incomplete alternative conformation \\
            that would need such a bond stops the program
  fill_occupancy_to_one = False
    .type = bool
    .help = Alternative conformations of a polymer atom whose occupancies \\
            add up to less than 1: if True the missing portion is added \\
            equally to each (structure factors change), if False they are \\
            left as they are. Sums above 1 are always scaled down to 1
  nproc = 1
    .type = int
    .help = Number of processes
}
"""

  def validate(self):
    self.data_manager.has_models(raise_sorry=True, expected_n=1,
      exact_count=True)

  def run(self):
    file_name = self.data_manager.get_default_model_name()
    model = self.data_manager.get_model()
    p = self.params.split_into_conformers
    self.result = split_into_conformers.run(
      model                  = model,
      clash_distance         = p.clash_distance,
      bond_stretch_threshold = p.bond_stretch_threshold,
      fill_occupancy_to_one  = p.fill_occupancy_to_one,
      nproc                  = p.nproc,
      log                    = self.logger)
    split_into_conformers.show(
      report                 = self.result.report,
      log                    = self.logger,
      clash_distance         = p.clash_distance,
      bond_stretch_threshold = p.bond_stretch_threshold)
    self.output_file_name = None
    if(self.result.report.nothing_done): return
    base, ext = os.path.splitext(os.path.basename(file_name))
    prefix = self.params.output.prefix
    if(prefix is None): prefix = base
    suffix = self.params.output.suffix
    if(suffix is None): suffix = "_split"
    if(ext.lower() in [".cif", ".mmcif"]):
      model_str, ext, format = self.result.model.model_as_mmcif(), ".cif", "cif"
    else:
      model_str, ext, format = self.result.model.model_as_pdb(), ".pdb", "pdb"
    self.output_file_name = self.data_manager.write_model_file(
      model_str, filename=prefix+suffix+ext, format=format, overwrite=True)
    print("Model written to:", self.output_file_name, file=self.logger)

  def get_results(self):
    return self.result
