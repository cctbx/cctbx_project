from __future__ import division
from mmtbx.ligands import rdkit_utils

def any_chemical_input_reader(params):
  pe=params.elbow
  if pe.input.chemical_component:
    molecule = rdkit_utils.mol_from_chemical_component(pe.input.chemical_component)
  elif pe.input.smiles:
    molecule = rdkit_utils.mol_from_smiles(pe.input.smiles)
  return molecule
