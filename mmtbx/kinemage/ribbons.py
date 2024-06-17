from __future__ import absolute_import, division, print_function
from iotbx.pdb import amino_acid_codes, nucleic_acid_codes
from scitbx.matrix import col

_amino_acid_resnames = sorted(amino_acid_codes.one_letter_given_three_letter.keys())
def _IsStandardResidue(resname):
  return resname.strip().upper() in _amino_acid_resnames

_nucleic_acid_resnames = set(nucleic_acid_codes.rna_one_letter_code_dict.keys()).union(
  set(nucleic_acid_codes.dna_one_letter_code_dict.keys()))
def _IsNucleicAcidResidue(resname):
  return resname.strip().upper() in _nucleic_acid_resnames

def _FindContiguousResiduesByAtomDistances(structure, type_function, desired_atoms, distance_threshold):
  '''Return a list of contiguous nucleic acid residues in the structure based on the distance between P atoms,
  or O5* or O5' if the P atom is not found.
  :param structure: iotbx.pdb.hierarchy.root object holding the structure
  :param type_function: function, a function that takes a residue name and returns True if the residue is of the desired type.
  :param desired_atoms: list of strings, the names of the atoms to use for distance calculations.  The first
  one found in the residue will be used.
  :param distance_threshold: float, the maximum distance between desired atoms to consider them contiguous.
  '''
  contiguous_residues = []

  for model in structure.models():
    for chain in model.chains():
      current_contig_residues = []
      prev_atom = None

      for residue_group in chain.residue_groups():
        # Check if the residue is of the desired type, based on its name
        if not type_function(residue_group.unique_resnames()[0]):
          continue

        # Attempt to find the desired atom in the current residue.  Search all atoms for the
        # first one in the desired_atoms list, then later ones if it is not found.
        atom = None
        for name in desired_atoms:
          name = name.strip().upper()
          for a in residue_group.atoms():
            if a.name.strip().upper() == name:
              atom = a
              break
          if atom is not None:
            break

        # If no atom is found, skip to the next residue
        if atom is None:
          continue

        # If this is the first residue being examined, initialize the list
        if prev_atom is None:
          current_contig_residues.append(residue_group)
        else:
          # Calculate the distance between the current and previous atoms
          current_pos = col(atom.xyz)
          prev_pos = col(prev_atom.xyz)
          distance = (current_pos - prev_pos).length()

          if distance < distance_threshold:
            # If the distance is within the threshold, add to the current list
            current_contig_residues.append(residue_group)
          else:
            # If not, start a new list for the current residue if we have at least two residues
            if len(current_contig_residues) > 1:
              contiguous_residues.append(current_contig_residues)
            current_contig_residues = [residue_group]

        # Update the previous atom to the current one
        prev_atom = atom

      # After iterating through the chain, add any remaining contiguous residues to the main list
      # if there are at least two.  
      if len(current_contig_residues) > 1:
        contiguous_residues.append(current_contig_residues)

  return contiguous_residues

def find_contiguous_protein_residues(structure, distance_threshold=5.0):
  '''Return a list of contiguous protein residues in the structure based on the distance between CA atoms.
  :param structure: iotbx.pdb.hierarchy.root object holding the structure
  :param distance_threshold: float, the maximum distance between CA atoms to consider them contiguous.
  The empirical 5.0 default value comes from the Richardson lab's Prekin code; the ideal length is 3.80.
  '''
  return _FindContiguousResiduesByAtomDistances(structure, _IsStandardResidue, ["CA"], distance_threshold)

def find_contiguous_nucleic_acid_residues(structure, distance_threshold=10.0):
  '''Return a list of contiguous nucleic acid residues in the structure based on the distance between P atoms,
  or O5* or O5' if the P atom is not found.
  :param structure: iotbx.pdb.hierarchy.root object holding the structure
  :param distance_threshold: float, the maximum distance between CA atoms to consider them contiguous.
  The empirical 10.0 default value comes from the Richardson lab's Prekin code; the ideal length is ~7.
  '''
  return _FindContiguousResiduesByAtomDistances(structure, _IsNucleicAcidResidue, ["P", "O5*", "O5'"], distance_threshold)
