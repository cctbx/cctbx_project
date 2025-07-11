from __future__ import absolute_import, division, print_function
import iotbx.pdb
import iotbx.pdb.amino_acid_codes
import sys, os

def run(args):
  if (len(args) == 0):
    raise RuntimeError("Please specify one or more pdb file names.")
  aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter
  ala_atom_names = set([" N  ", " CA ", " C  ", " O  ", " CB "])
  for file_name in args:
    pdb_obj = iotbx.pdb.input(file_name=file_name)
    hierarchy = pdb_obj.construct_hierarchy()
    hierarchy.overall_counts().show()
    n_amino_acid_residues = 0
    n_other_residues = 0
    n_atoms_removed = 0
    for model in hierarchy.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          def have_amino_acid():
            for ag in rg.atom_groups():
              if (ag.resname in aa_resnames):
                return True
            return False
          if (not have_amino_acid()):
            n_other_residues += 1
          else:
            n_amino_acid_residues += 1
            for ag in rg.atom_groups():
              for atom in ag.atoms():
                if (atom.name not in ala_atom_names):
                  ag.remove_atom(atom=atom)
                  n_atoms_removed += 1
    print("Number of amino acid residues:", n_amino_acid_residues)
    print("Number of other residues:", n_other_residues)
    print("Number of atoms removed:", n_atoms_removed)
    if (n_atoms_removed != 0):
      output_pdb = "v4_truncated_to_ala_"+os.path.basename(file_name)
      if (output_pdb.endswith(".gz")): output_pdb = output_pdb[:-3]
      print("Writing file:", output_pdb)
      hierarchy.write_pdb_or_mmcif_file(
        target_filename=output_pdb,
        crystal_symmetry=pdb_obj.crystal_symmetry(),
        append_end=True)
    print()

if (__name__ == "__main__"):
  run(sys.argv[1:])
