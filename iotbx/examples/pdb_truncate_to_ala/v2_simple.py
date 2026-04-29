from __future__ import absolute_import, division, print_function
import iotbx.pdb
import iotbx.pdb.amino_acid_codes
import sys

def run(args):
  if (len(args) == 0):
    raise RuntimeError("Please specify one or more pdb file names.")
  aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter
  ala_atom_names = set([" N  ", " CA ", " C  ", " O  ", " CB "])
  for file_name in args:
    pdb_obj = iotbx.pdb.input(file_name=file_name)
    hierarchy = pdb_obj.construct_hierarchy()
    hierarchy.overall_counts().show()
    for model in hierarchy.models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for ag in rg.atom_groups():
            if (ag.resname in aa_resnames):
              for atom in ag.atoms():
                if (atom.name not in ala_atom_names):
                  ag.remove_atom(atom=atom)
    output_pdb = "v2_truncated_to_ala_"+file_name
    hierarchy.write_pdb_or_mmcif_file(target_filename=output_pdb)

if (__name__ == "__main__"):
  run(sys.argv[1:])
