
import libtbx.load_env
from libtbx import group_args
import os

def extract_phi_psi (pdb_hierarchy, atom_selection=None) :
  from iotbx.pdb.amino_acid_codes import one_letter_given_three_letter
  from cctbx import geometry_restraints
  if (atom_selection is None) :
    from scitbx.array_family import flex
    atom_selection = flex.bool(pdb_hierarchy.atoms().size(), True)
  rama_angles = []
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for conformer in chain.conformers() :
        residues = conformer.residues()
        for i, residue in enumerate(residues) :
          altloc = ""
          if (not residue.resname in one_letter_given_three_letter) :
            continue
          next_res, prev_res = None, None
          resseq2 = residue.resseq_as_int()
          resseq1, resseq3 = None, None
          if (i > 0):
            resseq1 = residues[i-1].resseq_as_int()
            if (resseq2 != (resseq1 + 1)) :
              continue
            prev_res = residues[i-1]
          if (i < (len(residues) - 1)) :
            resseq3 = residues[i+1].resseq_as_int()
            if (resseq2 != (resseq3 - 1)) :
              continue
            next_res = residues[i+1]
          if (next_res is not None) and (prev_res is not None) :
            c1, n2, ca2, c2, n3 = None, None, None, None, None
            for atom in prev_res.atoms() :
              if (atom.name == " C  ") :
                c1 = atom
                break
            for atom in residue.atoms() :
              if (atom.name == " N  ") :
                n2 = atom
              elif (atom.name == " CA ") :
                ca2 = atom
              elif (atom.name == " C  ") :
                c2 = atom
            for atom in next_res.atoms() :
              if (atom.name == " N  ") :
                n3 = atom
            if (None in [c1, n2, ca2, c2, n3]) :
              #print >> log, "  incomplete backbone for %s %d-%d, skipping." % \
              #  (chain.id, resseq1, resseq3)
              continue
            altloc = ca2.fetch_labels().altloc
            i_seqs = [c1.i_seq,n2.i_seq,ca2.i_seq,c2.i_seq,n3.i_seq]
            for i_seq in i_seqs :
              if (not atom_selection[i_seq]) :
                continue
            pep1 = geometry_restraints.bond(
              sites=[c1.xyz,n2.xyz],
              distance_ideal=1,
              weight=1)
            pep2 = geometry_restraints.bond(
              sites=[c2.xyz,n3.xyz],
              distance_ideal=1,
              weight=1)
            if (pep1.distance_model > 4) or (pep2.distance_model > 4) :
              continue
            residue_name = residue.resname
            if (residue_name == "PRO") :
              residue_type = "pro"
            elif (residue_name == "GLY") :
              residue_type = "gly"
            elif (residues[i+1].resname == "PRO") :
              residue_type = "prepro"
            else :
              residue_type = "ala"
            phi_psi = group_args(
              residue_name=residue_name,
              residue_type=residue_type,
              chain_id=chain.id,
              altloc=altloc,
              resid=residue.resid(),
              i_seqs=i_seqs)
            rama_angles.append(phi_psi)
  return rama_angles

def exercise () :
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if (file_name is None) :
    print "PDB file not found, skipping test."
    return
  from iotbx import file_reader
  pdb_in = file_reader.any_file(file_name, force_type="pdb").file_object
  hierarchy = pdb_in.construct_hierarchy()
  atoms = hierarchy.atoms()
  atoms.reset_i_seq()
  angles = extract_phi_psi(hierarchy)
  assert (len(angles) == 237)
  i_seqs = angles[59].i_seqs
  assert (i_seqs == [457, 466, 467, 468, 477])

if (__name__ == "__main__") :
  exercise()
  print "OK"
