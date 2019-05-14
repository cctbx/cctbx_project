from __future__ import absolute_import, division, print_function
from mmtbx.cablam import cablam_fingerprints
#How to write and format a motif fingerprint for cablam

#Put "from __future__ import absolute_import, division, print_function
#Import the cablam_fingerprints module (see above)

#Write a short description of the motif being coded.

#Create an instance of the motif class
replace_this_with_name_of_motif = cablam_fingerprints.motif(
  motif_name = "replace_this_with_name_of_motif",
  residue_names = {"a":"residue1","b":"residue"},
  superpose_order = {"b":["CA","N"],"c":["CA"],"d":["CA"],"e":["CA","OH"]})
#Pass the class a name for the motif (as a string). This will be used in
#  printing, filenameing, etc.
#motif_name is an attribute of the motif class, not something to replace with
#  the name of the motif
#superpose_order defines the atoms from each indexed residue to be used for
#  automated superposition with superpose_pdbs.  Optional unless that feature is
#  desired for the motif.
#Pass the class a dictionary of names for the residues in the motif
#  The keys for this dictionary must correspond to the indices used to identify
#  residues later in the fingerprint.  Some functions may .sort() these keys for
#  printing, so using alphabetization-friendly keys is advised

#Add the first residue to the motif.  The add_residue() method also returns the
#  new residue for easy access.  Here it's named residue1.
residue1 = replace_this_with_name_of_motif.add_residue(
  allowed_resname = [],
  banned_resname = [],
  sequence_move=None,
  bond_move='',
  end_of_motif=False,
  index='')
#The default values for the function parameters are shown.
#
#allowed_resname and banned_resname accept list of 3-character amino acid names
#  e.g. ['GLY','PRO'].  Neither is required.
#The default [] value (empty list) allows all residue types.
#Residue types in banned_resname are disallowed at this position
#
#sequence_move and bond_move describe how to get to the next residue in the
#  motif.  One or the other, but not both, is *required* for all but the final
#  residue.
#sequence_move accepts a position or negative integer. This number represents
#  the sequence relationship to the next residue.  E.g. sequence_move=2 would
#  result in a move two residues forward (towards C-term) in sequence
#bond_move can be used when the sequence relationship of the destination residue
#  is not known.  bond_move accepts a character string.  This string *must*
#  match the index of a residue already found by the motif (see add_bond())
#
#end_of_motif is a necessary flag for the last residue in the motif.  Set it to
#  True in that case and only that case.
#index accepts a character string matching a key of the residue_names dictionary
#  in the motif object.  Different residues cannot have the same index, unless
#  that index is the default empty string''.

#Add a bond to the residue.  The add_bond() method also returns the new bond for
#  easy access.
bond1 = residue1.add_bond(
  required=True,
  banned=False,
  allow_bifurcated=False,
  src_atom='',
  trg_index='')
#The default values for the function parameters are shown.
#
#required is a boolean flag.  If True, this bond must be present in the motif
#banned is a boolean flag.  If True, this bond must be absent in the motif
#
#allow_bifurcated is a boolean flag.  If true, this Hbond is allowed to be
#  bifurcated.  By default, bifurcated bonds are disallowed.  If both targets of
#  a bifurcated bond are of interest, another add_bond() must be performed for
#  the other target.
#
#src_atom accepts a 4-character string, formatted as a pdb-style atom name. It
#  represents the atom in this residue from which the bond of interest
#  originates. This is *required*.  ' O  ' and ' H  ' are the strings for the
#  protein backbone carbonyl oxygen and amide hydrogen most often involved in
#  hydrogen bonding.
#
#trg_index accepts a character string matching a key of the residue_names
#  dictionary in the motif object.  This index will be assigned to the residue
#  on the target end of this bond (this may be relevant for bond_move). Trying
#  to assign one index to different residues or trying to assign different
#  indices to the same residue is considered a failure to match the motif.
#  trg_index may be left as the default empty string '' without provoking this
#  failure to match, however.

#Add target atoms to the bond.  This is the final step and does *not* return
#  anything.
bond1.add_target_atom(
  atomname=None,
  anyatom=False,
  seqdist=None,
  anyseqdist=False)
#The default values for the function parameters are shown.
#
#A bond may have multiple target atoms.  Each one represents an atom at which
#  the bond in question may terminate.  Typically, there will only be one
#  option, but multiple possible bonding partners can be represented by running
#  add_target_atom multiple times
#
#atomname accepts a 4-character string, formatted as a pdb-style atom name.
#  ' O  ' and ' H  ' are the strings for the protein backbone carbonyl oxygen
#  and amide hydrogen most often involved in hydrogen bonding.
#Alternatively, a bond to any atom, regardless of name, may be allowed by
#  setting anyatom to True
#One of either atomname or anyatom is *required*
#
#seqdist accepts a positive or negative integer and represents the sequence
#  distance to the end of the bond in question.
#Alternatively, a bond to an atom at any sequence separation may be allowed by
#  setting anyseqdist to True.  This is particularly useful in beta structure.
#One of either seqdist or anyseqdist is *required*

#Any number of possible target atoms may be added for each bond.
#Any number of bonds may be required or banned for each residue.
#Any number of residues may be defined for each motif.
#  Order matters when defining residues.  New residues are append()ed to a list,
#  and so have order.  Cablam will search for them in the order they are added.
#  Movement instructions in the form of sequence_move or bond_move are required
#  to get from one residue in this list to the next.

#When finished:
#Check that the last residue has end_of_motif=True
#Check that atom names use pdb-format 4-character names, including whitespace
#Check that the motif_name defined at the top of the motif is non-redundant with
#  any existing motif in cctbx/mmtbx/cablam/fingerprints that you do not wish to
#  overwrite.

#Add this to the bottom of the code:
if __name__ == "__main__":
  cablam_fingerprints.make_pickle(replace_this_with_name_of_motif)
  #Multiple motifs can be stored in the same code, but each needs its own call
  #  to cablam_fingerprints.make_pickle, like so:
  cablam_fingerprints.make_pickle(another_motif)
  #Place all make_pickle calls at the end in a if __name__==__main__ to prevent
  #  spontaneous generation of pickle files if the definitions are imported

#Run the code from the commandline using phenix.python.  This will generate
#  motif files in cctbx/mmtbx/cablam/fingerprints for cablam to use.
#If you make any changes to the motif definition, you must re-run the code to
#  regenerate the pickled motif files with updated information.
