"""One-letter and three-letter protein codes"""
from __future__ import absolute_import, division, print_function
one_letter_given_three_letter = {
"ALA": "A",
"ARG": "R",
"ASN": "N",
"ASP": "D",
"CYS": "C",
"GLN": "Q",
"GLU": "E",
"GLY": "G",
"HIS": "H",
"ILE": "I",
"LEU": "L",
"LYS": "K",
"MET": "M",
"MSE": "M",
"PHE": "F",
"PRO": "P",
'PYL': 'O',
'SEC': 'U',
"SER": "S",
"THR": "T",
"TRP": "W",
"TYR": "Y",
"VAL": "V",
"UNK": "X", # Not described in standard (pdb does not define one letter code)
}

one_letter_given_three_letter_modified_aa  = {
# modified or unusual AA
"CSO" : "C", # oxidized Cys
"LLP" : "K", # Lys + PLP
"MLY" : "K", # dimethyllysine
"PTR" : "Y", # phosphotyrosine
"SEP" : "S", # phosphoserine
"TPO" : "T", # phosphothreonine
"TYS" : "Y", # sulfonated tyrosine
# XXX https://lists.sdsc.edu/pipermail/pdb-l/2014-January/005899.html
#"PYL" : "O", # pyrrolysine
#"SEC" : "U", # selenocysteine
}

three_letter_given_one_letter = {
"A": "ALA",
"C": "CYS",
"D": "ASP",
"E": "GLU",
"F": "PHE",
"G": "GLY",
"H": "HIS",
"I": "ILE",
"K": "LYS",
"L": "LEU",
"M": "MET",
"N": "ASN",
'O': 'PYL',
"P": "PRO",
"Q": "GLN",
"R": "ARG",
"S": "SER",
"T": "THR",
'U': 'SEC',
"V": "VAL",
"W": "TRP",
"Y": "TYR",
"X": "UNK", # Not described in standard (pdb does not define one letter code)
}

three_letter_l_given_three_letter_d = {
"DAL": "ALA",
"DAR": "ARG",
"DAS": "ASP",
"DCY": "CYS",
"DGL": "GLU",
"DGN": "GLN",
"DHI": "HIS",
"DIL": "ILE",
"DLE": "LEU",
"DLY": "LYS",
"DPN": "PHE",
"DPR": "PRO",
"DSG": "ASN",
"DSN": "SER",
"DTH": "THR",
"DTR": "TRP",
"DTY": "TYR",
"DVA": "VAL",
"MED": "MET"}

three_letter_d_given_three_letter_l = {
"ALA": "DAL",
"ARG": "DAR",
"ASN": "DSG",
"ASP": "DAS",
"CYS": "DCY",
"GLN": "DGN",
"GLU": "DGL",
"HIS": "DHI",
"ILE": "DIL",
"LEU": "DLE",
"LYS": "DLY",
"MET": "MED",
"PHE": "DPN",
"PRO": "DPR",
"SER": "DSN",
"THR": "DTH",
"TRP": "DTR",
"TYR": "DTY",
"VAL": "DVA"}

def validate_sequence(sequence=None,
                      protein=True, strict_protein=True,
                      nucleic_acid=False, strict_nucleic_acid=True):
  '''

  =============================================================================
  Function for checking if a sequence conforms to the FASTA format.
  (http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml)

  Parameters:
  -----------
  sequence - str (None) - the sequence to be checked
  protein - bool (True) - check for protein letters if True
  strict_protein - bool (True) - only check for the 20 amino acids if True
  nucleic_acid - bool (False) - check for nucleic acid letters if True
  strict_nucleic_acid - bool (True) - only check for the 5 base pairs if True

  Return:
  -------
  The function returns a set of unknown letters. If no unknown letters are
  found, the set is of size 0.

  Notes:
  ------
  There is overlap between the letters used to represent amino aicds and the
  letters used to represent nucleic acids. So if both protein and nucleic_acid
  are set to True, the set of valid letters is the union of both. This will
  make the overall validation less strict.

  =============================================================================

  '''

  # construct set of FASTA letters to test against
  fasta_format = set()
  if (protein):
    fasta_format = fasta_format.union(
      set(three_letter_given_one_letter.keys()))
    #fasta_format.remove('U')     # non-standard letter
    if (not strict_protein):
      fasta_format = fasta_format.union(set(['B', 'U', 'Z', 'X', '*', '-']))
  if (nucleic_acid):
    fasta_format = fasta_format.union(set(['A', 'T', 'C', 'G', 'U']))
    if (not strict_nucleic_acid):
      fasta_format = fasta_format.union(set(['N', 'K', 'M', 'B', 'V', 'S', 'W',
                                             'D', 'Y', 'R', 'H', '-']))

  # test sequence
  unknown_letters = set(sequence.upper()).difference(fasta_format)

  return unknown_letters
