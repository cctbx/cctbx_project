from __future__ import division, print_function
from iotbx.pdb import modified_aa_names
from iotbx.pdb.amino_acid_codes import three_letter_given_one_letter

def get_aa_parent(code):
  one = modified_aa_names.lookup.get(code.upper(), False)
  if not one: return code
  return three_letter_given_one_letter.get(one, None)
