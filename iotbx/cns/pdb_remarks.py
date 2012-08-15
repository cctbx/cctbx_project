from __future__ import division
from iotbx.cns.crystal_symmetry_utils import \
  re_sg_uc, crystal_symmetry_from_re_match
import re

def extract_symmetry(pdb_record):
  m = re.match(r'REMARK\s+' + re_sg_uc , pdb_record)
  if (not m): return None
  return crystal_symmetry_from_re_match(m=m)
