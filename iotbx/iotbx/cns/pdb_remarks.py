from cctbx import crystal
import re

def extract_symmetry(pdb_record):
  m = re.match(
      r'REMARK\s+sg=\s*(\S+)\s*a=\s*(\S+)\s*b=\s*(\S+)\s*c=\s*(\S+)'
    + r'\s*alpha=\s*(\S+)\s*beta=\s*(\S+)\s*gamma=\s*(\S+)', pdb_record)
  if (not m): return None
  try: unit_cell = [float(m.group(i+2)) for i in xrange(6)]
  except ValueError: return None
  try:
    return crystal.symmetry(
      unit_cell=unit_cell,
      space_group_symbol=m.group(1))
  except RuntimeError:
    return None
