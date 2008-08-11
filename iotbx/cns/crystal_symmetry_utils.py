from iotbx.cns.space_group_symbols import cns_format
from cctbx import crystal
import re

re_sg_uc = r'sg=\s*(\S+)\s*a=\s*(\S+)\s*b=\s*(\S+)\s*c=\s*(\S+)' \
         + r'\s*alpha=\s*(\S+)\s*beta=\s*(\S+)\s*gamma=\s*(\S+)'

def crystal_symmetry_from_re_match(m):
  try: unit_cell = [float(m.group(i+2)) for i in xrange(6)]
  except ValueError: return None
  try:
    return crystal.symmetry(
      unit_cell=unit_cell,
      space_group_symbol=m.group(1))
  except RuntimeError:
    return None

def crystal_symmetry_as_sg_uc(crystal_symmetry):
  u = crystal_symmetry.unit_cell()
  s = crystal_symmetry.space_group_info()
  if (u is None):
    uc = "a=None b=None c=None alpha=None beta=None gamma=None"
  else:
    uc = "a=%.6g b=%.6g c=%.6g alpha=%.6g beta=%.6g gamma=%.6g" % u.parameters()
  sg = cns_format(space_group_info=s)
  if (sg is None): sg = str(s).replace(" ","")
  return "sg="+sg+" "+uc
