from __future__ import division

from mmtbx.conformation_dependent_library.cdl_setup import \
  before_pro_groups, not_before_pro_groups

def distance2(a,b):
  d2 = 0
  for i in range(3):
    d2 += (a.xyz[i]-b.xyz[i])**2
  return d2

def get_c_ca_n(atom_group):
  tmp = []
  outl = []
  for name in ["C", "CA", "N"]:
    for atom in atom_group.atoms():
      if atom.name.strip()==name:
        tmp.append(atom)
        break
    else:
      for atom in atom_group.atoms():
        outl.append(atom.format_atom_record())
      tmp = None
      break
  return tmp, outl

def round_to_ten(d):
  t = int(round((float(d))/10))*10
  if t==180: return -180
  return t

def get_res_type_group(resname1, resname2):
  if resname2=="PRO":
    lookup = before_pro_groups
  else:
    lookup = not_before_pro_groups
  for key in lookup:
    if resname1 in lookup[key]:
      return key
  return None




