from iotbx.cns.space_group_symbols import cns_format
from cctbx import crystal

uc_param_names = "a b c alpha beta gamma".split()

re_sg_uc = r'sg=\s*(\S+)\s*a=\s*(\S+)\s*b=\s*(\S+)\s*c=\s*(\S+)' \
         + r'\s*alpha=\s*(\S+)\s*beta=\s*(\S+)\s*gamma=\s*(\S+)'

re_uc_sg = r'a=\s*(\S+)\s*b=\s*(\S+)\s*c=\s*(\S+)' \
         + r'\s*alpha=\s*(\S+)\s*beta=\s*(\S+)\s*gamma=\s*(\S+)' \
         + r'\s*sg=\s*(\S+)'

def crystal_symmetry_from_re_match(m, i_sg=1, i_uc=2):
  assert m is not None
  try: unit_cell = [float(m.group(i+i_uc)) for i in xrange(6)]
  except ValueError: return None
  try:
    return crystal.symmetry(
      unit_cell=unit_cell,
      space_group_symbol=m.group(i_sg))
  except RuntimeError:
    return None

def crystal_symmetry_as_sg_uc_list(crystal_symmetry):
  u = crystal_symmetry.unit_cell()
  s = crystal_symmetry.space_group_info()
  if (u is None): uc = ["None"]*6
  else:           uc = ["%.6g" % v for v in u.parameters()]
  sg = cns_format(space_group_info=s)
  if (sg is None): sg = str(s).replace(" ","")
  return sg, uc

def crystal_symmetry_as_sg_uc(crystal_symmetry):
  sg, uc = crystal_symmetry_as_sg_uc_list(crystal_symmetry=crystal_symmetry)
  return "sg="+sg+" "+" ".join(["%s=%s" % kv
    for kv in zip(uc_param_names, uc)])

def crystal_symmetry_as_cns_inp_defines(crystal_symmetry):
  sg, uc = crystal_symmetry_as_sg_uc_list(crystal_symmetry=crystal_symmetry)
  result = ['{===>} sg="%s";' % sg]
  for kv in zip(uc_param_names, uc):
    result.append('{===>} %s=%s;' % kv)
  return result
