from iotbx.shelx.crystal_symmetry_from_ins import read_shelx_latt
from cctbx import sgtbx
import sys

def convert(file_object):
  space_group = None
  for line in file_object:
    l = line.rstrip().split("!")[0]
    if (l.startswith("LATT ")):
      assert space_group is None
      latt = read_shelx_latt(l)
      space_group = sgtbx.space_group()
      if (latt.centric):
        space_group.expand_inv(sgtbx.tr_vec((0,0,0)))
      space_group.expand_conventional_centring_type(latt.z)
    elif (l.startswith("SYMM ")):
      assert space_group is not None
      s = sgtbx.rt_mx(l[5:])
      space_group.expand_smx(s)
  return sgtbx.space_group_info(group=space_group)

def run(args):
  if (len(args) == 0):
    space_group_info = convert(file_object=sys.stdin)
    print space_group_info.type().lookup_symbol()
  else:
    for file_name in args:
      space_group_info = convert(file_object=open(file_name))
      print space_group_info.type().lookup_symbol()

if (__name__ == "__main__"):
  run(sys.argv[1:])
