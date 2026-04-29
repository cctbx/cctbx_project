"""Example of how to convert SHELX lattice symbol to a space group symbol"""
from __future__ import absolute_import, division, print_function
from iotbx import shelx
from cctbx import sgtbx
import sys

def convert(file_object):
  """ Examplify the direct use of the tool from shelx.lexer

  In practice, one is strongly encouraged to make use of the tools
  from shelx.parsers: that is to say, for the task handled here,
  crystal_symmetry_parser (the code to follow just parrots
  the implementation of crystal_symmetry_parser).
  """
  space_group = None
  for command, line in shelx.command_stream(file=file_object):
    cmd, args = command[0], command[-1]
    if cmd == "LATT":
      assert space_group is None
      assert len(args) == 1
      space_group = sgtbx.space_group()
      n = int(args[0])
      if n > 0:
        space_group.expand_inv(sgtbx.tr_vec((0,0,0)))
      z = "*PIRFABC"[abs(n)]
      space_group.expand_conventional_centring_type(z)
    elif cmd == "SYMM":
      assert space_group is not None
      assert len(args) == 1
      s = sgtbx.rt_mx(args[0])
      space_group.expand_smx(s)
    elif cmd == "SFAC":
      return sgtbx.space_group_info(group=space_group)

def run(args):
  if (len(args) == 0):
    space_group_info = convert(file_object=sys.stdin)
    print(space_group_info.type().lookup_symbol())
  else:
    for file_name in args:
      space_group_info = convert(file_object=open(file_name))
      print(space_group_info.type().lookup_symbol())

if (__name__ == "__main__"):
  run(sys.argv[1:])
