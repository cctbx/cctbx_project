from cctbx.sgtbx import cosets
from cctbx import sgtbx
import sys

def run(g,h,out=None):
  if out is None:
    out = sys.stdout

  g = sgtbx.space_group_info( g ).group()
  h = sgtbx.space_group_info( h ).group()

  success = False
  try:
    cc = cosets.left_decomposition( h,g)
    cc.show(out=out)
    success = True
  except: pass

  if not success:
    print >> out, "Coset decomposition not successfull."
    print >> out, "Group %s might not be a subgroup of %s"%(sgtbx.space_group_info( group=g ), sgtbx.space_group_info( group=h )   )
    print >> out, "Sorry....."


if __name__ == "__main__":
  run(sys.argv[1], sys.argv[2])
