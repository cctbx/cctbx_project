from cctbx import xray
from crys3d.qttbx import xray_structure_viewer
try:
  import durham_structures
except ImportError:
  durham_structures = None
import sys, os

name = sys.argv[1]
filename = os.path.expanduser(name)
if os.path.exists(filename):
  xs = xray.structure.from_shelx(filename=filename, strictly_shelxl=False)
elif durham_structures:
  xs = durham_structures.some([name]).next().xray_structure
else:
  print "%s not found" % name
  sys.exit(1)
xray_structure_viewer.display(xray_structure=xs, name=name)
