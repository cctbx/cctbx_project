from __future__ import absolute_import, division, print_function

import sys
from mmtbx.refinement.real_space import flipbase

if __name__ == "__main__":
  flipbase.run(sys.argv[1:])
