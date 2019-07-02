
from __future__ import absolute_import, division, print_function
import sys
from mmtbx.validation.molprobity import mp_geo

if __name__ == "__main__":
  mp_geo.run(sys.argv[1:])
