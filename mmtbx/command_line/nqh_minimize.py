
from __future__ import absolute_import, division, print_function
import sys
from mmtbx.validation.molprobity import nqh_minimize

if __name__ == "__main__":
  nqh_minimize.run(sys.argv[1:])
