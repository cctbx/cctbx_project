from __future__ import division

# chain_comparison.py
# a tool to compare main-chain from two structures with or without crystal
# symmetry
#

import os,sys
from mmtbx.validation.chain_comparison import run

if __name__=="__main__":
  args=sys.argv[1:]
  run(args=args,out=sys.stdout)
