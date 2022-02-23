from __future__ import absolute_import, division, print_function
from mmtbx.pair_interaction import pair_interaction
import os
import libtbx.load_env

def run():
  folder = libtbx.env.find_in_repositories("qrefine/plugin/yoink/dat")
  if folder is None: return # skip test
  assert not os.path.isfile("p__lda.wfc")
  o = pair_interaction.load_wfc("p")
  assert os.path.isfile("p__lda.wfc")

if __name__=="__main__":
  run()
