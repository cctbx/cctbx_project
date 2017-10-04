from __future__ import division
import iotbx.pdb.utils

def exercise_all_chain_ids():
  ids = iotbx.pdb.utils.all_chain_ids()
  assert len(ids)==3844

def run():
  exercise_all_chain_ids()
  print "OK"

if __name__ == '__main__':
  run()
