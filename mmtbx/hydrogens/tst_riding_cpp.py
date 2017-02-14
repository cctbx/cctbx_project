from __future__ import division
import time
from mmtbx_hydrogens_ext import *

# -----------------------------------------------------------------
# This test is temporary, to make sure that intermediate steps
# during conversion to c++ work fine
# -----------------------------------------------------------------

def run():
  print 'Test if c++ class is available in python code.'
# fill in values in c++ object
  rc = riding_coefficients(htype='flat_2neigbs', a0=5, a1=2, a2=3, a3=6,
    a=3.467, b=5.4, h=3.58, n=2, disth=0.887)
# print the values
  print rc.htype, rc.a0, rc.a1, rc.a2, rc.a3, rc.a, rc.b, rc.h,\
    rc.n, rc.disth

# for later tests
#def run():
#  #for use_ideal_bonds_angles in [True, False]:
#  for pdb_str, str_name in zip(pdb_list,pdb_list_name):
#      #print 'pdb_string:', str_name, 'idealize =', idealize
#    exercise(pdb_str=pdb_str)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print "OK. Time:", round(time.time()-t0, 2), "seconds"
