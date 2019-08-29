from __future__ import absolute_import, division, print_function
import sys, os
from libtbx import easy_pickle
from scitbx.matrix import sqr, col
from cctbx import crystal # dependency for integration pickles

"""
Script that examines a set of cctbx.xfel integration pickles and writes out their basis vectors
in reciprocal space in gnuplot format.  Useful to test for crystal alignment in a liquid jet.

Similar code was used to make figure 5 in Brewster (2015) Acta D.

Usage: run this script with the path to some directories filled with integration pickles to create
and arrows.p file. Run gnuplot, then enter load "arrows.p".
"""

f = open("arrows.p",'w')

for dirname in sys.argv[1:]:
  for filename in os.listdir(dirname):
    print("Reading", os.path.join(dirname, filename))
    try:
      data = easy_pickle.load(os.path.join(dirname, filename))
    except Exception as e:
      print("Couldn't read", filename)
      continue
    ori = data['current_orientation'][0]
    A = sqr(ori.reciprocal_matrix())
    abasis = A * col((1,0,0))
    bbasis = A * col((0,1,0))
    cbasis = A * col((0,0,1))

    f.write("set arrow to % 6.6f, % 6.6f, % 6.6f lc rgb 'red'  \n"%(abasis[0],abasis[1],abasis[2]))
    f.write("set arrow to % 6.6f, % 6.6f, % 6.6f lc rgb 'green'\n"%(bbasis[0],bbasis[1],bbasis[2]))
    f.write("set arrow to % 6.6f, % 6.6f, % 6.6f lc rgb 'blue' \n"%(cbasis[0],cbasis[1],cbasis[2]))
f.close()
