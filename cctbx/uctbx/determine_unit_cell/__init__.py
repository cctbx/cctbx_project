from __future__ import division, print_function
import boost.python
try:
    boost.python.import_ext("determine_unit_cell_ext")
except ImportError:
    print("Cannot import the boost-bound NCDist module. Are you sure that NCDist.h is in the source tree, and that you have rebuilt?")
from determine_unit_cell_ext import *
