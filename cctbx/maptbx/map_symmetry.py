# -*- coding: utf-8 -*-

from __future__ import division, print_function
import sys
from libtbx import  adopt_init_args
from libtbx.utils import Sorry

#  map_symmetry
#  tool to identify and evaluate reconstruction symmetry in a map

class map_symmetry:

  def __init__(self,params=None,
      map_data=None,
      crystal_symmetry=None,
      ncs_object=None,
      log=sys.stdout):
    adopt_init_args(self, locals())

    self.ncs_cc=None
    self.ncs_object=None

    if self.params and self.params.control.verbose:
      self.local_log=log
    else:
      self.local_log=null_out()

  def clean_up(self):
    pass

  def run(self):
    print ("Finding symmetry in map",file=self.log)
