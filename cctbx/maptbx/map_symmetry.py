# -*- coding: utf-8 -*-

from __future__ import division, print_function
import sys
from libtbx import adopt_init_args
from libtbx.utils import null_out
from scitbx.matrix import col

#  map_symmetry
#  tool to identify and evaluate reconstruction symmetry in a map

class map_symmetry:

  def __init__(self,params=None,
      map_data=None,
      crystal_symmetry=None,
      ncs_object=None,
      log=sys.stdout):
    adopt_init_args(self, locals())
    self.cc=None
    self.original_ncs_object=None
    if ncs_object:
      self.original_ncs_object=ncs_object.deep_copy()
    self.ncs_object=ncs_object

    if self.params and self.params.control.verbose:
      self.local_log=log
    else:
      self.local_log=null_out()

  def get_results(self):
    from libtbx import group_args
    return group_args(
     cc = self.cc,
     ncs_object = self.ncs_object
    )

  def clean_up(self):
    pass

  def run(self):

    # Shift the map origin if necessary

    self.shift_origin()
    self.get_resolution()

    print ("Finding symmetry in map",file=self.log)

    from cctbx.maptbx.segment_and_split_map import run_get_ncs_from_map

    new_ncs_obj,ncs_cc,ncs_score=run_get_ncs_from_map(params=self.params,
      map_data=self.map_data,
      crystal_symmetry=self.crystal_symmetry,
      ncs_obj=self.ncs_object,
      out=self.log)

    if not new_ncs_obj:
      print ("\nNo symmetry found..",file=self.log)
      return

    # Now shift back if necessary
    ncs_object=new_ncs_obj.coordinate_offset(
        coordinate_offset=col(self.origin_shift_cart))

    print ("\nFinal symmetry obtained:",file=self.log)
    print ("\nCorrelation of symmetry-related regions: %.2f   Copies: %d " %(
       ncs_cc,ncs_object.max_operators()), file=self.log)

    if self.params.control.verbose:
      ncs_object.display_all(log=self.log)
    # write to output file
    ncs_object.format_all_for_group_specification(
         file_name=self.params.output_files.symmetry_out)
    print ("\nWrote operators in .ncs_spec format to %s" %(
      self.params.output_files.symmetry_out),file=self.log)

    # Final results
    self.ncs_object=ncs_object
    self.cc=ncs_cc

  def get_resolution(self):
    if not self.params.crystal_info.resolution:
      from cctbx.maptbx import d_min_from_map
      self.params.crystal_info.resolution=\
        d_min_from_map(map_data=self.map_data,
         unit_cell=self.crystal_symmetry.unit_cell())
      print ("\nResolution estimated from map is %.1f A " %(
        self.params.crystal_info.resolution),file=self.log)

  def shift_origin(self):

     origin_shift=(
         self.map_data.origin()[0]/self.map_data.all()[0],
         self.map_data.origin()[1]/self.map_data.all()[1],
         self.map_data.origin()[2]/self.map_data.all()[2])
     origin_shift_cart=\
       self.crystal_symmetry.unit_cell().orthogonalize(origin_shift)

     acc=self.map_data.accessor()
     shift_needed = not \
        (self.map_data.focus_size_1d() > 0 and self.map_data.nd() == 3 and
         self.map_data.is_0_based())
     if(shift_needed):
       self.map_data = self.map_data.shift_origin()
       self.origin_frac=origin_shift
       self.origin_shift_cart=origin_shift_cart
     else:
       self.origin_frac=(0.,0.,0.)
       self.origin_shift_cart=(0,0,0)

     if self.ncs_object:
        self.ncs_object=self.ncs_object.coordinate_offset(
        coordinate_offset=-1*col(self.origin_shift_cart))
