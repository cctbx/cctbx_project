"""Tools to identify and evaluate reconstruction symmetry in a map"""

# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
import sys, os
from libtbx import adopt_init_args
from libtbx.utils import null_out
from scitbx.matrix import col


class map_symmetry:
  """Identify and evaluate reconstruction symmetry in a map"""
  def __init__(self,params=None,
      map_data=None,
      map_coeffs=None,
      crystal_symmetry=None,
      ncs_object=None,
      fourier_filter = None,
      log=sys.stdout):
    adopt_init_args(self, locals())
    self.cc=None
    self.score=None
    self.original_ncs_object=None
    if ncs_object:
      self.original_ncs_object=ncs_object.deep_copy()
    self.ncs_object=ncs_object

    if self.params and self.params.control.verbose:
      self.local_log=log
    else:
      self.local_log=null_out()
    if self.params.reconstruction_symmetry.find_ncs_directly and \
       not self.map_coeffs:
      from cctbx.maptbx.segment_and_split_map import get_f_phases_from_map
      self.map_coeffs,dummy=get_f_phases_from_map(map_data=self.map_data,
        crystal_symmetry=self.crystal_symmetry,
        d_min=self.params.crystal_info.resolution,
        return_as_map_coeffs=True, # required
        out=self.log)


  def get_results(self):
    """Get results from map_symmetry"""
    from libtbx import group_args
    if not self.ncs_object:
      from mmtbx.ncs.ncs import ncs
      self.ncs_object=ncs()
      self.score=None
      self.cc=None
    return group_args(
     cc = self.cc,
     ncs_object = self.ncs_object,
     ncs_name = str(self.ncs_object.get_ncs_name()),
     ncs_operators = self.ncs_object.max_operators(),
     score=self.score,
    )

  def clean_up(self):
    """Does nothing"""
    pass

  def run(self):
    """Identify symmetry from map"""

    # Print out values of parameters
    import iotbx.phil
    from phenix.programs.map_symmetry import master_phil_str
    master_phil=iotbx.phil.parse(master_phil_str)
    print ("\nInput parameters for map_symmetry:\n",file=self.log)
    master_phil.format(python_object=self.params).show(out=self.log)

    # Shift the map origin if necessary
    self.shift_origin()
    self.get_resolution()

    print ("Finding symmetry in map",file=self.log)

    if self.params.reconstruction_symmetry.find_ncs_directly:
      new_ncs_obj,ncs_cc,ncs_score=self.find_ncs_from_density()

    else:  # usual
      from cctbx.maptbx.segment_and_split_map import run_get_ncs_from_map
      new_ncs_obj,ncs_cc,ncs_score=run_get_ncs_from_map(params=self.params,
        map_data=self.map_data,
        crystal_symmetry=self.crystal_symmetry,
        ncs_obj=self.ncs_object,
        fourier_filter = self.fourier_filter,
        out=self.log)

    if not new_ncs_obj:
      print ("\nNo symmetry found..",file=self.log)
      return

    # Now shift back if necessary
    ncs_object=new_ncs_obj.coordinate_offset(
        coordinate_offset=col(self.origin_shift_cart))

    print ("\nFinal symmetry obtained:",file=self.log)
    if ncs_object.get_ncs_name():
      print ("NCS type: %s" %(ncs_object.get_ncs_name()),file=self.log)
    print ("Correlation of symmetry-related regions: %.2f   Copies: %d " %(
       ncs_cc,ncs_object.max_operators()), file=self.log)

    if self.params.control.verbose:
      ncs_object.display_all(log=self.log)
    # write to output file
    if self.params.output_files.symmetry_out:
      ncs_object.format_all_for_group_specification(
         file_name=self.params.output_files.symmetry_out)
      print ("\nWrote operators in .ncs_spec format to %s" %(
        self.params.output_files.symmetry_out),file=self.log)

    # Final results
    self.ncs_object=ncs_object
    self.cc=ncs_cc
    self.score=ncs_score

  def find_ncs_from_density(self):
    """Find local symmetry (ncs) from density in map"""
    # First test to make sure the function is available
    try:
      from phenix.command_line.find_ncs_from_density import \
       find_ncs_from_density as find_ncs
    except Exception as e:
      ncs_cc=None
      ncs_object=None
      ncs_score=None
      return ncs_object,ncs_cc,ncs_score

    print ("Finding symmetry in map by search for matching density",
       file=self.log)
    # Write out mtz file to search in...
    temp_dir=self.params.output_files.temp_dir
    if not os.path.isdir(temp_dir):
      os.mkdir(temp_dir)
    map_coeffs_file=os.path.join(temp_dir,"map_coeffs.mtz")
    self.map_coeffs.as_mtz_dataset(
       column_root_label='FWT').mtz_object().write(file_name=map_coeffs_file)


    args=["%s" %(map_coeffs_file),"map_operators_inside_unit_cell=True"]
    if self.params.crystal_info.resolution:
      args.append("resolution=%s" %(self.params.crystal_info.resolution))
    find_ncs_from_density=find_ncs( args,out=self.log)
    if hasattr(find_ncs_from_density,'ncs_object') and \
        find_ncs_from_density.ncs_object:
      ncs_object=find_ncs_from_density.ncs_object
      ncs_cc=ncs_object.overall_cc()
      ncs_score=ncs_cc*(ncs_object.max_operators())**0.5
    else:
      ncs_cc=None
      ncs_object=None
      ncs_score=None

    return ncs_object,ncs_cc,ncs_score
  def get_resolution(self):
    """Use specified resolution or estimate it from map"""
    if not self.params.crystal_info.resolution:
      from cctbx.maptbx import d_min_from_map
      self.params.crystal_info.resolution=\
        d_min_from_map(map_data=self.map_data,
         unit_cell=self.crystal_symmetry.unit_cell())
      print ("\nResolution estimated from map is %.1f A " %(
        self.params.crystal_info.resolution),file=self.log)

  def shift_origin(self):
     """Shift origin to put it at (0,0,0)"""
     origin_shift_grid_units=self.map_data.origin()
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
       # Adjust likely center position
       new_location=[]
       for xx,osc in zip(self.params.reconstruction_symmetry.symmetry_center,
         origin_shift_cart):
         new_location.append(xx-osc)
       self.params.reconstruction_symmetry.symmetry_center=tuple(new_location)
       print("Shifted guess for symmetry center is at: (%.2f,%.2f,%.2f) A "%(
        self.params.reconstruction_symmetry.symmetry_center),file=self.log)
     else:
       self.origin_frac=(0.,0.,0.)
       self.origin_shift_cart=(0,0,0)

     if self.ncs_object:
        self.ncs_object=self.ncs_object.coordinate_offset(
        coordinate_offset=-1*col(self.origin_shift_cart))
