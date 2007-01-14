from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses
import libtbx.phil.command_line
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
import sys, os

class twin_data(object):
  def __init__(self,
               miller_array,
               phil_block,
               out=None):
    if out is None:
      out=sys.stdout

    print >> out, "Artifically twinning the data"
    print >> out

    self.miller_array = miller_array.deep_copy().set_observation_type(
      miller_array ).map_to_asu()

    self.twin_law = phil_block.scaling.input.optional.twinning.twin_law
    assert (self.twin_law is not None)
    self.twin_law=sgtbx.rt_mx(self.twin_law, r_den=24,t_den=288 )
    if self.twin_law.r().determinant() != 1:
      raise Sorry("The determinant of the provided twin law is not equal to unity")
    self.alpha = phil_block.scaling.input.optional.twinning.fraction
    assert (self.alpha is not None)
    assert self.alpha<0.5
    assert self.alpha>=0.0

    self.mtzout = phil_block.scaling.input.optional.hklout
    assert (self.mtzout is not None)

    # make sure we have intensities
    if self.miller_array.is_real_array():
      if not self.miller_array.is_xray_intensity_array():
        self.miller_array = self.miller_array.f_as_f_sq()
    assert self.miller_array.is_xray_intensity_array()

    cb_op = sgtbx.change_of_basis_op( self.twin_law )

    self.new_miller = self.miller_array.change_basis( cb_op ).map_to_asu()

    xa,xb = self.miller_array.common_sets( self.new_miller )

    new_data = (1.0-self.alpha)*xa.data() + self.alpha*xb.data()
    xa = xa.customized_copy(data=new_data,
                            sigmas=new_data/100.0)

    mtz_dataset = xa.as_mtz_dataset(
      column_root_label='I_TWIN')
    mtz_dataset.mtz_object().write(
      file_name=phil_block.scaling.input.optional.hklout)


















class detwin_data(object):
  def __init__(self,
               miller_array,
               phil_block,
               out=None,
               b_wilson=30.0
               ):
    if out is None:
      out=sys.stdout

    print >> out, "Detwinning data with:"
    print >> out, "  - twin law:      %s"%(phil_block.scaling.input.optional.twinning.twin_law)
    print >> out, "  - twin fraction: %3.2f"%(phil_block.scaling.input.optional.twinning.fraction)
    print >> out
    print >> out, "BE WARNED! DETWINNING OF DATA DOES NOT SOLVE YOUR TWINNING PROBLEM!"
    print >> out, "PREFERABLY, REFINEMENT SHOULD BE CARRIED OUT AGAINST ORIGINAL DATA "
    print >> out, "ONLY USING A TWIN SPECIFIC TARGET FUNCTION!"
    print >> out

    self.miller_array = miller_array.deep_copy().set_observation_type(
      miller_array )

    self.twin_law = phil_block.scaling.input.optional.twinning.twin_law
    assert (self.twin_law is not None)
    self.twin_law=sgtbx.rt_mx(self.twin_law, r_den=24,t_den=288 )
    if self.twin_law.r().determinant() != 1:
      raise Sorry("The determinant of the provided twin law is not equal to unity")

    self.alpha = phil_block.scaling.input.optional.twinning.fraction
    assert (self.alpha is not None)
    assert self.alpha<0.5
    assert self.alpha>=0.0

    self.mtzout = phil_block.scaling.input.optional.hklout
    assert (self.mtzout is not None)

    # make sure we have intensities
    if self.miller_array.is_real_array():
      if not self.miller_array.is_xray_intensity_array():
        self.miller_array = self.miller_array.f_as_f_sq()
    assert self.miller_array.is_xray_intensity_array()

    # add b value
    factor = flex.exp( -b_wilson*self.miller_array.d_star_sq().data()/2.0 )
    self.miller_array = self.miller_array.customized_copy(
                          data = self.miller_array.data()*factor,
                          sigmas = self.miller_array.sigmas()*factor
                        ).set_observation_type( self.miller_array )

    detwin_object = mmtbx.scaling.detwin(self.miller_array.indices(),
                                         self.miller_array.data(),
                                         self.miller_array.sigmas(),
                                         self.miller_array.space_group(),
                                         self.miller_array.anomalous_flag(),
                                         self.twin_law.r().as_double() )

    detwin_object.detwin_with_alpha( self.alpha )

    new_intensities = detwin_object.detwinned_i()
    new_sigmas = detwin_object.detwinned_sigi()
    new_hkl = detwin_object.detwinned_hkl()


    self.miller_array =  self.miller_array.customized_copy(
      indices = new_hkl,
      data =  new_intensities,
      sigmas = new_sigmas ).set_observation_type( self.miller_array )


    mtz_dataset = self.miller_array.as_mtz_dataset(
      column_root_label='I_DETWIN')
    mtz_dataset.mtz_object().write(
      file_name=phil_block.scaling.input.optional.hklout)
