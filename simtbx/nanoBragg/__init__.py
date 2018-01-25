from __future__ import division
import boost.python
import cctbx.uctbx # possibly implicit
ext = boost.python.import_ext("simtbx_nanoBragg_ext")
from simtbx_nanoBragg_ext import *

class _(boost.python.injector, ext.nanoBragg):

  def __getattr__(self,name):
    """assemble miller array of structure factors used to compute spot intensities from the internal C cube array
       how do we specify docstrings for individual overriden members? """
    if name == "Fhkl":
      from cctbx.crystal import symmetry
      cs = symmetry(unit_cell = self.unit_cell_Adeg,space_group="P 1")
      from cctbx.miller import set, array
      indices,data = self.Fhkl_tuple
      mset = set(crystal_symmetry=cs, anomalous_flag=True, indices=indices)
      return array(mset, data=data).set_observation_type_xray_amplitude()

  def __setattr__(self,name,value):
    """use a P1 anomalous=True miller array to initialize the internal C cube array with structure factors for the spot intensities
       how do we specify docstrings for individual overriden members? """
    if name in ["Fhkl"]:
      value=value.expand_to_p1()
      value=value.generate_bijvoet_mates()
      assert value.space_group_info().type().lookup_symbol() == "P 1"
      # handle exception by expanding to P1
      assert value.anomalous_flag() == True
      # handle exception by copying all F(hkl) to F(-h-k-l)
      #assert values are amplitudes # not sure how to guarantee this
      self.unit_cell_Adeg = value.unit_cell()
      #self.mock_up_group = value.space_group()
      #self.mock_up_anomalous_flag = value.anomalous_flag()
      self.Fhkl_tuple = (value.indices(),value.data())
    else:
      super(ext.nanoBragg,self).__setattr__(name,value)
