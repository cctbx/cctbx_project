from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args

class monitors(object):
  def __init__(
        self,
        params,
        model,
        fmodels,
        log = None,
        neutron_refinement = None,
        call_back_handler = None):
    adopt_init_args(self, locals())
    from mmtbx.refinement import print_statistics
    self.monitor_xray = print_statistics.refinement_monitor(
      params    = params,
      neutron_refinement = self.neutron_refinement,
      out       = log,
      call_back_handler = call_back_handler,
      is_neutron_monitor = False)
    self.monitor_neutron = None
    if(fmodels.fmodel_n is not None):
      self.monitor_neutron = print_statistics.refinement_monitor(
        params    = params,
        neutron_refinement = self.neutron_refinement,
        out       = log,
        call_back_handler = call_back_handler,
        is_neutron_monitor = True)

  def collect(
        self,
        step,
        fmodels = None,
        model = None,
        rigid_body_shift_accumulator = None):
    from mmtbx import utils
    if(model is not None): self.model = model
    if(fmodels is not None): self.fmodels = fmodels
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.get_xray_structure())
    self.monitor_xray.collect(
      model                        = self.model,
      fmodel                       = self.fmodels.fmodel_xray(),
      step                         = step,
      wilson_b                     = self.fmodels.fmodel_xray().wilson_b(),
      rigid_body_shift_accumulator = rigid_body_shift_accumulator)
    if(self.monitor_neutron is not None):
      utils.assert_xray_structures_equal(
        x1 = self.fmodels.fmodel_neutron().xray_structure,
        x2 = self.model.get_xray_structure())
      self.monitor_neutron.collect(
        model                        = self.model,
        fmodel                       = self.fmodels.fmodel_neutron(),
        step                         = step,
        wilson_b                     = self.fmodels.fmodel_neutron().wilson_b(),
        rigid_body_shift_accumulator = rigid_body_shift_accumulator)
