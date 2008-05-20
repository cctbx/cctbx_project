from libtbx import adopt_init_args
import cctbx.array_family.flex
from mmtbx.refinement import print_statistics
from mmtbx import utils

class monitors(object):
  def __init__(self, params,
                     model,
                     fmodels,
                     model_ref = None,
                     log       = None,
                     call_back_after_collect = None):
    adopt_init_args(self, locals())
    self.monitor_xray = print_statistics.refinement_monitor(
      params    = params,
      model_ref = model_ref,
      out       = log,
      call_back_after_collect = call_back_after_collect)
    self.monitor_neutron = None
    if(fmodels.fmodel_n is not None):
      self.monitor_neutron = print_statistics.refinement_monitor(
        params    = params,
        model_ref = model_ref,
        out       = log,
        call_back_after_collect = call_back_after_collect)
    self.target_weights = None

  def collect(self, step,
                    fmodels = None,
                    model = None,
                    target_weights = None,
                    rigid_body_shift_accumulator = None):
    if(target_weights is not None):
      self.target_weights = target_weights
    if(model is not None): self.model = model
    if(fmodels is not None): self.fmodels = fmodels
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.xray_structure)
    self.monitor_xray.collect(
      model                        = self.model,
      fmodel                       = self.fmodels.fmodel_xray(),
      step                         = step,
      tan_b_iso_max                = 0, # XXX clean later
      wilson_b                     = self.model.wilson_b,
      target_weights               = self.target_weights,
      rigid_body_shift_accumulator = rigid_body_shift_accumulator)
    if(self.monitor_neutron is not None):
      utils.assert_xray_structures_equal(
        x1 = self.fmodels.fmodel_neutron().xray_structure,
        x2 = self.model.xray_structure)
      self.monitor_neutron.collect(
        model                        = self.model,
        fmodel                       = self.fmodels.fmodel_neutron(),
        step                         = step,
        tan_b_iso_max                = 0, # XXX clean later
        wilson_b                     = self.model.wilson_b,
        target_weights               = self.target_weights,
        rigid_body_shift_accumulator = rigid_body_shift_accumulator)
