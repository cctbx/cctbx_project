from libtbx import adopt_init_args
import cctbx.array_family.flex
from mmtbx.refinement import print_statistics

class monitors(object):
  def __init__(self, params,
                     model,
                     fmodels,
                     model_ref = None,
                     log       = None):
    adopt_init_args(self, locals())
    self.monitor_xray = print_statistics.refinement_monitor(
      params    = params,
      model_ref = model_ref,
      out       = log)
    self.monitor_neutron = None
    if(fmodels.fmodel_n is not None):
      self.monitor_neutron = print_statistics.refinement_monitor(
        params    = params,
        model_ref = model_ref,
        out       = log)

  def collect(self, step,
                    target_weights_xray,
                    target_weights_neutron,
                    rigid_body_shift_accumulator = None):
    self.monitor_xray.collect(
      model                        = self.model,
      fmodel                       = self.fmodels.fmodel_xray(),
      step                         = step,
      tan_b_iso_max                = 0, # XXX clean later
      wilson_b                     = self.model.wilson_b,
      target_weights               = target_weights_xray,
      rigid_body_shift_accumulator = rigid_body_shift_accumulator)
    if(self.monitor_neutron is not None):
      self.monitor_neutron.collect(
        model                        = self.model,
        fmodel                       = self.fmodels.fmodel_neutron(),
        step                         = step,
        tan_b_iso_max                = 0, # XXX clean later
        wilson_b                     = self.model.wilson_b,
        target_weights               = target_weights_neutron,
        rigid_body_shift_accumulator = rigid_body_shift_accumulator)
