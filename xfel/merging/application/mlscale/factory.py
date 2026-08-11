from __future__ import absolute_import, division, print_function

"""
Factory for the mlscale workers.

cctbx.xfel.merge splits a step name on '_', uses the first token to locate
'xfel.merging.application.<token>.factory', and passes the remaining tokens as
additional_info.  So 'mlscale_step1' and 'mlscale_step2' both resolve here,
with additional_info ['step1'] and ['step2'] respectively.
"""

from xfel.merging.application.worker import factory as factory_base


class factory(factory_base):
  """Factory class for maximum-likelihood scaling and merging."""

  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    if not additional_info:
      raise ValueError(
          "The mlscale factory requires a step suffix. Use 'mlscale_step1' "
          "(annotate observations with partiality geometry) and "
          "'mlscale_step2' (joint ML scaling and merging) in "
          "dispatch.step_list.")

    step = additional_info[0]

    if step == 'step1':
      from xfel.merging.application.mlscale.mlscale_step1 import rho_annotator
      return [rho_annotator(params, mpi_helper, mpi_logger)]

    if step == 'step2':
      from xfel.merging.application.mlscale.mlscale_step2 import ml_merger
      return [ml_merger(params, mpi_helper, mpi_logger)]

    raise ValueError("Unknown mlscale step '%s'. Expected 'step1' or 'step2'."
                     % step)
