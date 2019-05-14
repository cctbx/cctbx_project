from __future__ import absolute_import, division, print_function
from xfel.merging.application.modify.polarization import polarization
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for modification of intensites. """
  @staticmethod
  def from_parameters(params, additional_info=None):
    """ Presently, only apply polarization correction """
    return [polarization(params)]
