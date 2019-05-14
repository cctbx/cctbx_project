from __future__ import absolute_import, division, print_function
from xfel.merging.application.input.file_loader import simple_file_loader
from xfel.merging.application.worker import factory as factory_base

""" Factory class for file loading """

class factory(factory_base):
  @staticmethod
  def from_parameters(params, additional_info=None):
    """ Only one kind of loading supported at present, so construct a simple file loader """
    return [simple_file_loader(params)]
