from __future__ import absolute_import, division, print_function
from xfel.merging.application.edit.reflection_table_editor import reflection_table_editor
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """Factory class for editing reflection table."""
  @staticmethod
  def from_parameters(params, additional_info=None):
    """Create a column with reduced HKLs; remove unnecessary columns"""
    return [reflection_table_editor(params)]
