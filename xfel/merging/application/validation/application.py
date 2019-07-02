from __future__ import absolute_import, division, print_function

class application:
  def __init__(self,params):

    self.params = params
    self.application_level_validation()

  def application_level_validation(self):
    assert [self.params.scaling.model, self.params.scaling.unit_cell].count(None) == 1, 'Provide only the model or a unit cell, not both'
    assert [self.params.scaling.model, self.params.scaling.space_group].count(None) == 1, 'Provide only the model or a space group, not both'
