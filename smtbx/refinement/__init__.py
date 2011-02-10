from smtbx.refinement import constraints, least_squares
import smtbx.utils

class model(object):

  def from_shelx(cls, *args, **kwds):
    from iotbx.shelx import smtbx_refinement_model_from as _
    return _(cls, *args, **kwds)
  from_shelx = classmethod(from_shelx)

  def __init__(self, fo_sq, xray_structure,
               constraints, restraints_manager, weighting_scheme):
    self.fo_sq = fo_sq
    self.xray_structure = xray_structure
    self.constraints = constraints
    self.restraints_manager = restraints_manager
    self.weighting_scheme = weighting_scheme
    self.connectivity_table = smtbx.utils.connectivity_table(
      self.xray_structure)

  def least_squares(self):
    reparametrisation = constraints.reparametrisation(
      self.xray_structure,
      self.constraints,
      self.connectivity_table)
    return least_squares.normal_equations(
      self.fo_sq,
      reparametrisation,
      restraints_manager=self.restraints_manager,
      weighting_scheme=self.weighting_scheme)
