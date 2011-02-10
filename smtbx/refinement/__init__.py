from cctbx import xray
from smtbx.refinement import constraints, least_squares, restraints
import smtbx.utils

class model(object):

  def from_shelx(cls, *args, **kwds):
    from iotbx.shelx import smtbx_refinement_model_from as _
    return _(cls, *args, **kwds)
  from_shelx = classmethod(from_shelx)

  def from_cif(cls, model, reflections):
    """
    We could try to read in the weighting scheme.
    As for constraints and restraints, the CIF format does not support them
    yet
    """
    from iotbx.reflection_file_reader import any_reflection_file
    xs = xray.structure.from_cif(file_path=model)
    ma = any_reflection_file(reflections).as_miller_arrays()[0]
    assert ma.is_xray_intensity_array()
    return cls(fo_sq=ma, xray_structure=xs,
               constraints=[],
               restraints_manager=restraints.manager(),
               weighting_scheme=least_squares.sigma_weighting())
  from_cif=classmethod(from_cif)

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
    return least_squares.crystallographic_ls(
      self.fo_sq,
      reparametrisation,
      restraints_manager=self.restraints_manager,
      weighting_scheme=self.weighting_scheme)
