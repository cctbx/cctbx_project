from __future__ import absolute_import, division, print_function
from cctbx import xray
from smtbx.refinement import constraints, least_squares, restraints
import smtbx.utils

class model(object):

  @classmethod
  def from_shelx(cls, *args, **kwds):
    from iotbx.shelx import _smtbx_refinement_model_from as _
    return _(cls, *args, **kwds)

  @classmethod
  def from_cif(cls, model, reflections):
    """
    We could try to read in the weighting scheme.
    As for constraints and restraints, the CIF format does not support them
    yet
    """
    from iotbx.reflection_file_reader import any_reflection_file
    xs_dict = xray.structure.from_cif(file_path=model)
    assert len(xs_dict) == 1, "CIF should contain only one xray structure"
    xs = list(xs_dict.values())[0]
    mas = any_reflection_file(reflections).as_miller_arrays(crystal_symmetry=xs)
    fo_sq = None
    for ma in mas:
      if ma.is_xray_intensity_array() and ma.sigmas() is not None:
        fo_sq = ma.as_xray_observations()
        break
    assert fo_sq is not None
    return cls(fo_sq=fo_sq, xray_structure=xs,
               constraints=[],
               restraints_manager=restraints.manager(),
               weighting_scheme=least_squares.sigma_weighting(),
               wavelength=xs.wavelength)

  def __init__(self, fo_sq, xray_structure,
               constraints, restraints_manager, weighting_scheme,
               temperature_in_celsius=20,
               conformer_indices=None,
               wavelength=None):
    self.fo_sq = fo_sq
    self.xray_structure = xray_structure
    self.constraints = constraints
    self.restraints_manager = restraints_manager
    self.weighting_scheme = weighting_scheme
    self.connectivity_table = smtbx.utils.connectivity_table(
      self.xray_structure,
      conformer_indices=conformer_indices
    )
    self.temperature_in_celsius = temperature_in_celsius
    self.wavelength = wavelength

  def make_anisotropic(self):
    self.xray_structure.convert_to_anisotropic()
    for sc in self.xray_structure.scatterers():
      sc.flags.set_grad_u_aniso(True)
      sc.flags.set_grad_u_iso(False)

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
