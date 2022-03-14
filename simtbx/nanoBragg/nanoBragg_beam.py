"""
Organizer for nanoBragg beam properties
"""
from __future__ import print_function, division
from dxtbx.model.beam import BeamFactory
from dxtbx_model_ext import flex_Beam


class NBbeam(object):

  def __init__(self):
    self.spectrum = [(1.8, 1e12)] # angstroms, photons per pulse
    self.unit_s0 = 1, 0, 0
    self.polarization_fraction = 1
    self.divergence = 0
    self.size_mm = 0.001
    self._undo_nanoBragg_norm_by_nbeams = True # we undo it by default

  @property
  def size_mm(self):
    return self._size_mm

  @size_mm.setter
  def size_mm(self, val):
    self._size_mm = val

  @property
  def spectrum(self):
    """ list of (wavelength, flux) defining the energy spectrum"""
    return self._spectrum

  @spectrum.setter
  def spectrum(self, val):
    self._spectrum = val

  @property
  def unit_s0(self):
    return self._unit_s0

  @unit_s0.setter
  def unit_s0(self, val):
    self._unit_s0 = val

  @property
  def divergence(self):
    return self._divergence

  @divergence.setter
  def divergence(self, val):
    self._divergence = val

  @property
  def polarization_fraction(self):
    return self._polarization_fraction

  @polarization_fraction.setter
  def polarization_fraction(self, val):
    self._polarization_fraction = val

  @property
  def xray_beams(self):
    self._xray_beams = flex_Beam()
    norm = 1
    if self._undo_nanoBragg_norm_by_nbeams:
      norm = float(len(self.spectrum))

    for wavelen, flux in self.spectrum:
      beam = BeamFactory.simple(wavelen * 1e-10)
      beam.set_flux(flux / norm)
      beam.set_unit_s0(self.unit_s0)
      beam.set_polarization_fraction(self.polarization_fraction)
      beam.set_divergence(self.divergence)
      self._xray_beams.append(beam)

    return self._xray_beams

  @property
  def nanoBragg_constructor_beam(self):
    """dumb necessity FIXME please"""

    beam = BeamFactory.from_dict(self.xray_beams[0].to_dict())

    # set the nominal beam to have the average wavelength and the total flux
    num = 0
    den = 0
    flux = 0
    for b in self.xray_beams:
      wave = b.get_wavelength()
      wt = b.get_flux()
      num += wave * wt
      den += wt
      flux += b.get_flux()
    ave_wave = num / den

    beam.set_wavelength(ave_wave * 1e10)
    beam.set_flux(flux)
    return beam

  @property
  def number_of_sources(self):
    return len(self.spectrum)
