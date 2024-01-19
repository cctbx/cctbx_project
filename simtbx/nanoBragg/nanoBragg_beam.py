"""
Organizer for nanoBragg beam properties
"""
from __future__ import print_function, division
from dxtbx.model.beam import BeamFactory
from dxtbx_model_ext import flex_Beam
import numpy as np
from copy import deepcopy
from mmtbx_reduce_ext import RotatePointDegreesAroundAxisDir


def rotate_axis(v, axis, phi):
    """
    :param v: vector to rotate
    :param axis: axis of rotation
    :param phi: angle in radians
    """
    new_v = RotatePointDegreesAroundAxisDir((0,0,0), axis, v, phi*180/np.pi)
    return new_v


class NBbeam(object):

  def __init__(self):
    self.spectrum = [(1.8, 1e12)] # angstroms, photons per pulse
    self.unit_s0 = 1, 0, 0  # forward beam direction
    self.polarization_fraction = 1  # defines horizontal and vertical polarization fraction
    self.divergence_mrad = 0  # set the divergence cone angle
    self.divsteps = 0  # number of divergence steps, will be squared (one per horizontal, vertical directions)
    self.size_mm = 0.001 # beam spot size
    self._undo_nanoBragg_norm_by_nbeams = True # we undo it by default
    self.prev_xray_beams = None  # used to cache most recent xray_beams property for efficiency
    self.num_div_angles_within_cone = 0  # used to count how manydivergence angles we sample within the cone of divergence

  @property
  def divsteps(self):
    return self._divsteps

  @divsteps.setter
  def divsteps(self, val):
    if val > 0:
      assert val % 2 == 0, "divsteps must be even"
    self._divsteps = val

  @property
  def divergences(self):
    divrange = self.divergence_mrad/1000.
    if self.divsteps==0:
      return [(0,0)]
    else:
      all_divs = np.arange(0, divrange+1e-7, divrange / self.divsteps) - divrange / 2
      return [(hdiv, vdiv) for vdiv in all_divs for hdiv in all_divs]

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

    divs = self.divergences

    wavelen = self.spectrum[0][0]
    nominal_beam = BeamFactory.simple(wavelen * 1e-10)
    nominal_beam.set_unit_s0(self.unit_s0)
    nominal_beam.set_polarization_fraction(self.polarization_fraction)
    beam_vector = nominal_beam.get_sample_to_source_direction()
    beam_vector /= np.linalg.norm(beam_vector)
    vert_vector = nominal_beam.get_polarization_normal()
    polar_vector = np.cross(beam_vector, vert_vector)
    polar_vector /= np.linalg.norm(polar_vector)

    self.num_div_angles_within_cone = 0
    beams = []
    for hdiv, vdiv in divs:
      vec_xyz = rotate_axis(-beam_vector, polar_vector, vdiv)
      unit_s0 = rotate_axis(vec_xyz, vert_vector, hdiv)
      div_ang = np.arccos(np.dot(unit_s0, -beam_vector))
      if hdiv == 0 and vdiv == 0:
        assert div_ang == 0
        assert np.allclose(unit_s0, nominal_beam.get_unit_s0())
      if div_ang > 1.1*(self.divergence_mrad / 1000. / 2.):
        continue
      self.num_div_angles_within_cone += 1
      for wavelen, flux in self.spectrum:
        beam = deepcopy(nominal_beam)
        beam.set_wavelength(wavelen*1e-10)
        beam.set_flux(flux)
        beam.set_polarization_fraction(self.polarization_fraction)
        beam.set_unit_s0(unit_s0)
        beam.set_divergence(div_ang)
        beams.append(beam)

    # set normalization
    norm = 1
    if self._undo_nanoBragg_norm_by_nbeams:
      norm = len(beams)
    for beam in beams:
      beam.set_flux(beam.get_flux()/norm)
      self._xray_beams.append(beam)

    self.prev_xray_beams = self._xray_beams
    return self._xray_beams

  @property
  def nanoBragg_constructor_beam(self):
    """dumb necessity for instantiating nanoBragg."""
    if self.prev_xray_beams is None:
      self.prev_xray_beams = self.xray_beams

    beam = BeamFactory.from_dict(self.prev_xray_beams[0].to_dict())

    # set the nominal beam to have the average wavelength and the total flux
    num = 0
    den = 0
    flux = 0
    u0 = []
    div = 0
    count = 0
    for b in self.prev_xray_beams:
      wave = b.get_wavelength()
      wt = b.get_flux()
      num += wave * wt
      den += wt
      flux += b.get_flux()
      div += b.get_divergence()
      u0.append( b.get_unit_s0())
      count += 1
    u0 = np.mean(u0, 0)
    ave_wave = num / den
    div = div / count

    beam.set_divergence(div)
    beam.set_wavelength(ave_wave * 1e10)
    beam.set_flux(flux)
    beam.set_unit_s0(u0)
    return beam

  @property
  def number_of_sources(self):
    return len(self.spectrum)
