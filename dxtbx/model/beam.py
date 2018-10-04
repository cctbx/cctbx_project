from __future__ import absolute_import, division
#!/usr/bin/env python
# beam.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A model for the beam for the "updated experimental model" project documented
# in internal ticket #1555. This is not designed to be used outside of the
# XSweep classes.

import math
import pycbf
from dxtbx_model_ext import Beam
import libtbx.phil

beam_phil_scope = libtbx.phil.parse('''
  beam
    .expert_level = 1
    .short_caption = "Beam overrides"
  {

    wavelength = None
      .type = float
      .help = "Override the beam wavelength"

    direction = None
      .type = floats(size=3)
      .help = "Override the sample to source direction"
      .short_caption = "Sample to source direction"

    polarization_normal = None
      .type = floats(size=3)
      .help = "Override the polarization normal"
      .short_caption = "Polarization normal"

    polarization_fraction = None
      .type = float(value_min=0.0, value_max=1.0)
      .help = "Override the polarization fraction"
      .short_caption = "Polarization fraction"

  }
''')


class BeamFactory:
  '''A factory class for beam objects, which encapsulate standard beam
  models. In cases where a full cbf description is available this
  will be used, otherwise simplified descriptions can be applied.'''

  def __init__(self):
    pass

  @staticmethod
  def from_phil(params, reference=None):
    '''
    Convert the phil parameters into a beam model

    '''
    # Check the input
    if reference is None:
      beam = Beam()
    else:
      beam = reference

    # Set the parameters
    if params.beam.wavelength is not None:
      beam.set_wavelength(params.beam.wavelength)
    elif reference is None:
      raise RuntimeError("No wavelength set")
    if params.beam.direction is not None:
      beam.set_direction(params.beam.direction)
    elif reference is None:
      raise RuntimeError("No beam direction set")
    if params.beam.polarization_normal is not None:
      beam.set_polarization_normal(params.beam.polarization_normal)
    if params.beam.polarization_fraction is not None:
      beam.set_polarization_fraction(params.beam.polarization_fraction)

    # Return the model
    return beam

  @staticmethod
  def from_dict(d, t=None):
    ''' Convert the dictionary to a beam model

    Params:
        d The dictionary of parameters
        t The template dictionary to use

    Returns:
        The beam model

    '''
    from dxtbx.model import Beam

    # If None, return None
    if d == None:
      if t == None: return None
      else: return from_dict(t, None)
    elif t != None:
      d = dict(t.items() + d.items())

    # Create the model from the dictionary
    return Beam.from_dict(d)


  @staticmethod
  def make_beam(sample_to_source=None, wavelength=None,
                s0=None, unit_s0=None,
                divergence=None, sigma_divergence=None):

    if divergence is None or sigma_divergence is None:
      divergence = 0.0
      sigma_divergence = 0.0

    if sample_to_source:
      assert(wavelength)
      return Beam(
          tuple(map(float, sample_to_source)),
          float(wavelength),
          float(divergence),
          float(sigma_divergence))
    elif unit_s0:
      assert(wavelength)
      return Beam(
          tuple(map(lambda x: - x, map(float, unit_s0))),
          float(wavelength),
          float(divergence),
          float(sigma_divergence))
    else:
      assert(s0)
      return Beam(tuple(map(float, s0)))

  @staticmethod
  def make_polarized_beam(sample_to_source=None, wavelength=None,
                          s0=None, unit_s0=None,
                          polarization=None, polarization_fraction=None,
                          divergence=None, sigma_divergence=None,
                          flux=None, transmission=None):
    assert(polarization)
    assert polarization_fraction >= 0.0
    assert polarization_fraction <= 1.0

    if divergence == None or sigma_divergence == None:
      divergence = 0.0
      sigma_divergence = 0.0

    if flux is None:
      flux = 0
    if transmission is None:
      transmission = 1.0

    if sample_to_source:
      assert(wavelength)
      return Beam(
          tuple(map(float, sample_to_source)),
          float(wavelength),
          float(divergence),
          float(sigma_divergence),
          tuple(map(float, polarization)),
          float(polarization_fraction),
          float(flux),
          float(transmission))
    elif unit_s0:
      assert(wavelength)
      return Beam(
          tuple(map(lambda x: - x, map(float, unit_s0))),
          float(wavelength),
          float(divergence),
          float(sigma_divergence),
          tuple(map(float, polarization)),
          float(polarization_fraction),
          float(flux),
          float(transmission))
    else:
      assert(s0)
      return Beam(
          tuple(map(float, s0)),
          float(divergence),
          float(sigma_divergence),
          tuple(map(float, polarization)),
          float(polarization_fraction),
          float(flux),
          float(transmission))

  @staticmethod
  def simple(wavelength):
    '''Construct a beam object on the principle that the beam is aligned
    with the +z axis, as is quite normal. Also assume the beam has
    polarization fraction 0.999 and is polarized in the x-z plane.'''

    return BeamFactory.make_beam(
        sample_to_source=(0.0, 0.0, 1.0),
        wavelength=wavelength)

  @staticmethod
  def simple_directional(sample_to_source, wavelength):
    '''Construct a beam with direction and wavelength.'''

    return BeamFactory.make_beam(
        sample_to_source=sample_to_source,
        wavelength=wavelength)

  @staticmethod
  def complex(sample_to_source, polarization_fraction,
              polarization_plane_normal, wavelength):
    '''Full access to the constructor for cases where we do know everything
    that we need...'''

    return BeamFactory.make_polarized_beam(
        sample_to_source=sample_to_source,
        wavelength=wavelength,
        polarization=polarization_plane_normal,
        polarization_fraction=polarization_fraction)

  @staticmethod
  def imgCIF(cif_file):
    '''Initialize a detector model from an imgCIF file. N.B. the
    definition of the polarization plane is not completely helpful
    in this - it is the angle between the polarization plane and the
    +Y laboratory frame vector.'''

    d2r = math.pi / 180.0

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

    cbf_handle.find_category('axis')

    # find record with equipment = source
    cbf_handle.find_column('equipment')
    cbf_handle.find_row('source')

    # then get the vector and offset from this
    direction = []

    for j in range(3):
      cbf_handle.find_column('vector[%d]' % (j + 1))
      direction.append(cbf_handle.get_doublevalue())

    # and the wavelength
    wavelength = cbf_handle.get_wavelength()

    # and information about the polarization - FIXME this should probably
    # be a rotation about the beam not about the Z axis.

    try:
      polar_fraction, polar_angle = cbf_handle.get_polarization()
    except Exception:
      polar_fraction = 0.999
      polar_angle = 0.0

    polar_plane_normal = (
        math.sin(polar_angle * d2r), math.cos(polar_angle * d2r), 0.0)

    return BeamFactory.make_polarized_beam(
            sample_to_source=direction,
            wavelength=wavelength,
            polarization=polar_plane_normal,
            polarization_fraction=polar_fraction)

  @staticmethod
  def imgCIF_H(cbf_handle):
    '''Initialize a detector model from an imgCIF file. N.B. the
    definition of the polarization plane is not completely helpful
    in this - it is the angle between the polarization plane and the
    +Y laboratory frame vector. This example works from a cbf_handle,
    which is already configured.'''

    d2r = math.pi / 180.0

    cbf_handle.find_category('axis')

    # find record with equipment = source
    cbf_handle.find_column('equipment')
    cbf_handle.find_row('source')

    # then get the vector and offset from this
    direction = []

    for j in range(3):
      cbf_handle.find_column('vector[%d]' % (j + 1))
      direction.append(cbf_handle.get_doublevalue())

    # and the wavelength
    wavelength = cbf_handle.get_wavelength()

    # and information about the polarization - FIXME this should probably
    # be a rotation about the beam not about the Z axis.

    try:
      polar_fraction, polar_angle = cbf_handle.get_polarization()
    except Exception:
      polar_fraction = 0.999
      polar_angle = 0.0

    polar_plane_normal = (
        math.sin(polar_angle * d2r), math.cos(polar_angle * d2r), 0.0)

    return BeamFactory.make_polarized_beam(
        sample_to_source=direction,
        wavelength=wavelength,
        polarization=polar_plane_normal,
        polarization_fraction=polar_fraction)
