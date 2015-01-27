from __future__ import division
#!/usr/bin/env python
# goniometer.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A model for the goniometer for the "updated experimental model" project
# documented in internal ticket #1555. This is not designed to be used outside
# of the XSweep classes.

import pycbf
from dxtbx_model_ext import Goniometer, KappaGoniometer

from goniometer_helpers import cbf_gonio_to_effective_axis_fixed

class goniometer_factory:
  '''A factory class for goniometer objects, which will encapsulate
  some standard goniometer designs to make it a little easier to get
  started with all of this - for cases when we are not using a CBF.
  When we have a CBF just use that factory method and everything will be
  peachy.'''

  def __init__(self):
    pass

  @staticmethod
  def make_goniometer(rotation_axis, fixed_rotation):
    return Goniometer(
        tuple(map(float, rotation_axis)),
        tuple(map(float, fixed_rotation)))

  @staticmethod
  def make_kappa_goniometer(alpha, omega, kappa, phi, direction, scan_axis):
    return KappaGoniometer(
        float(alpha),
        float(omega),
        float(kappa),
        float(phi),
        str(direction),
        str(scan_axis))

  @staticmethod
  def single_axis():
    '''Construct a single axis goniometer which is canonical in the
    CBF reference frame.'''

    axis = (1, 0, 0)
    fixed = (1, 0, 0, 0, 1, 0, 0, 0, 1)

    return goniometer_factory.make_goniometer(axis, fixed)

  @staticmethod
  def single_axis_reverse():
    '''Construct a single axis goniometer which is canonical in the
    CBF reference frame, but reversed in rotation.'''

    axis = (-1, 0, 0)
    fixed = (1, 0, 0, 0, 1, 0, 0, 0, 1)

    return goniometer_factory.make_goniometer(axis, fixed)

  @staticmethod
  def known_axis(axis):
    '''Return an goniometer instance for a known rotation axis, assuming
    that nothing is known about the fixed element of the rotation axis.'''

    assert(len(axis) == 3)

    fixed = (1, 0, 0, 0, 1, 0, 0, 0, 1)

    return Goniometer(axis, fixed)

  @staticmethod
  def kappa(alpha, omega, kappa, phi, direction, scan_axis):
    '''Return a kappa goniometer where omega is the primary axis (i,e.
    aligned with X in the CBF coordinate frame) and has the kappa arm
    with angle alpha attached to it, aligned with -z, +y, +z or -y at
    omega = 0, that being the direction, which in turn has phi fixed to it
    which should initially be coincident with omega. We also need to know
    which axis is being used for the scan i.e. phi or omega. All angles
    should be given in degrees. This will work by first constructing the
    rotation axes and then composing them to the scan axis and fixed
    component of the rotation.'''

    return goniometer_factory.make_kappa_goniometer(alpha, omega, kappa,
        phi, direction, scan_axis)

  @staticmethod
  def imgCIF(cif_file):
    '''Initialize a goniometer model from an imgCIF file.'''

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

    cbf_gonio = cbf_handle.construct_goniometer()

    axis, fixed = cbf_gonio_to_effective_axis_fixed(cbf_gonio)

    cbf_gonio.__swig_destroy__(cbf_gonio)
    del(cbf_gonio)

    return goniometer_factory.make_goniometer(axis, fixed)

  @staticmethod
  def imgCIF_H(cbf_handle):
    '''Initialize a goniometer model from an imgCIF file handle, where
    it is assumed that the file has already been read.'''

    cbf_gonio = cbf_handle.construct_goniometer()

    axis, fixed = cbf_gonio_to_effective_axis_fixed(cbf_gonio)

    cbf_gonio.__swig_destroy__(cbf_gonio)
    del(cbf_gonio)

    return goniometer_factory.make_goniometer(axis, fixed)
