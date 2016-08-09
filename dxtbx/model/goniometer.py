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
from dxtbx_model_ext import Goniometer, KappaGoniometer, MultiAxisGoniometer

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
  def make_multi_axis_goniometer(axes, angles, scan_axis):
    return MultiAxisGoniometer(axes, angles, scan_axis)

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
  def multi_axis(axes, angles, scan_axis):
    '''
    '''

    return goniometer_factory.make_multi_axis_goniometer(
      axes, angles, scan_axis)

  @staticmethod
  def imgCIF(cif_file):
    '''Initialize a goniometer model from an imgCIF file.'''

    # FIXME in here work out how to get the proper setting matrix if != 1

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

    return goniometer_factory.imgCIF_H(cbf_handle)

  @staticmethod
  def imgCIF_H(cbf_handle):
    '''Initialize a goniometer model from an imgCIF file handle, where
    it is assumed that the file has already been read.'''

    # find the goniometer axes and dependencies
    from scitbx.array_family import flex
    axis_names = flex.std_string()
    depends_on = flex.std_string()
    axes = flex.vec3_double()
    angles = flex.double()
    scan_axis = None
    cbf_handle.find_category("axis")
    for i in range(cbf_handle.count_rows()):
      cbf_handle.find_column("equipment")
      if cbf_handle.get_value() == "goniometer":
        cbf_handle.find_column("id")
        axis_names.append(cbf_handle.get_value())
        axis = []
        for i in range(3):
          cbf_handle.find_column("vector[%i]" %(i+1))
          axis.append(float(cbf_handle.get_value()))
        axes.append(axis)
        cbf_handle.find_column("depends_on")
        depends_on.append(cbf_handle.get_value())
        cbf_handle.next_row()

    # find the starting angles of each goniometer axis and figure out which one
    # is the scan axis (i.e. non-zero angle_increment)
    cbf_handle.find_category("diffrn_scan_axis")
    for i in range(cbf_handle.count_rows()):
      cbf_handle.find_column("axis_id")
      axis_name = cbf_handle.get_value()
      if axis_name not in axis_names: continue
      cbf_handle.find_column("angle_start")
      axis_angle = float(cbf_handle.get_value())
      cbf_handle.find_column("angle_increment")
      increment = float(cbf_handle.get_value())
      angles.append(axis_angle)
      if abs(increment) > 0:
        assert scan_axis is None, "More than one scan axis is defined: not currently supported"
        scan_axis = flex.first_index(axis_names, axis_name)
      cbf_handle.next_row()
    assert axes.size() == angles.size()
    assert scan_axis is not None

    # figure out the order of the axes from the depends_on values
    order = flex.size_t()
    for i in range(axes.size()):
      if depends_on[i] == '.':
        o = 0
      else:
        o = flex.first_index(axis_names, depends_on[i])+1
      assert o not in order
      order.append(o)

    # multi-axis gonio requires axes in order as viewed from crystal to gonio base
    # i.e. the reverse of the order we have from cbf header
    order = order.reversed()
    axes = axes.select(order)
    angles = angles.select(order)
    scan_axis = axes.size() - scan_axis - 1

    # construct a multi-axis goniometer
    gonio = goniometer_factory.multi_axis(axes, angles, scan_axis)
    return gonio
