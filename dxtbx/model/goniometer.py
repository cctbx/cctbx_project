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

import math
import pycbf
from scitbx import matrix

from goniometer_helpers import cbf_gonio_to_effective_axis_fixed

class goniometer:
    '''A class to represent the rotation axis for a standard rotation
    geometry diffraction data set.'''

    def __init__(self, axis, fixed):
        '''Initialize the goniometer, with the real rotation axis (in the CBF
        coordinate frame) and a fixed additional rotation as a matrix which
        represents the additional fixed rotation of the sample attached to the
        rotating axis - for example the effects of kappa and phi for an omega
        scan on a kappa goniometer. This should be a list of 9 floating point
        values:

        fixed = (f11, f12, f13, f21, f22, f23, f31, f31, f33) =

        f11 f12 f13
        f21 f22 f23
        f31 f32 f33

        where this is a rotation matrix which should be applied as

        A = [R][F][U][B]

        in a standard orientation matrix.'''

        assert(len(axis) == 3)
        assert(len(fixed) == 9)

        self._axis = matrix.col(axis)
        self._fixed = matrix.sqr(fixed)

        return

    def __repr__(self):
        '''Generate a useful-to-print representation.'''

        f_axis = '%6.3f %6.3f %6.3f\n'
        f_fixed = 3 * f_axis

        return f_axis % self._axis.elems + f_fixed % self._fixed.elems

    def __cmp__(self, other):
        '''Compare this rotation axis with another.'''

        angle = self._axis.angle(other.get_axis_c())

        if angle < -1.0e-6:
            return -1
        elif angle > 1.0e-6:
            return 1

        return 0

    def get_axis(self):
        '''Get the values for the rotation axis.'''

        return self._axis.elems

    def get_axis_c(self):
        '''Return a cctbx vector for the rotation axis.'''

        return self._axis

    def get_fixed(self):
        '''Return the elements for the fixed rotation matrix.'''

        return self._fixed.elems

    def get_fixed_c(self):
        '''Return the cctbx matrix for the fixed rotation.'''

        return self._fixed

class goniometer_factory:
    '''A factory class for goniometer objects, which will encapsulate
    some standard goniometer designs to make it a little easier to get
    started with all of this - for cases when we are not using a CBF.
    When we have a CBF just use that factory method and everything will be
    peachy.'''

    def __init__(self):
        pass

    @staticmethod
    def single_axis():
        '''Construct a single axis goniometer which is canonical in the
        CBF reference frame.'''

        axis = (1, 0, 0)
        fixed = (1, 0, 0, 0, 1, 0, 0, 0, 1)

        return goniometer(axis, fixed)

    @staticmethod
    def known_axis(axis):
        '''Return an goniometer instance for a known rotation axis, assuming
        that nothing is known about the fixed element of the rotation axis.'''

        assert(len(axis) == 3)

        fixed = (1, 0, 0, 0, 1, 0, 0, 0, 1)

        return goniometer(axis, fixed)

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

        assert(direction in ['-z', '+y', '+z', '-y'])
        assert(scan_axis in ['phi', 'omega'])

        _omega = matrix.col((1, 0, 0))

        c = math.cos(alpha * math.pi / 180.0)
        s = math.sin(alpha * math.pi / 180.0)

        if direction == '-z':
            _kappa = matrix.col((c, 0.0, -s))
        elif direction == '+z':
            _kappa = matrix.col((c, 0.0, s))
        elif direction == '-y':
            _kappa = matrix.col((c, -s, 0.0))
        elif direction == '+y':
            _kappa = matrix.col((c, s, 0.0))

        _phi = matrix.col((1, 0, 0))

        if scan_axis == 'omega':

            K = _kappa.axis_and_angle_as_r3_rotation_matrix(kappa, deg = True)
            P = _phi.axis_and_angle_as_r3_rotation_matrix(phi, deg = True)

            return goniometer(_omega.elems, (K * P).elems)

        elif scan_axis == 'phi':

            O = _omega.axis_and_angle_as_r3_rotation_matrix(omega, deg = True)
            K = _kappa.axis_and_angle_as_r3_rotation_matrix(kappa, deg = True)
            I = (1, 0, 0, 0, 1, 0, 0, 0, 1)

            return goniometer(O * K * _phi.elems, I)

        return

    @staticmethod
    def imgCIF(cif_file):
        '''Initialize a goniometer model from an imgCIF file.'''

        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

        cbf_gonio = cbf_handle.construct_goniometer()

        axis, fixed = cbf_gonio_to_effective_axis_fixed(cbf_gonio)

        cbf_gonio.__swig_destroy__(cbf_gonio)
        del(cbf_gonio)

        return goniometer(axis.elems, fixed.elems)

    @staticmethod
    def imgCIF_H(cbf_handle):
        '''Initialize a goniometer model from an imgCIF file handle, where
        it is assumed that the file has already been read.'''

        cbf_gonio = cbf_handle.construct_goniometer()

        axis, fixed = cbf_gonio_to_effective_axis_fixed(cbf_gonio)

        cbf_gonio.__swig_destroy__(cbf_gonio)
        del(cbf_gonio)

        return goniometer(axis.elems, fixed.elems)
