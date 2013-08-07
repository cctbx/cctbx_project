from __future__ import division
#!/usr/bin/env python
#
# dxtbx.serialize.xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Authors: James Parkhurst, Richard Gildea, Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

import sys
from scitbx import matrix
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from dxtbx.model.detector_helpers_types import detector_helpers_types

def to_imageset(input_filename, extra_filename=None):
    '''Get an image set from the xds input filename plus an extra filename

    Params:
        input_filename The XDS.INP file
        extra_filename A (G)XPARM.XDS, INTGRATE.HKL or XDS_ASCII.HKL file

    Returns:
        The imageset

    '''
    from iotbx.xds import xds_inp
    from dxtbx.imageset import ImageSetFactory
    import dxtbx

    # Read the input filename
    handle = xds_inp.reader()
    handle.read_file(input_filename)

    # Get the template
    template = handle.name_template_of_data_frames.replace('?', '#')

    # Create the imageset
    imageset = ImageSetFactory.from_template(template)[0]

    # If an extra filename has been specified, try to load models
    if extra_filename:
        models = dxtbx.load(extra_filename)
        imageset.set_beam(models.get_beam())
        imageset.set_detector(models.get_detector())
        imageset.set_goniometer(models.get_goniometer())

    # Return the imageset
    return imageset


def xds_detector_name(dxtbx_name):
    '''Translate from a xia2 name from the detector library to an XDS detector
    name.'''

    if 'pilatus' in dxtbx_name:
        return 'PILATUS'
    if 'rayonix' in dxtbx_name:
        return 'CCDCHESS'
    if 'adsc' in dxtbx_name:
        return 'ADSC'
    if 'saturn' in dxtbx_name:
        return 'SATURN'
    if 'raxis' in dxtbx_name:
        return 'RAXIS'

    raise RuntimeError, 'detector %s unknown' % dxtbx_name


class to_xds(object):
    '''A class to export contents of a Sweep as XDS.INP or XPARM.XDS.'''

    def __init__(self, sweep):
        self._sweep = sweep

        # detector dimensions in pixels
        self.detector_size = map(int, self.get_detector().get_image_size())
        self.fast, self.slow = self.detector_size

        R = align_reference_frame(
            self.get_detector().get_fast_axis(), (1,0,0),
            self.get_detector().get_slow_axis(), (0,1,0))

        self.imagecif_to_xds_transformation_matrix = R

        self.detector_x_axis = (
            R * matrix.col(self.get_detector().get_fast_axis())).elems
        self.detector_y_axis = (
            R * matrix.col(self.get_detector().get_slow_axis())).elems

        F = R * matrix.col(self.get_detector().get_fast_axis())
        S = R * matrix.col(self.get_detector().get_slow_axis())
        N = F.cross(S)
        self.detector_normal = N.elems

        origin = R * matrix.col(self.get_detector().get_origin())

        centre = -(origin - origin.dot(N) * N)
        x = centre.dot(F)
        y = centre.dot(S)

        self.pixel_size = self.get_detector().get_pixel_size()
        f, s = self.pixel_size
        self.detector_distance = origin.dot(N)
        # Need to add 0.5 because XDS seems to do centroids in fortran coords
        self.detector_origin = x/f + 0.5, y/f + 0.5

        # Beam stuff
        self.wavelength = self.get_beam().get_wavelength()
        self.beam_vector = R * matrix.col(self.get_beam().get_direction())
        # just to make sure it is the correct length
        self.beam_vector = self.beam_vector.normalize() / self.wavelength
        self.beam_vector = (- self.beam_vector).elems

        # Scan and goniometer stuff
        self.starting_frame = self.get_scan().get_image_range()[0]
        self.starting_angle = self.get_scan().get_oscillation()[0]
        self.oscillation_range = self.get_scan().get_oscillation()[1]
        self.rotation_axis = (
            R * matrix.col(self.get_goniometer().get_rotation_axis())).elems
        return

    def get_detector(self):
        return self._sweep.get_detector()

    def get_goniometer(self):
        return self._sweep.get_goniometer()

    def get_beam(self):
        return self._sweep.get_beam()

    def get_scan(self):
        return self._sweep.get_scan()

    def get_template(self):
        return self._sweep.get_template()

    def XDS_INP(self, out=None,
                job_card="XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT"):
        if out is None:
            out = sys.stdout

        sensor = self.get_detector().get_type()
        fast, slow = self.detector_size
        f, s = self.pixel_size
        df = int(1000 * f)
        ds = int(1000 * s)

        # FIXME probably need to rotate by pi about the X axis

        detector = xds_detector_name(
            detector_helpers_types.get(sensor, fast, slow, df, ds))
        trusted = self.get_detector().get_trusted_range()

        print >> out, 'DETECTOR=%s MINIMUM_VALID_PIXEL_VALUE=%d OVERLOAD=%d' % \
              (detector, trusted[0] + 1, trusted[1])

        if detector == 'PILATUS':
            print >> out, 'SENSOR_THICKNESS= 0.32'

        print >> out, 'DIRECTION_OF_DETECTOR_X-AXIS= %.3f %.3f %.3f' % \
              self.detector_x_axis

        print >> out, 'DIRECTION_OF_DETECTOR_Y-AXIS= %.3f %.3f %.3f' % \
              self.detector_y_axis

        print >> out, 'NX=%d NY=%d QX=%.4f QY=%.4f' % (fast, slow, f, s)

        print >> out, 'DETECTOR_DISTANCE= %.3f' % self.detector_distance
        print >> out, 'ORGX= %.1f ORGY= %.1f' % self.detector_origin
        print >> out, 'ROTATION_AXIS= %.3f %.3f %.3f' % \
              self.rotation_axis
        print >> out, 'STARTING_ANGLE= %.3f' % \
              self.starting_angle
        print >> out, 'OSCILLATION_RANGE= %.3f' % \
              self.oscillation_range
        print >> out, 'X-RAY_WAVELENGTH= %.5f' % \
              self.wavelength
        print >> out, 'INCIDENT_BEAM_DIRECTION= %.3f %.3f %.3f' % \
              self.beam_vector

        # FIXME LATER
        if hasattr(self.get_beam(), "get_polarization_fraction"):
            print >> out, 'FRACTION_OF_POLARIZATION= %.3f' % \
                self.get_beam().get_polarization_fraction()
            print >> out, 'POLARIZATION_PLANE_NORMAL= %.3f %.3f %.3f' % \
                self.get_beam().get_polarization_normal()
        print >> out, 'NAME_TEMPLATE_OF_DATA_FRAMES= %s' % \
            self.get_template().replace('#', '?')
        print >> out, 'TRUSTED_REGION= 0.0 1.41'
        for f0, s0, f1, s1 in self.get_detector().get_mask():
            print >> out, 'UNTRUSTED_RECTANGLE= %d %d %d %d' % \
                  (f0 - 1, f1 + 1, s0 - 1, s1 + 1)

        start_end = self.get_scan().get_image_range()

        if start_end[0] == 0:
            start_end = (1, start_end[1])

        print >> out, 'DATA_RANGE= %d %d' % start_end
        print >> out, 'JOB=%s' %job_card

    def xparm_xds(self, real_space_a, real_space_b, real_space_c,
                  space_group, out=None):
        from cctbx import uctbx
        R = self.imagecif_to_xds_transformation_matrix
        unit_cell_a_axis = R * matrix.col(real_space_a)
        unit_cell_b_axis = R * matrix.col(real_space_b)
        unit_cell_c_axis = R * matrix.col(real_space_c)
        A_inv = matrix.sqr(unit_cell_a_axis.elems +
                           unit_cell_b_axis.elems +
                           unit_cell_c_axis.elems)
        metrical_matrix = (A_inv * A_inv.transpose()).as_sym_mat3()
        unit_cell = uctbx.unit_cell(metrical_matrix=metrical_matrix)
        from iotbx.xds import xparm
        writer = xparm.writer(
            self.starting_frame,
            self.starting_angle,
            self.oscillation_range,
            self.rotation_axis,
            self.wavelength,
            self.beam_vector,
            space_group,
            unit_cell.parameters(),
            unit_cell_a_axis.elems,
            unit_cell_b_axis.elems,
            unit_cell_c_axis.elems,
            None, # num_segments
            self.detector_size,
            self.pixel_size,
            self.detector_origin,
            self.detector_distance,
            self.detector_x_axis,
            self.detector_y_axis,
            self.detector_normal,
            segments=None,
            orientation=None)
        writer.show(out=out)
