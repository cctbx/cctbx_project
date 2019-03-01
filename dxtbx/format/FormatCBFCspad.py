#!/usr/bin/env python
# FormatCBFCspad.py
#
# Contains class methods specific to interacting with CSPAD images
#
# $Id:
#

from __future__ import absolute_import, division

from dxtbx.format.FormatCBFMultiTileHierarchy import FormatCBFMultiTileHierarchyStill
from dxtbx.format.FormatCBFFull import FormatCBFFullStill
from dxtbx.model import ParallaxCorrectedPxMmStrategy
from scitbx.matrix import col, sqr
import pycbf


class FormatCBFCspad(FormatCBFMultiTileHierarchyStill):
    """An image reading class CSPAD CBF files"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CSPD CBF format image, i.e. we can
        make sense of it."""

        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_widefile(image_file, pycbf.MSG_DIGEST)

        cbf_handle.find_category("diffrn_detector")
        if cbf_handle.count_rows() > 1:
            return False  # support 1 detector per file for now

        cbf_handle.find_column("type")

        return cbf_handle.get_value() == "CS PAD"

    def _detector(self):
        d = FormatCBFMultiTileHierarchyStill._detector(self)

        try:
            # a header only CBF file will not have a beam object
            beam = self._beam()
        except Exception as e:
            if "CBF_NOTFOUND" not in str(e):
                raise e
            return d

        # take into consideration here the thickness of the sensor also the
        # wavelength of the radiation (which we have in the same file...)
        wavelength = beam.get_wavelength()
        thickness = 0.5  # mm, see Hart et al. 2012
        from cctbx.eltbx import attenuation_coefficient

        table = attenuation_coefficient.get_table("Si")
        # mu_at_angstrom returns cm^-1
        mu = table.mu_at_angstrom(wavelength) / 10.0  # mu: mm^-1
        t0 = thickness

        for panel in d:
            panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
        return d

    def sync_detector_to_cbf(self, detector=None):
        """If the dectector object has been changed, due to refinment or manual repositioning
        in a gui, call this function to synchronize these changes to the underlying CBF handle"""

        def recursive_sync(cbf, group, cbf_detectors=None, root=False):
            """ Walks the hierarchy and synchronizes the dxtbx matrices to the CBF axes."""

            d_mat = group.get_local_d_matrix()
            fast = col((d_mat[0], d_mat[3], d_mat[6])).normalize()
            slow = col((d_mat[1], d_mat[4], d_mat[7])).normalize()
            orig = col((d_mat[2], d_mat[5], d_mat[8]))

            if group.is_panel():
                if cbf.has_sections():
                    # use the pre-mapping
                    cbf_detector = cbf_detectors[group.get_name()]
                else:
                    # figure out which panel number this panel is by finding it in diffrn_data_frame
                    cbf.find_category("diffrn_data_frame")
                    cbf.find_column("array_id")
                    cbf.find_row(group.get_name())
                    cbf.find_column("binary_id")
                    array_num = int(cbf.get_value()) - 1

                    cbf_detector = cbf.construct_detector(array_num)

                axis0 = cbf_detector.get_detector_surface_axes(0)
                axis1 = cbf_detector.get_detector_surface_axes(1)
                assert cbf.get_axis_depends_on(axis0) == axis1
                name = cbf.get_axis_depends_on(axis1)

                # account for the sign change when moving from dxtbx to ImageCIF
                slow = col((slow[0], -slow[1], slow[2]))

                # account for the offset from the center of the sensor to the center
                # of the ASIC.
                orig -= col(cbf.get_axis_offset(axis1))
            else:
                name = group.get_name()

            v3 = fast.cross(slow).normalize()

            if root:
                # use the X, Y, Z settings table to record root detector position instead of axis offsets
                dx = orig[0]
                dy = orig[1]
                distance = orig[2]
                orig = col((0, 0, 0))

            cbf.find_category("axis")
            cbf.find_column("id")
            while True:
                if cbf.get_value() == name:
                    axis_name = cbf.get_value()
                    break
                cbf.next_row()

            cbf.find_column("offset[1]")
            cbf.set_value(str(orig[0]))
            cbf.find_column("offset[2]")
            cbf.set_value(str(orig[1]))
            cbf.find_column("offset[3]")
            cbf.set_value(str(orig[2]))

            r3 = sqr(
                (
                    fast[0],
                    slow[0],
                    v3[0],
                    fast[1],
                    slow[1],
                    v3[1],
                    fast[2],
                    slow[2],
                    v3[2],
                )
            )

            angle, axis = r3.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(
                deg=True
            )
            if angle == 0:
                axis = col((0, 0, 1))

            assert axis.length() > 0

            cbf.find_column("vector[1]")
            cbf.set_value(str(axis[0]))
            cbf.find_column("vector[2]")
            cbf.set_value(str(axis[1]))
            cbf.find_column("vector[3]")
            cbf.set_value(str(axis[2]))

            cbf.set_axis_setting(axis_name, angle, 0)

            if root:
                # synchronize the new root origin
                for axis_id, setting_value in zip(
                    ["_X", "_Y", "_Z"], [dx, dy, distance]
                ):
                    while not axis_id in axis_name:
                        axis_name = cbf.get_axis_depends_on(axis_name)
                    assert cbf.get_axis_type(axis_name) == "translation"
                    cbf.set_axis_setting(axis_name, setting_value, 0)

            if group.is_group():
                for subgroup in group:
                    recursive_sync(cbf, subgroup, cbf_detectors)

        if detector is None:
            detector = self.get_detector()

        root = detector.hierarchy()
        assert len(root) == 4

        cbf = self._get_cbf_handle()
        if cbf.has_sections():
            """
            when using sections, the order of the panels in the detector object doesn't match the
            order returned by cbflib's construct_detector.  The mapping between them is found
            in the table array_structure list, which maps an array_section_id to axis_set_id,
            which can then be matchedup with the axes of the object returned by construct_detector
            """
            all_cbfdetectors = [
                cbf.construct_detector(i) for i in xrange(len(detector))
            ]
            all_panelnames = [panel.get_name() for panel in detector]

            # map the array_section_ids, which match the panel names, to their root axis names
            panel_name_mapping = {}
            cbf.find_category("array_structure_list")
            for i in xrange(cbf.count_rows()):
                cbf.find_column("array_section_id")
                name = cbf.get_value()
                if name in all_panelnames and name not in panel_name_mapping.values():
                    cbf.find_column("axis_set_id")
                    panel_name_mapping[cbf.get_value()] = name
                cbf.next_row()

            # map the panel names (same as array_section_ids) to specific cbf detector objects
            mapped_detectors = {}
            for cbf_d in all_cbfdetectors:
                root_axis = cbf_d.get_detector_surface_axes(0)
                mapped_detectors[panel_name_mapping[root_axis]] = cbf_d
        else:
            mapped_detectors = None

        recursive_sync(cbf, root, mapped_detectors, True)


class FormatCBFFullStillInMemory(FormatCBFFullStill):
    """Overrides the Format object's init method to accept a cbf handle instead
    of a file name. Used with XFELs when it is desirable to never write
    a file to disk, but to process it only in memory.
    """

    @staticmethod
    def understand(image_file):
        """This class shouldn't be found by the dxtbx Registry. Instead it should
        be instantiated directly as needed
        """
        return False

    def __init__(self, cbf_handle, **kwargs):
        """ @param cbf_handle In memory cbf_handle, alredy initialized """
        from dxtbx.model.detector import Detector, DetectorFactory
        from dxtbx.model.beam import Beam, BeamFactory

        self._goniometer_instance = None
        self._scan_instance = None
        self._detector_instance = None
        self._detector_factory = DetectorFactory
        self._beam_factory = BeamFactory

        self._cbf_handle = cbf_handle
        self._raw_data = None
        try:
            detector_instance = self._detector()
            assert isinstance(detector_instance, Detector)
            self._detector_instance = detector_instance

            beam_instance = self._beam()
            assert isinstance(beam_instance, Beam)
            self._beam_instance = beam_instance

        except Exception:
            # FIXME ideally should not squash the errors here...
            pass
        finally:
            self._end()

    def __getstate__(self):
        # For pickling and copying, don't carry the cbf_handle
        return self._detector_instance, self._beam_instance, self.get_raw_data()

    def __setstate__(self, state):
        # For pickling and copying, don't carry the cbf_handle
        self._goniometer_instance = None
        self._scan_instance = None
        self._detector_instance = state[0]
        self._beam_instance = state[1]
        self._raw_data = state[2]

    def __del__(self):
        # If copied or pickled, may not have _cbf_handle any more
        if hasattr(self, "_cbf_handle") and self._cbf_handle is not None:
            self._cbf_handle.__swig_destroy__(self._cbf_handle)


class FormatCBFMultiTileHierarchyStillInMemory(
    FormatCBFFullStillInMemory, FormatCBFMultiTileHierarchyStill
):
    def get_raw_data(self):
        # need to specify the get_raw_data function needed
        return FormatCBFMultiTileHierarchyStill.get_raw_data(self)


class FormatCBFCspadInMemory(FormatCBFMultiTileHierarchyStillInMemory, FormatCBFCspad):
    """ mixin class for in memory cspad images"""
