# Reads a multi-tile CBF image, discovering its detector geometry
# automatically

from __future__ import absolute_import, division, print_function

import pycbf

from dxtbx.format.FormatCBFFull import FormatCBFFull
from dxtbx.format.FormatStill import FormatStill
from dxtbx.model.detector import Detector


class cbf_wrapper(pycbf.cbf_handle_struct):
    """ Wrapper class that provids convience functions for working with cbflib"""

    def add_category(self, name, columns):
        """ Create a new category and populate it with column names """
        self.new_category(name)
        for column in columns:
            self.new_column(column)

    def add_row(self, data):
        """Add a row to the current category.  If data contains more entries than
        there are columns in this category, then the remainder is truncated
        Use '.' for an empty value in a row."""
        self.new_row()
        self.rewind_column()
        for item in data:
            self.set_value(item)
            if item == ".":
                self.set_typeofvalue("null")
            try:
                self.next_column()
            except Exception:
                break

    def has_sections(self):
        """True if the cbf has the array_structure_list_section table, which
        changes how its data is stored in the binary sections
        """
        try:
            self.find_category("array_structure_list_section")
        except Exception as e:
            if "CBF_NOTFOUND" not in str(e):
                raise e
            return False
        return True


class FormatCBFMultiTile(FormatCBFFull):
    """An image reading class multi-tile CBF files"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        try:
            cbf_handle = pycbf.cbf_handle_struct()
            cbf_handle.read_widefile(image_file, pycbf.MSG_DIGEST)
        except Exception as e:
            if "CBFlib Error" in str(e):
                return False

        # check if multiple arrays
        try:
            return cbf_handle.count_elements() > 1
        except Exception as e:
            if "CBFlib Error" in str(e):
                return False

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        FormatCBFFull.__init__(self, image_file, **kwargs)

        self._raw_data = None

    def _start(self):
        """Open the image file as a cbf file handle, and keep this somewhere
        safe."""
        from dxtbx.format.FormatCBF import FormatCBF

        FormatCBF._start(self)  # Note, skip up an inhieritance level

    def detectorbase_start(self):
        pass

    def _get_cbf_handle(self):
        try:
            return self._cbf_handle
        except AttributeError:
            self._cbf_handle = cbf_wrapper()
            self._cbf_handle.read_widefile(self._image_file, pycbf.MSG_DIGEST)
            return self._cbf_handle

    def _detector(self):
        """Return a working detector instance."""

        cbf = self._get_cbf_handle()

        d = Detector()

        for i in xrange(cbf.count_elements()):
            ele_id = cbf.get_element_id(i)
            cbf.find_category("diffrn_data_frame")
            cbf.find_column("detector_element_id")
            cbf.find_row(ele_id)
            cbf.find_column("array_id")
            array_id = cbf.get_value()

            cbf_detector = cbf.construct_detector(i)

            p = d.add_panel()
            p.set_name(array_id)

            # code adapted below from dxtbx.model.detector.DetectorFactory.imgCIF_H
            pixel = (
                cbf_detector.get_inferred_pixel_size(1),
                cbf_detector.get_inferred_pixel_size(2),
            )

            fast = cbf_detector.get_detector_axes()[0:3]
            slow = cbf_detector.get_detector_axes()[3:6]
            origin = cbf_detector.get_pixel_coordinates_fs(0, 0)

            size = tuple(reversed(cbf.get_image_size(0)))

            try:
                cbf.find_category("array_intensities")
                cbf.find_column("undefined_value")
                underload = cbf.get_doublevalue()
                overload = cbf.get_overload(0)
                trusted_range = (underload, overload)
            except Exception:
                trusted_range = (0.0, 0.0)

            cbf_detector.__swig_destroy__(cbf_detector)
            del cbf_detector

            p.set_local_frame(fast, slow, origin)

            p.set_pixel_size(tuple(map(float, pixel)))
            p.set_image_size(size)
            p.set_trusted_range(tuple(map(float, trusted_range)))
            # p.set_px_mm_strategy(px_mm) FIXME

        return d

    def _beam(self):
        """Return a working beam instance."""

        return self._beam_factory.imgCIF_H(self._get_cbf_handle())

    def get_raw_data(self):
        if self._raw_data is None:
            self._raw_data = []

            cbf = self._get_cbf_handle()

            # find the data
            cbf.select_category(0)
            while cbf.category_name().lower() != "array_data":
                try:
                    cbf.next_category()
                except Exception:
                    return None
            cbf.select_column(0)
            cbf.select_row(0)

            d = self.get_detector()

            import numpy
            from scitbx.array_family import flex

            for panel in d:
                name = panel.get_name()
                cbf.find_column("array_id")
                assert name == cbf.get_value()

                cbf.find_column("data")
                assert cbf.get_typeofvalue().find("bnry") > -1

                image_string = cbf.get_realarray_as_string()
                image = flex.double(numpy.fromstring(image_string, numpy.float))

                parameters = cbf.get_realarrayparameters_wdims_fs()
                image_size = (parameters[6], parameters[5])

                image.reshape(flex.grid(*image_size))

                self._raw_data.append(image)

                try:
                    cbf.next_row()
                except Exception:
                    break
            assert len(d) == len(self._raw_data)

        return tuple(self._raw_data)


class FormatCBFMultiTileStill(FormatStill, FormatCBFMultiTile):
    """An image reading class for full CBF format images i.e. those from
    a variety of cameras which support this format. Custom derived from
    the FormatStill to handle images without a gonimeter or scan"""

    @staticmethod
    def understand(image_file):
        """Check to see if this looks like an CBF format image, i.e. we can
        make sense of it."""

        header = FormatCBFMultiTile.get_cbf_header(image_file)

        # According to ImageCIF, "Data items in the DIFFRN_MEASUREMENT_AXIS
        # category associate axes with goniometers."
        # http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Cdiffrn_measurement_axis.html
        if "diffrn_measurement_axis" in header:
            return False
        return True


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatCBFMultiTile.understand(arg))
