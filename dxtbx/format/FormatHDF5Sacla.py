from __future__ import absolute_import, division
from dxtbx.format.Format import Format
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.format.FormatStill import FormatStill


class FormatHDF5Sacla(FormatHDF5, FormatStill):
    """
  Class for reading SACLA images created by the DataConvert SACLA
  script (this script lives on the SACLA hpc).

  This assumes the argument -reconstr was passed to
  DataConvert in order to Reconstruct the image.

  Also, this processes only a single run's worth of data in the hdf5
  """

    @staticmethod
    def understand(image_file):
        import h5py

        h5_handle = h5py.File(image_file, "r")
        understood = False
        if "file_info" in h5_handle and "run_number_list" in h5_handle["file_info"]:
            run_grp = FormatHDF5Sacla._get_run_h5group(h5_handle)
            if "detector_2d_assembled_1" in list(run_grp.keys()):
                understood = True
        return understood

    def __init__(self, image_file, **kwargs):
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatHDF5.__init__(self, image_file, **kwargs)

    def _start(self):
        import h5py

        self._h5_handle = h5py.File(self.get_image_file(), "r")
        self._run = FormatHDF5Sacla._get_run_h5group(self._h5_handle)
        event_info = self._run["event_info"]
        tag_number_list = event_info["tag_number_list"]
        self._images = ["tag_%d" % tag for tag in tag_number_list]

    @staticmethod
    def _get_run_h5group(h5_handle):
        """returns the first run group found"""
        file_info = h5_handle["file_info"]
        run_number_list = file_info["run_number_list"]
        run_str = "run_%d" % run_number_list[0]
        return h5_handle[run_str]

    def _detector(self, index=None):
        from scitbx import matrix

        # Get the pixel and image size
        detector_2d_assembled_1 = self._run["detector_2d_assembled_1"]
        detector_info = detector_2d_assembled_1["detector_info"]
        pixel_size = (
            detector_info["pixel_size_in_micro_meter"][0] / 1000,
            detector_info["pixel_size_in_micro_meter"][1] / 1000,
        )
        tag = detector_2d_assembled_1[self._images[0]]
        data = tag["detector_data"][()]

        # detector_image_size is fast-, slow- , just in case the dataset is ever non-square
        image_size = (data.shape[1], data.shape[0])
        trusted_range = (0, 200000)

        # Initialise detector frame
        fast = matrix.col((1.0, 0.0, 0.0))
        slow = matrix.col((0.0, -1.0, 0.0))
        orig = matrix.col(
            (
                -image_size[0] * pixel_size[0] / 2,
                image_size[1] * pixel_size[1] / 2,
                -100.0,
            )
        )
        # Make the detector
        return self._detector_factory.make_detector(
            "", fast, slow, orig, pixel_size, image_size, trusted_range
        )

    def _beam(self, index=None):
        run_info = self._run["run_info"]
        sacla_config = run_info["sacla_config"]
        eV = sacla_config["photon_energy_in_eV"].value

        return self._beam_factory.simple(12398.4 / eV)

    def get_num_images(self):
        return len(self._images)

    def get_raw_data(self, index=0):
        from scitbx.array_family import flex
        import numpy as np

        detector_2d_assembled_1 = self._run["detector_2d_assembled_1"]
        tag = detector_2d_assembled_1[self._images[index]]
        return flex.double(tag["detector_data"].value.astype(np.float64))

    def get_detectorbase(self, index=None):
        raise NotImplementedError

    def get_image_file(self, index=None):
        return Format.get_image_file(self)

    def get_detector(self, index=None):
        return self._detector_instance

    def get_beam(self, index=None):
        return self._beam_instance


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatHDF5Sacla.understand(arg))
