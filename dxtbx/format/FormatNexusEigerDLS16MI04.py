from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from dxtbx.format.FormatNexus import FormatNexus


class FormatNexusEigerDLS16MI04(FormatNexus):
    @staticmethod
    def understand(image_file):
        """Check to see if this format class can understand the image file.

    Args:
      image_file (str): The file path of the image file to check.

    Returns:
      bool: Returns ``True`` if the image_file is understood by this format class,
      else returns ``False``.

    """

        # this depends on DIALS for the goniometer shadow model; if missing
        # simply return False

        try:
            from dials.util.masking import GoniometerShadowMaskGenerator  # test import
        except ImportError:
            return False

        import h5py

        # Get the file handle
        handle = h5py.File(image_file, "r")
        if (
            "short_name" not in handle["/entry/instrument"].attrs
            or handle["/entry/instrument"].attrs["short_name"].lower() != "i04"
        ):
            return False

        return FormatNexus.understand(image_file)

    @staticmethod
    def has_dynamic_shadowing(**kwargs):
        import libtbx

        dynamic_shadowing = kwargs.get("dynamic_shadowing", False)
        if dynamic_shadowing in (libtbx.Auto, "Auto"):
            return True
        return dynamic_shadowing

    def __init__(self, image_file, **kwargs):
        """Initialise the image structure from the given file."""

        import libtbx
        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)

        self._dynamic_shadowing = self.has_dynamic_shadowing(**kwargs)
        super(FormatNexusEigerDLS16MI04, self).__init__(image_file, **kwargs)

    def get_mask(self, index, goniometer=None):
        mask = super(FormatNexusEigerDLS16MI04, self).get_mask()
        if mask is None:
            # XXX mask really shouldn't be None
            # https://jira.diamond.ac.uk/browse/SCI-8308
            mask = tuple(
                flex.bool(flex.grid(reversed(panel.get_image_size())), True)
                for panel in self.get_detector()
            )
            panel = self.get_detector()[0]
        if self._dynamic_shadowing:
            gonio_masker = self.get_goniometer_shadow_masker(goniometer=goniometer)
            scan = self.get_scan()
            detector = self.get_detector()
            shadow_mask = gonio_masker.get_mask(
                detector, scan.get_angle_from_image_index(index)
            )
            assert len(mask) == len(shadow_mask)
            for m, sm in zip(mask, shadow_mask):
                if sm is not None:
                    m &= sm
        return mask

    def get_goniometer_shadow_masker(self, goniometer=None):
        if goniometer is None:
            goniometer = self.get_goniometer()

        assert goniometer is not None

        if goniometer.get_names()[1] == "chi":
            # SmarGon
            from dxtbx.format.SmarGonShadowMask import SmarGonShadowMaskGenerator

            return SmarGonShadowMaskGenerator(goniometer)

        else:
            raise RuntimeError(
                "Don't understand this goniometer: %s" % list(goniometer.get_names())
            )


if __name__ == "__main__":

    import sys

    for arg in sys.argv[1:]:
        print(FormatNexusEigerDLS16MI04.understand(arg))
