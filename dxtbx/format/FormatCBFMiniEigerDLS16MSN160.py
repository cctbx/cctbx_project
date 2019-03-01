from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from dxtbx.format.FormatCBFMiniEiger import FormatCBFMiniEiger


class FormatCBFMiniEigerDLS16MSN160(FormatCBFMiniEiger):
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

        header = FormatCBFMiniEiger.get_cbf_header(image_file)

        for record in header.split("\n"):
            if (
                "# detector" in record.lower()
                and "eiger" in record.lower()
                and "S/N 160-0001" in header
            ):
                return True

        return False

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
        super(FormatCBFMiniEigerDLS16MSN160, self).__init__(image_file, **kwargs)

    def _goniometer(self):
        """Return a model for a multi-axis goniometer.

    This should be checked against the image header, though for miniCBF
    there are limited options for this.

    Returns:
      dxtbx.model.Goniometer.MultiAxisGoniometer: The goniometer model for
      this detector.

    """

        if "Phi" in self._cif_header_dictionary:
            phi = float(self._cif_header_dictionary["Phi"].split()[0])
        else:
            phi = 0.0

        if "Chi" in self._cif_header_dictionary:
            chi = float(self._cif_header_dictionary["Chi"].split()[0])
        else:
            chi = 0.0

        if "Omega" in self._cif_header_dictionary:
            omega = float(self._cif_header_dictionary["Omega"].split()[0])
        else:
            omega = 0.0

        phi_axis = (1, 0, 0)
        chi_axis = (0, 0, -1)
        omega_axis = (1, 0, 0)
        axes = flex.vec3_double((phi_axis, chi_axis, omega_axis))
        angles = flex.double((phi, chi, omega))
        names = flex.std_string(("GON_PHI", "GON_CHI", "GON_OMEGA"))
        return self._goniometer_factory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    def get_mask(self, goniometer=None):
        mask = super(FormatCBFMiniEigerDLS16MSN160, self).get_mask()
        if self._dynamic_shadowing:
            gonio_masker = self.get_goniometer_shadow_masker(goniometer=goniometer)
            scan = self.get_scan()
            detector = self.get_detector()
            shadow_mask = gonio_masker.get_mask(detector, scan.get_oscillation()[0])
            assert len(mask) == len(shadow_mask)
            for m, sm in zip(mask, shadow_mask):
                if sm is not None:
                    m &= sm
        return mask

    def get_goniometer_shadow_masker(self, goniometer=None):
        if goniometer is None:
            goniometer = self.get_goniometer()

        assert goniometer is not None

        if goniometer.get_names()[1] == "GON_CHI":
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
        print(FormatCBFMiniEigerDLS16MSN160.understand(arg))
