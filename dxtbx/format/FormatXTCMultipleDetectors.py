from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatXTC import FormatXTC, locator_str
from dxtbx.format.FormatXTCRayonix import FormatXTCRayonix, rayonix_locator_str
from dxtbx.format.FormatXTCCspad import FormatXTCCspad, cspad_locator_str
from dxtbx.format.FormatXTCJungfrau import FormatXTCJungfrau, jungfrau_locator_str
from libtbx.phil import parse

multiple_locator_scope = parse(
    rayonix_locator_str + cspad_locator_str + jungfrau_locator_str + locator_str,
    process_includes=True,
)


class FormatXTCMultipleDetectors(FormatXTCRayonix, FormatXTCCspad, FormatXTCJungfrau):
    def __init__(self, image_file, **kwargs):
        assert self.understand(image_file)
        FormatXTC.__init__(
            self, image_file, locator_scope=multiple_locator_scope, **kwargs
        )
        self._ds = FormatXTCMultipleDetectors._get_datasource(image_file, self.params)
        self._env = self._ds.env()
        self.populate_events()
        self.n_images = len(self.times)

        if any(["rayonix" in src.lower() for src in self.params.detector_address]):
            FormatXTCRayonix.__init__(self, image_file, **kwargs)
        if any(["cspad" in src.lower() for src in self.params.detector_address]):
            FormatXTCCspad.__init__(self, image_file, **kwargs)
        if any(["jungfrau" in src.lower() for src in self.params.detector_address]):
            FormatXTCJungfrau.__init__(self, image_file, **kwargs)
        FormatXTC.__init__(
            self, image_file, locator_scope=multiple_locator_scope, **kwargs
        )

    @staticmethod
    def understand(image_file):
        try:
            params = FormatXTC.params_from_phil(multiple_locator_scope, image_file)
        except Exception:
            return False
        ds = FormatXTC._get_datasource(image_file, params)

        if params.detector_address is None or len(params.detector_address) <= 1:
            return False
        return [
            "rayonix" in src.lower()
            or "cspad" in src.lower()
            or "jungfrau" in src.lower()
            for src in params.detector_address
        ].count(True) >= 2

    def get_raw_data(self, index=None):
        if index is None:
            index = 0

        all_addresses = self.params.detector_address

        raw_data = []

        for address in all_addresses:
            self.params.detector_address = [address]
            sub_d = None
            if "rayonix" in address.lower():
                data = FormatXTCRayonix.get_raw_data(self, index)
            elif "cspad" in address.lower():
                data = FormatXTCCspad.get_raw_data(self, index)
            elif "jungfrau" in address.lower():
                data = FormatXTCJungfrau.get_raw_data(self, index)
            assert data is not None, address
            if not isinstance(data, tuple):
                data = (data,)
            raw_data.extend(data)
        self.params.detector_address = all_addresses

        self._raw_data = raw_data
        return tuple(raw_data)

    def get_detector(self, index=None):
        return self._detector(index)

    def _detector(self, index=None):
        from dxtbx.model import Detector

        if index is None:
            index = 0
        d = Detector()
        root = d.hierarchy()

        all_addresses = self.params.detector_address

        def recursive_add_node(a, b):
            # add a to b
            if a.is_panel():
                b.add_panel(a)
            else:
                g = b.add_group(a)
                for child in a:
                    recursive_add_node(child, g)

        for address in all_addresses:
            self.params.detector_address = [address]
            sub_d = None
            try:
                if "rayonix" in address.lower():
                    sub_d = FormatXTCRayonix._detector(self)
                elif "cspad" in address.lower():
                    sub_d = FormatXTCCspad._detector(self, index)
                elif "jungfrau" in address.lower():
                    sub_d = FormatXTCJungfrau._detector(self, index)
            except Exception:
                continue
            assert sub_d is not None, address
            recursive_add_node(sub_d.hierarchy(), root)
        self.params.detector_address = all_addresses

        return d


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatXTCMultipleDetectors.understand(arg))
