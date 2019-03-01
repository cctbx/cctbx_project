from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatXTC import FormatXTC, locator_str
from libtbx.phil import parse

try:
    from xfel.cxi.cspad_ana import cspad_tbx
    from xfel.cftbx.detector import cspad_cbf_tbx
except ImportError:
    # xfel not configured
    pass

cspad_locator_str = """
  cspad {
    detz_offset = None
      .type = float
      .help = Distance from back of detector rail to sample interaction region (CXI) \
              or actual detector distance (XPP/MFX)
    apply_gain_mask = True
      .type = bool
      .help = flag to indicate if gain should be applied to cspad data
    dark_correction = True
      .type = bool
      .help = flag to decide if dark correction should be done
    use_psana_calib = False
      .type = bool
      .help = Use the psana calibration
    common_mode = default
      .type = str
      .help = Common mode correction "default", "cspad_default", or "unbonded"\
              see https://confluence.slac.stanford.edu/display/PSDM/Common+mode+correction+algorithms\
              default means no common mode corrections... the other two are psana corrections
    }
"""

cspad_locator_scope = parse(cspad_locator_str + locator_str, process_includes=True)


class FormatXTCCspad(FormatXTC):
    def __init__(self, image_file, locator_scope=cspad_locator_scope, **kwargs):
        assert self.understand(image_file)
        FormatXTC.__init__(self, image_file, locator_scope=locator_scope, **kwargs)
        assert (
            self.params.cspad.detz_offset is not None
        ), "Supply a detz_offset for the cspad"
        self._ds = FormatXTC._get_datasource(image_file, self.params)
        self._psana_runs = FormatXTC._get_psana_runs(self._ds)
        self._cache_psana_det()  # NOTE: move to base FormatXTC class
        self._cache_psana_pedestals()  # NOTE: move to base FormatXTC class
        self._cache_psana_gain()
        self.populate_events()
        self.n_images = len(self.times)

    @staticmethod
    def understand(image_file):
        try:
            params = FormatXTC.params_from_phil(cspad_locator_scope, image_file)
        except Exception:
            return False
        ds = FormatXTC._get_datasource(image_file, params)
        return any(["cspad" in src.lower() for src in params.detector_address])

    def _cache_psana_gain(self):
        """
    checks if user wants gain applied and caches a gain map per run
    """
        run_numbers = self._psana_runs.keys()
        self._gain_masks = {}
        for r in run_numbers:
            if self.params.cspad.apply_gain_mask:
                self._gain_masks[r] = self._psana_det[r].gain_mask(r) > 0
            else:
                self._gain_masks[r] = None

    def _cache_psana_det(self):
        """Store a psana detector instance for each run"""
        assert len(self.params.detector_address) == 1
        import psana

        self._psana_det = {}
        for run_number, run in self._psana_runs.items():
            self._psana_det[run_number] = psana.Detector(
                self.params.detector_address[0], run.env()
            )

    def _cache_psana_pedestals(self):
        """Store a pedestal for each psana detector instance"""
        self._pedestals = {}
        for run_number, run in self._psana_runs.iteritems():
            det = self._psana_det[run_number]
            self._pedestals[run_number] = det.pedestals(run)

    def get_raw_data(self, index):
        from scitbx.array_family import flex
        import numpy as np

        assert len(self.params.detector_address) == 1
        d = self.get_detector(index)
        event = self._get_event(index)
        run_number = event.run()
        det = self._psana_det[run_number]
        data = cspad_cbf_tbx.get_psana_corrected_data(
            det,
            event,
            use_default=self.params.cspad.use_psana_calib,
            dark=self._pedestals[run_number],
            common_mode=self.params.cspad.common_mode,
            apply_gain_mask=self.params.cspad.apply_gain_mask,
            gain_mask_value=None,
            per_pixel_gain=False,
            gain_mask=self._gain_masks[run_number],
        )
        data = data.astype(np.float64)
        self._raw_data = []
        for quad_count, quad in enumerate(d.hierarchy()):
            for sensor_count, sensor in enumerate(quad):
                for asic_count, asic in enumerate(sensor):
                    fdim, sdim = asic.get_image_size()
                    asic_data = data[
                        sensor_count + quad_count * 8,
                        :,
                        asic_count * fdim : (asic_count + 1) * fdim,
                    ]  # 8 sensors per quad
                    self._raw_data.append(flex.double(np.array(asic_data)))
        assert len(d) == len(self._raw_data)
        return tuple(self._raw_data)

    def get_num_images(self):
        return self.n_images

    def get_detector(self, index=None):
        return FormatXTCCspad._detector(self, index)

    def get_beam(self, index=None):
        return self._beam(index)

    def _beam(self, index=None):
        """Returns a simple model for the beam """
        if index is None:
            index = 0
        evt = self._get_event(index)
        wavelength = cspad_tbx.evt_wavelength(evt)
        if wavelength is None:
            return None
        return self._beam_factory.simple(wavelength)

    def get_goniometer(self, index=None):
        return None

    def get_scan(self, index=None):
        return None

    # XXX Implement recursive version
    def _detector(self, index=None):
        import psana
        from xfel.cftbx.detector.cspad_cbf_tbx import read_slac_metrology
        from dxtbx.model import Detector
        from scitbx.matrix import col
        from dxtbx.model import ParallaxCorrectedPxMmStrategy
        from xfel.cxi.cspad_ana.cspad_tbx import env_distance

        if index is None:
            index = 0

        ev = self._get_event(index)
        run_number = ev.run()
        run = self._psana_runs[run_number]
        det = self._psana_det[run_number]
        geom = det.pyda.geoaccess(run_number)
        cob = read_slac_metrology(geometry=geom, include_asic_offset=True)
        distance = env_distance(
            self.params.detector_address[0], run.env(), self.params.cspad.detz_offset
        )
        d = Detector()
        pg0 = d.hierarchy()
        # first deal with D0
        det_num = 0
        origin = col((cob[(0,)] * col((0, 0, 0, 1)))[0:3])
        fast = col((cob[(0,)] * col((1, 0, 0, 1)))[0:3]) - origin
        slow = col((cob[(0,)] * col((0, 1, 0, 1)))[0:3]) - origin
        origin += col((0.0, 0.0, -distance))
        pg0.set_local_frame(fast.elems, slow.elems, origin.elems)
        pg0.set_name("D%d" % (det_num))
        for quad_num in xrange(4):
            # Now deal with Qx
            pg1 = pg0.add_group()
            origin = col((cob[(0, quad_num)] * col((0, 0, 0, 1)))[0:3])
            fast = col((cob[(0, quad_num)] * col((1, 0, 0, 1)))[0:3]) - origin
            slow = col((cob[(0, quad_num)] * col((0, 1, 0, 1)))[0:3]) - origin
            pg1.set_local_frame(fast.elems, slow.elems, origin.elems)
            pg1.set_name("D%dQ%d" % (det_num, quad_num))
            for sensor_num in xrange(8):
                # Now deal with Sy
                pg2 = pg1.add_group()
                origin = col((cob[(0, quad_num, sensor_num)] * col((0, 0, 0, 1)))[0:3])
                fast = (
                    col((cob[(0, quad_num, sensor_num)] * col((1, 0, 0, 1)))[0:3])
                    - origin
                )
                slow = (
                    col((cob[(0, quad_num, sensor_num)] * col((0, 1, 0, 1)))[0:3])
                    - origin
                )
                pg2.set_local_frame(fast.elems, slow.elems, origin.elems)
                pg2.set_name("D%dQ%dS%d" % (det_num, quad_num, sensor_num))
                # Now deal with Az
                for asic_num in xrange(2):
                    val = "ARRAY_D0Q%dS%dA%d" % (quad_num, sensor_num, asic_num)
                    p = pg2.add_panel()
                    origin = col(
                        (cob[(0, quad_num, sensor_num, asic_num)] * col((0, 0, 0, 1)))[
                            0:3
                        ]
                    )
                    fast = (
                        col(
                            (
                                cob[(0, quad_num, sensor_num, asic_num)]
                                * col((1, 0, 0, 1))
                            )[0:3]
                        )
                        - origin
                    )
                    slow = (
                        col(
                            (
                                cob[(0, quad_num, sensor_num, asic_num)]
                                * col((0, 1, 0, 1))
                            )[0:3]
                        )
                        - origin
                    )
                    p.set_local_frame(fast.elems, slow.elems, origin.elems)
                    p.set_pixel_size(
                        (cspad_cbf_tbx.pixel_size, cspad_cbf_tbx.pixel_size)
                    )
                    p.set_image_size(cspad_cbf_tbx.asic_dimension)
                    p.set_trusted_range(
                        (
                            cspad_tbx.cspad_min_trusted_value,
                            cspad_tbx.cspad_saturated_value,
                        )
                    )
                    p.set_name(val)

        try:
            beam = self._beam(index)
        except Exception:
            print(
                "No beam object initialized. Returning CSPAD detector without parallax corrections"
            )
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


if __name__ == "__main__":
    import sys

    for arg in sys.argv[1:]:
        print(FormatXTCCspad.understand(arg))
