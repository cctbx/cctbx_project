from __future__ import division

class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_beam()
        self.tst_detector()
        self.tst_goniometer()
        self.tst_scan()

    def tst_beam(self):
        from dxtbx.serialize import beam
        from dxtbx.model import Beam
        b1 = Beam((1, 0, 0), 2, 0.1, 0.1)
        d = beam.to_dict(b1)
        b2 = beam.from_dict(d)
        assert(d['direction'] == (1, 0, 0))
        assert(d['wavelength'] == 2)
        assert(abs(d['divergence'] - 0.1) <= 1e-7)
        assert(abs(d['sigma_divergence'] - 0.1) <= 1e-7)
        assert(b1 == b2)
        print 'OK'

    def tst_detector(self):
        self.tst_panel()
        print 'OK'

    def tst_panel(self):
        from dxtbx.serialize.detector import panel_to_dict, panel_from_dict
        from dxtbx.model import Panel
        p1 = Panel("UNKNOWN",
                   (1, 0, 0), (0, 1, 0), (0, 0, 1),
                   (0.072, 0.072),
                   (1000, 2000), (0, 100))
        d = panel_to_dict(p1)
        p2 = panel_from_dict(d)
        assert(d['type'] == "UNKNOWN")
        assert(d['fast_axis'] == (1, 0, 0))
        assert(d['slow_axis'] == (0, 1, 0))
        assert(d['origin'] == (0, 0, 1))
        assert(d['pixel_size'] == (0.072, 0.072))
        assert(d['image_size'] == (1000, 2000))
        assert(d['trusted_range'] == (0, 100))
        assert(p1 == p2)
        print 'OK'

    def tst_goniometer(self):
        from dxtbx.serialize import goniometer
        from dxtbx.model import Goniometer
        g1 = Goniometer((1, 0, 0), (1, 0, 0, 0, 1, 0, 0, 0, 1))
        d = goniometer.to_dict(g1)
        g2 = goniometer.from_dict(d)
        assert(d['rotation_axis'] == (1, 0, 0))
        assert(d['fixed_rotation'] == (1, 0, 0, 0, 1, 0, 0, 0, 1))
        assert(g1 == g2)
        print 'OK'

    def tst_scan(self):
        from dxtbx.serialize import scan
        from dxtbx.model import Scan
        from scitbx.array_family import flex
        s1 = Scan((1, 3), (1.0, 0.2), 0.1, flex.double([0.1, 0.2, 0.3]))
        d = scan.to_dict(s1)
        s2 = scan.from_dict(d)
        assert(d['image_range'] == (1, 3))
        assert(d['oscillation'] == (1.0, 0.2))
        assert(d['exposure_time'] == 0.1)
        assert(d['epochs'] == [0.1, 0.2, 0.3])
        assert(s1 == s2)
        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()
