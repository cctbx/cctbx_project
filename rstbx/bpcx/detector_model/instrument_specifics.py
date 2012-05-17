# Hardware specific code that makes use of the abstract definitions in
# module detector

from rstbx.bpcx import sensor
from rstbx.bpcx.detector_model.detector import detector
from rstbx.cftbx.coordinate_frame_converter import coordinate_frame_converter
from libtbx.test_utils import approx_equal

class pilatus(detector):
    '''Trivial implementation of a PILATUS type detector of unspecified
    size. Single sensor surface fixed with pixel fast and slow directions
    parallel to the basis directions on that abstract surface. The origin
    is currently chosen in the corner of the detector frame'''

    # For this trivial example we require that there is only a single
    # sensor surface. The idea is that this class advertises its
    # suitability to a factory function.

    @staticmethod
    def understand(cfc):
        '''Ultimately prefer to understand an image file directly, but
        for the Use Case we use coordinate_frame_converter as an agent
        for obtaining information about the experiment. We will also
        want a more nuanced report of understanding. Here it is just
        boolean'''

        #TODO Implementation
        return True

    def __init__(self, origin, dir1, dir2, lim1, lim2, pixel_size_fast,
                 pixel_size_slow, npx_fast, npx_slow):

        self._px_size_fast = pixel_size_fast
        self._px_size_slow = pixel_size_slow
        self._npx_fast = npx_fast
        self._npx_slow = npx_slow

        # set up a single sensor
        s = sensor(origin, dir1, dir2, lim1, lim2)

        # set up the abstract detector
        detector.__init__(self, [s])

    def mm_to_px(self, intersection):
        '''Take an intersection in mm in the abstract coordinate frame
        attached to the sensor and convert it to coordinates in units of
        pixels in image space. This basic implementation assumes that the
        two origins coincide, the fast axis is aligned along dir1 and
        the slow axis along dir2. This does no bounds checking, it is
        purely a units conversion function'''

        fast = intersection[0] / self._px_size_fast
        slow = intersection[1] / self._px_size_slow

        return (fast, slow)

    # getters
    def px_size_fast(self):
        return self._px_size_fast

    def px_size_slow(self):
        return self._px_size_slow

    def npx_fast(self):
        return self._npx_fast

    def npx_slow(self):
        return self._npx_slow


class detector_factory_from_cfc:
    '''Identify and return a suitable detector class based on the
    information in a coordinate_frame_converter object'''

    def __init__(self, cfc):

        self._cfc = cfc
        self._known_detectors = [pilatus]

    def choose(self):

        choice = [d for d in self._known_detectors if d.understand(self._cfc)]
        return choice[0]

    def build(self):

        d = self.choose()
        origin = self._cfc.get_c('detector_origin')
        dir1 = self._cfc.get_c('detector_fast')
        dir2 = self._cfc.get_c('detector_slow')
        px_lim = self._cfc.get('detector_size_fast_slow')
        px_size_fast, px_size_slow = self._cfc.get('detector_pixel_size_fast_slow')

        lim1 = (0, px_lim[0] * px_size_fast)
        lim2 = (0, px_lim[1] * px_size_slow)

        return d(origin, dir1, dir2, lim1, lim2, px_size_fast,
                 px_size_slow, px_lim[0], px_lim[1])


if __name__ == '__main__':

    import sys

    if len(sys.argv) != 2:
        raise RuntimeError, '%s configuration-file' % sys.argv[0]

    configuration_file = sys.argv[1]

    cfc = coordinate_frame_converter(configuration_file)

    df = detector_factory_from_cfc(cfc)
    d = df.build()

    # A little test. Intersection at centre of the abstract frame should
    # be at the centre of the pixel grid
    s = d.sensors()[0]
    intersection = ((s.lim1[1] - s.lim1[0]) / 2., (s.lim2[1] - s.lim2[0]) / 2.)

    assert approx_equal(d.mm_to_px(intersection), (d.npx_fast() / 2, d.npx_slow() / 2))
    print "OK"
