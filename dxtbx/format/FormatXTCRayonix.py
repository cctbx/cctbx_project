from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatXTC import FormatXTC, locator_str
from libtbx.phil import parse

try:
  from xfel.cxi.cspad_ana import cspad_tbx, rayonix_tbx
except ImportError:
  # xfel not configured
  pass

rayonix_locator_str = """
  rayonix {
    bin_size = 2
      .type = int
      .help = Detector binning mode
  }
"""

rayonix_locator_scope = parse(rayonix_locator_str+locator_str, process_includes=True)

class FormatXTCRayonix(FormatXTC):

  def __init__(self, image_file, **kwargs):
    assert(self.understand(image_file))
    FormatXTC.__init__(self, image_file, locator_scope = rayonix_locator_scope, **kwargs)
    self._ds = self._get_datasource(image_file)
    self._env = self._ds.env()
    self.events_list = []
    self.populate_events()
    self.n_images = len(self.events_list)

  def _get_event(self,index):
    return self.events_list[index]

  def populate_events(self):
    for nevent,evt in enumerate(self._ds.events()):
      self.events_list.append(evt)

  @staticmethod
  def understand(image_file):
    import psana
    try:
      if FormatXTC._src is None:
        FormatXTC._src = names[int(raw_input("Please Enter name of detector numbered 1 through %d : "%(len(names))))-1][0]
      if 'rayonix' in FormatXTC._src.lower():
        return True
      return False
    except Exception:
      return False

  def get_raw_data(self,index):
    import psana
    from scitbx.array_family import flex
    det = psana.Detector(self._src, self._env)
    data = rayonix_tbx.get_data_from_psana_event(self._get_event(index), self._src)
    return flex.double(data)

  def get_num_images(self):
    return self.n_images

  def get_detector(self, index=None):
    return self._detector()

  def get_beam(self, index=None):
    return self._beam(index)

  def _beam(self, index=None):
    '''Returns a simple model for the beam '''
    if index is None: index=0
    evt = self.events_list[index]
    return self._beam_factory.simple(cspad_tbx.evt_wavelength(evt))

  def get_goniometer(self, index=None):
      return None

  def get_scan(self, index=None):
      return None

  def _detector(self):
    import psana
    pixel_size = rayonix_tbx.get_rayonix_pixel_size(self.params.rayonix.bin_size)
    return self._detector_factory.simple(
      sensor = 'UNKNOWN',
      distance = 100.0,
      beam_centre = (50., 50.),
      fast_direction = '+x',
      slow_direction = '-y',
      pixel_size = (pixel_size, pixel_size),
      image_size = rayonix_tbx.get_rayonix_detector_dimensions(self.params.rayonix.bin_size),
      trusted_range = (rayonix_tbx.rayonix_min_trusted_value, rayonix_tbx.rayonix_saturated_value),
      mask = [])

if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print(FormatXTCRayonix.understand(arg))
