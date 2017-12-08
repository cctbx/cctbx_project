from __future__ import absolute_import, division
from dxtbx.format.FormatXTC import FormatXTC
from xfel.cxi.cspad_ana import cspad_tbx, rayonix_tbx

class FormatXTCRayonix(FormatXTC):

  def __init__(self, image_file, **kwargs):
    assert(self.understand(image_file))
    self._ds = self._get_datasource(image_file)
    self.events_list = []
    self.populate_events()
    FormatXTC.__init__(self, image_file, **kwargs)

  def _get_event(self,index):
    return self.events_list[index]

  def _start(self):
    self.n_images = len(self.events_list)

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
    except Exception,e:
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
    self._env = self._ds.env()
    self._det = psana.Detector(self._src,self._env)
    return self._detector_factory.simple(
      sensor = 'UNKNOWN',
      distance = 100.0,
      beam_centre = (50., 50.),
      fast_direction = '+x',
      slow_direction = '-y',
      pixel_size = (rayonix_tbx.get_rayonix_pixel_size(2), rayonix_tbx.get_rayonix_pixel_size(2)),
      image_size = self._det.shape(),
      trusted_range = (rayonix_tbx.rayonix_min_trusted_value, rayonix_tbx.rayonix_saturated_value),
      mask = [])

if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print FormatXTCRayonix.understand(arg)
