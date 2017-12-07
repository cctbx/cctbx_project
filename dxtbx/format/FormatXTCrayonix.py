from __future__ import absolute_import, division
from dxtbx.format.FormatXTC import FormatXTC
from xfel.cxi.cspad_ana import cspad_tbx, rayonix_tbx
from libtbx.phil import parse


locator_scope = parse("""
  locator = None
    .type = str
    .help = file format as specified at LCLS,eg. exp=mfxo1916:run=20:smd
  detector_address = None
    .type = str
    .help = detector used for collecting the data at LCLS
""")


class FormatXTCrayonix(FormatXTC):

  def __init__(self, image_file, **kwargs):
    assert(self.understand(image_file))
    self._ds = FormatXTCrayonix._get_datasource(image_file,'xtc')
    self.indices = []
    self.populate_events()
    FormatXTC.__init__(self, image_file, **kwargs)

  @staticmethod
  def _get_event(class_instance,index):
    return class_instance.indices[index]

  def _start(self):
    retEvt = None
    self.n_images = len(self.indices)

  def populate_events(self):
    for nevent,evt in enumerate(self._ds.events()):
      self.indices.append(evt)

  @staticmethod
  def _get_datasource(image_file, ds_type='xtc'):
    import psana
    from psana import DataSource

    user_input = parse(file_name=image_file)
    working_phil = locator_scope.fetch(sources=[user_input])
    params = working_phil.extract()
    if params.locator is None:
      return False
    else:
      img = params.locator
    return DataSource(img)

  @staticmethod
  def understand(image_file):
    import psana
    try:
      user_input = parse(file_name=image_file)
      working_phil = locator_scope.fetch(sources=[user_input])
      params = working_phil.extract()
      if FormatXTC._src is None:
        FormatXTC._src = names[int(raw_input("Please Enter name of detector numbered 1 through %d : "%(len(names))))-1][0]
      if 'rayonix' in FormatXTC._src.lower():
        return True
      else:
        return False
      understood = True
    except IOError,e:
      return False
    return understood

  def get_raw_data(self,index):
    import psana
    from scitbx.array_family import flex
    det = psana.Detector(self._src, self._env)
    data = rayonix_tbx.get_data_from_psana_event(FormatXTCrayonix._get_event(self,index), self._src)
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
    evt = self.indices[index]
    return self._beam_factory.simple(cspad_tbx.evt_wavelength(evt))

  def get_goniometer(self, index=None):
      return None

  def get_scan(self, index=None):
      return None

  def get_mask(self, index=None, goniometer=None):
    return None

  def get_detectorbase(self, index=None):
    print 'get_detectorbase Overload!'

  def get_image_file(self, index=None):
    print 'get_image_file Overload!'

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
      trusted_range = (-1, 65535),
      mask = [])

if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print FormatXTCrayonix.understand(arg)
