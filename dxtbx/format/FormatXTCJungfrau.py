from __future__ import absolute_import, division,print_function
from dxtbx.format.FormatXTC import FormatXTC, locator_str
from libtbx.phil import parse
import numpy as np
try:
  from xfel.cxi.cspad_ana import cspad_tbx
except ImportError:
  # xfel not configured
  pass

jungfrau_locator_str = """
  jungfrau {
    dark = True
      .type = bool
      .help = Dictates if dark subtraction is done from raw data
    monolithic = False
      .type = bool
      .help = switch to FormatXTCJungfrauMonolithic if True. Used for LS49 image averaging
  }
"""

jungfrau_locator_scope = parse(jungfrau_locator_str+locator_str, process_includes=True)

class FormatXTCJungfrau(FormatXTC):

  def __init__(self, image_file, **kwargs):
    assert(self.understand(image_file))
    FormatXTC.__init__(self, image_file, locator_scope = jungfrau_locator_scope, **kwargs)
    self._ds = FormatXTCJungfrau._get_datasource(image_file, self.params)
    self._env = self._ds.env()
    self.populate_events()
    self.n_images = len(self.times)
    self._cached_detector = {}
    self._cached_psana_detectors = {}

  @staticmethod
  def understand(image_file):
    try:
      params = FormatXTC.params_from_phil(jungfrau_locator_scope,image_file)
    except Exception:
      return False
    ds = FormatXTC._get_datasource(image_file, params)
    return any(['jungfrau' in src.lower() for src in params.detector_address])

  def get_raw_data(self,index):
    import psana
    from scitbx.array_family import flex
    d = FormatXTCJungfrau.get_detector(self, index)
    evt = self._get_event(index)
    run = self.get_run_from_index(index)
    if run.run() not in self._cached_psana_detectors:
      assert len(self.params.detector_address) == 1
      self._cached_psana_detectors[run.run()] = psana.Detector(self.params.detector_address[0], self._env)
    det = self._cached_psana_detectors[run.run()]
    data = det.calib(evt)
    data = data.astype(np.float64)
    self._raw_data = []
    for quad_count, quad in enumerate(d.hierarchy()):
      for asic_count,asic in enumerate(quad):
        fdim,sdim =  asic.get_image_size()
        sensor_id = asic_count //4 # There are 2X4 asics per quadrant
        asic_in_sensor_id = asic_count%4 # this number will be 0,1,2 or 3
        asic_data = data[quad_count][sensor_id*sdim:(sensor_id+1)*sdim,asic_in_sensor_id*fdim:(asic_in_sensor_id+1)*fdim] # 8 sensors per quad
        self._raw_data.append(flex.double(np.array(asic_data)))
    assert len(d) == len(self._raw_data)
    return tuple(self._raw_data)

  def get_num_images(self):
    return self.n_images

  def get_detector(self, index=None):
    return FormatXTCJungfrau._detector(self, index)

  def get_beam(self, index=None):
    return self._beam(index)

  def _beam(self, index=None):
    '''Returns a simple model for the beam '''
    if index is None: index=0
    evt = self._get_event(index)
    wavelength = cspad_tbx.evt_wavelength(evt)
    if wavelength is None: return None
    return self._beam_factory.simple(wavelength)

  def get_goniometer(self, index=None):
      return None

  def get_scan(self, index=None):
      return None

  def _detector(self, index=None):
    run = self.get_run_from_index(index)
    if run.run() in self._cached_detector:
      return self._cached_detector[run.run()]
    import psana
    from dxtbx.model import Detector
    from scitbx.matrix import col
    if index is None: index = 0
    self._env = self._ds.env()
    assert len(self.params.detector_address) == 1
    self._det = psana.Detector(self.params.detector_address[0],self._env)
    geom=self._det.pyda.geoaccess(self._get_event(index).run())
    pixel_size = self._det.pixel_size(self._get_event(index))/1000.0 # convert to mm
    d = Detector()
    pg0 = d.hierarchy()
    # first deal with D0
    det_num = 0
    D0= geom.get_top_geo().get_list_of_children()[0]
    xx,yy,zz = D0.get_pixel_coords()
    xx = xx/1000.0 # to mm
    yy = yy/1000.0 # to mm
    zz = zz/1000.0 # to mm
    oriD0 = col((np.mean(xx),np.mean(yy), -np.mean(zz)))
    fp = col((xx[0][0][1],yy[0][0][1], zz[0][0][1]))
    sp = col((xx[0][1][0], yy[0][1][0], zz[0][1][0]))
    op = col((xx[0][0][0], yy[0][0][0], zz[0][0][0]))
    origin = oriD0
    fast   = (fp-op).normalize()
    slow   = (sp-op).normalize()
    pg0.set_local_frame(fast.elems,slow.elems,origin.elems)
    pg0.set_name('D%d'%(det_num))

    # Now deal with Qx
    for quad_num in xrange(2):
      pg1 = pg0.add_group()
      Qx = D0.get_list_of_children()[quad_num]
      xx,yy,zz = Qx.get_pixel_coords()
      xx = xx/1000.0 # to mm
      yy = yy/1000.0 # to mm
      zz = zz/1000.0 # to mm
      oriQx = col((np.mean(xx), np.mean(yy), np.mean(zz)))
      fp = col((xx[0][1],yy[0][1], zz[0][1]))
      sp = col((xx[1][0], yy[1][0], zz[1][0]))
      op = col((xx[0][0], yy[0][0], zz[0][0]))
      origin = oriQx
      fast   = (fp-op).normalize()
      slow   = (sp-op).normalize()
      pg1.set_local_frame(fast.elems,slow.elems,origin.elems)
      pg1.set_name('D%dQ%d'%(det_num, quad_num))

      # Now deal with Az
      for asic_num in xrange(8):
        val = 'ARRAY_D0Q%dA%d'%(quad_num,asic_num)
        p = pg1.add_panel()
        dim_slow = xx.shape[0]
        dim_fast = xx.shape[1]
        sensor_id = asic_num //4 # There are 2X4 asics per quadrant
        asic_in_sensor_id = asic_num%4 # this number will be 0,1,2 or 3
        id_slow = sensor_id*(dim_slow//2)
        id_fast = asic_in_sensor_id*(dim_fast//4)
        oriAy = col((xx[id_slow][id_fast],yy[id_slow][id_fast], zz[id_slow][id_fast]))
        fp = col((xx[id_slow][id_fast+1],yy[id_slow][id_fast+1], zz[id_slow][id_fast+1]))
        sp = col((xx[id_slow+1][id_fast],yy[id_slow+1][id_fast], zz[id_slow+1][id_fast]))
        origin = oriAy - oriQx
        fast   = (fp - oriAy).normalize()
        slow   = (sp - oriAy).normalize()
        p.set_local_frame(fast.elems,slow.elems,origin.elems)
        p.set_pixel_size((pixel_size,pixel_size))
        p.set_image_size((dim_fast//4,dim_slow//2))
        p.set_trusted_range((-1, 2e6))
        p.set_name(val)
    self._cached_detector[run.run()] = d
    return d

class FormatXTCJungfrauMonolithic(FormatXTCJungfrau):
  ''' Monolithic version of the Jungfrau, I.E. use the psana detector image function to assemble a monolithic image '''
  def __init__(self, image_file, **kwargs):
    assert(self.understand(image_file))
    FormatXTC.__init__(self, image_file, locator_scope = jungfrau_locator_scope, **kwargs)
    self._ds = FormatXTCJungfrauMonolithic._get_datasource(image_file, self.params)
    self._env = self._ds.env()
    self.populate_events()
    self.n_images = len(self.times)
    self._cached_detector = {}
    self._cached_psana_detectors = {}

  @staticmethod
  def understand(image_file):
    try:
      params = FormatXTC.params_from_phil(jungfrau_locator_scope,image_file)
      if params.jungfrau.monolithic:
        return True
      return False
    except Exception:
      return False

  def get_raw_data(self,index):
    import psana
    from scitbx.array_family import flex
    d = self.get_detector(index)
    evt = self._get_event(index)
    run = self.get_run_from_index(index)
    if run.run() not in self._cached_psana_detectors:
      assert len(self.params.detector_address) == 1
      self._cached_psana_detectors[run.run()] = psana.Detector(self.params.detector_address[0], self._env)
    det = self._cached_psana_detectors[run.run()]
    data = det.image(evt)
    data = data.astype(np.float64)
    self._raw_data = flex.double(data)
    return self._raw_data

  def get_detector(self, index=None):
    return self._detector(index)

  def _detector(self, index=None):
    if index is None: index = 0
    return self._detector_factory.simple(
      sensor = 'UNKNOWN',
      distance = 100.0,
      beam_centre = (50., 50.),
      fast_direction = '+x',
      slow_direction = '-y',
      pixel_size = (0.075, 0.075),
      image_size = (1030,1064),
      trusted_range = (-1,2e6),
      mask = [])

if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    # Bug, should call this part differently for understand method to work
    print(FormatXTCJungfrau.understand(arg))
