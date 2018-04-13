from __future__ import absolute_import, division
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
  }
"""

jungfrau_locator_scope = parse(jungfrau_locator_str+locator_str, process_includes=True)

class FormatXTCJungfrau(FormatXTC):

  def __init__(self, image_file, **kwargs):
    assert(self.understand(image_file))
    FormatXTC.__init__(self, image_file, locator_scope = jungfrau_locator_scope, **kwargs)
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
    try:
      if FormatXTC._src is None:
        FormatXTC._src = names[int(raw_input("Please Enter name of detector numbered 1 through %d : "%(len(names))))-1][0]
      if 'jungfrau' in FormatXTC._src.lower():
        return True
      return False
    except Exception:
      return False

  def get_raw_data(self,index):
    import psana
    from scitbx.array_family import flex
    det = psana.Detector(self._src, self._env)
    d = self.get_detector()
    evt = self._get_event(index)
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

  def _detector(self, index=None):
    import psana
    from dxtbx.model import Detector
    from scitbx.matrix import col
    if index is None: index = 0
    self._env = self._ds.env()
    self._det = psana.Detector(self._src,self._env)
    geom=self._det.pyda.geoaccess(self.events_list[index])
    pixel_size = self._det.pixel_size(self.events_list[index])/1000.0 # convert to mm
    d = Detector()
    pg0 = d.hierarchy()
    # first deal with D0
    det_num = 0
    D0= geom.get_top_geo().get_list_of_children()[0]
    yy,xx,zz = D0.get_pixel_coords() # weird convention by psana for this level
    xx = xx/1000.0 # to mm
    yy = yy/1000.0 # to mm
    zz = zz/1000.0 # to mm
    oriD0 = col((np.mean(xx),np.mean(yy), -np.mean(zz)))
    origin = oriD0
    fast   = col((1,0,0))
    slow   = col((0,-1,0))
    pg0.set_local_frame(fast.elems,slow.elems,origin.elems)
    pg0.set_name('D%d'%(det_num))
    for quad_num in xrange(2):
      # Now deal with Qx
      pg1 = pg0.add_group()
      Qx = D0.get_list_of_children()[quad_num]
      xx,yy,zz = Qx.get_pixel_coords()
      xx = xx/1000.0 # to mm
      yy = yy/1000.0 # to mm
      zz = zz/1000.0 # to mm
      oriQx = col((np.mean(xx), np.mean(yy), np.mean(zz) + oriD0[2])) # psana returns 0 for zz
      origin = oriQx - oriD0
      fast   = col((1,0,0))
      slow   = col((0,-1,0))
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
        id_slow = int(dim_slow/(sensor_id+1))-1
        id_fast = asic_in_sensor_id*(dim_fast//4)
        oriAy = col((xx[id_slow][id_fast],yy[id_slow][id_fast], zz[id_slow][id_fast] + oriD0[2]))  # psana returns 0 for zz
        fp = col((xx[id_slow][id_fast+1],yy[id_slow][id_fast+1], zz[id_slow][id_fast+1]+oriD0[2])) # psana returns 0 for zz
        sp = col((xx[id_slow-1][id_fast],yy[id_slow-1][id_fast], zz[id_slow-1][id_fast]+oriD0[2])) # psana returns 0 for zz
        origin = oriAy - oriQx
        fast   = (fp - oriAy).normalize()
        slow   = (sp - oriAy).normalize()
        p.set_local_frame(fast.elems,slow.elems,origin.elems)
        p.set_pixel_size((pixel_size,pixel_size))
        p.set_image_size((dim_fast//4,dim_slow//2))
        p.set_trusted_range((-1, 2e6))
        p.set_name(val)
    return d

if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print FormatXTCJungfrau.understand(arg)
