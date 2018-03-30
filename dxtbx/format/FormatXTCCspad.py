from __future__ import absolute_import, division
from dxtbx.format.FormatXTC import FormatXTC,locator_str
from libtbx.phil import parse
try:
  from xfel.cxi.cspad_ana import cspad_tbx
  from xfel.cftbx.detector import cspad_cbf_tbx
except ImportError:
  # xfel not configured
  pass

cspad_locator_str = """
  cspad {
    apply_gain_mask = True
      .type = bool
      .help = flag to indicate if gain should be applied to cspad data
    dark_correction = True
      .type = bool
      .help = flag to decide if dark correction should be done
    }
"""

cspad_locator_scope = parse(cspad_locator_str+locator_str, process_includes=True)

class FormatXTCCspad(FormatXTC):

  def __init__(self, image_file, **kwargs):
    assert(self.understand(image_file))
    self._ds = self._get_datasource(image_file)
    self.events_list = []
    self.populate_events()
    self.params = FormatXTC.params_from_phil(cspad_locator_scope, image_file)
    #from IPython import embed; embed(); exit()
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
      if 'cspad' in FormatXTC._src.lower():
        return True
      return False
    except Exception:
      return False

  def get_raw_data(self,index):
    import psana
    from scitbx.array_family import flex
    import numpy as np
    det = psana.Detector(self._src, self._env)
    d = self.get_detector()
    data = cspad_cbf_tbx.get_psana_corrected_data(det, self._get_event(index),
                                                  use_default=False,
                                                  dark=self.params.cspad.dark_correction,
                                                  common_mode=None,
                                                  apply_gain_mask=self.params.cspad.apply_gain_mask,
                                                  gain_mask_value=None,
                                                  per_pixel_gain=False)
    data = data.astype(np.float64)
    self._raw_data = []
    for quad_count, quad in enumerate(d.hierarchy()):
      for sensor_count,sensor in enumerate(quad):
        for asic_count,asic in enumerate(sensor):
          fdim,sdim =  asic.get_image_size()
          asic_data = data[sensor_count+quad_count*8, :, asic_count*fdim:(asic_count+1)*fdim] # 8 sensors per quad
          self._raw_data.append(flex.double(np.array(asic_data)))
    assert len(d) == len(self._raw_data)
    return tuple(self._raw_data)

  def get_num_images(self):
    return self.n_images


  def get_detector(self, index=None):
    return self._detector(index)

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

#XXX Implement recursive version
  def _detector(self, index=None):
    import psana
    from xfel.cftbx.detector.cspad_cbf_tbx import read_slac_metrology
    from dxtbx.model import Detector
    from scitbx.matrix import col
    from dxtbx.model import ParallaxCorrectedPxMmStrategy
    if index is None: index = 0
    self._env = self._ds.env()
    self._det = psana.Detector(self._src,self._env)
    geom=self._det.pyda.geoaccess(self.events_list[index])
    cob = read_slac_metrology(geometry=geom, include_asic_offset=True)
    d = Detector()
    pg0 = d.hierarchy()
    # first deal with D0
    det_num = 0
    origin = col((cob[(0,)] * col((0,0,0,1)))[0:3])
    fast   = col((cob[(0,)] * col((1,0,0,1)))[0:3]) - origin
    slow   = col((cob[(0,)] * col((0,1,0,1)))[0:3]) - origin
    origin += col((0., 0., -100.))
    pg0.set_local_frame(fast.elems,slow.elems,origin.elems)
    pg0.set_name('D%d'%(det_num))
    for quad_num in xrange(4):
      # Now deal with Qx
      pg1 = pg0.add_group()
      origin = col((cob[(0,quad_num)] * col((0,0,0,1)))[0:3])
      fast   = col((cob[(0,quad_num)] * col((1,0,0,1)))[0:3]) - origin
      slow   = col((cob[(0,quad_num)] * col((0,1,0,1)))[0:3]) - origin
      pg1.set_local_frame(fast.elems,slow.elems,origin.elems)
      pg1.set_name('D%dQ%d'%(det_num, quad_num))
      for sensor_num in xrange(8):
      # Now deal with Sy
        pg2=pg1.add_group()
        origin = col((cob[(0,quad_num,sensor_num)] * col((0,0,0,1)))[0:3])
        fast   = col((cob[(0,quad_num,sensor_num)] * col((1,0,0,1)))[0:3]) - origin
        slow   = col((cob[(0,quad_num,sensor_num)] * col((0,1,0,1)))[0:3]) - origin
        pg2.set_local_frame(fast.elems,slow.elems,origin.elems)
        pg2.set_name('D%dQ%dS%d'%(det_num,quad_num,sensor_num))
        # Now deal with Az
        for asic_num in xrange(2):
          val = 'ARRAY_D0Q%dS%dA%d'%(quad_num,sensor_num,asic_num)
          p = pg2.add_panel()
          origin = col((cob[(0,quad_num,sensor_num, asic_num)] * col((0,0,0,1)))[0:3])
          fast   = col((cob[(0,quad_num,sensor_num, asic_num)] * col((1,0,0,1)))[0:3]) - origin
          slow   = col((cob[(0,quad_num,sensor_num, asic_num)] * col((0,1,0,1)))[0:3]) - origin
          p.set_local_frame(fast.elems,slow.elems,origin.elems)
          p.set_pixel_size((cspad_cbf_tbx.pixel_size, cspad_cbf_tbx.pixel_size))
          p.set_image_size(cspad_cbf_tbx.asic_dimension)
          p.set_trusted_range((cspad_tbx.cspad_min_trusted_value, cspad_tbx.cspad_saturated_value))
          p.set_name(val)

    try:
      beam = self._beam(index)
    except Exception:
      print 'No beam object initialized. Returning CSPAD detector without parallax corrections'
      return d

    # take into consideration here the thickness of the sensor also the
    # wavelength of the radiation (which we have in the same file...)
    wavelength = beam.get_wavelength()
    thickness = 0.5 # mm, see Hart et al. 2012
    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    # mu_at_angstrom returns cm^-1
    mu = table.mu_at_angstrom(wavelength) / 10.0 # mu: mm^-1
    t0 = thickness
    for panel in d:
      panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
    return d



if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print FormatXTCCspad.understand(arg)
