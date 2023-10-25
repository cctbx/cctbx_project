from __future__ import division

# given value of rayonix detector saturation xppi6115
rayonix_saturated_value = 2**16 -1

# minimum value for rayonix data
rayonix_min_trusted_value = 0
rayonix_max_trusted_value = rayonix_saturated_value - 1

def get_rayonix_pixel_size(bin_size):
  ''' Given a bin size determine a pixel size.

Michael Blum from Rayonix said The pixel size is recorded in the header,
but can be derived trivially from the overall dimension of the corrected imaging
area (170mm) and the number of pixels. (3840 unbinned). The corrected image is
forced to this size.

unbinned 170/3840  = 0.04427

I believe the accuracy of the MEAN pixel size to be at least as good as 0.1%
which is the limit to which I can measure our calibration plate and exceeds the
 parallax error in our calibration station.

Note, the Rayonix MX340 has the same pixel size as the MX170:

unbinned 340/7680  = 0.04427

  @param bin_size rayonix bin size as an integer
  '''
  pixel_size=bin_size*170/3840
  return pixel_size

def get_rayonix_detector_dimensions(env):
  ''' Given a psana env object, find the detector dimensions
      @param env psana environment object
  '''
  import psana
  cfgs = env.configStore()
  rayonix_cfg = cfgs.get(psana.Rayonix.ConfigV2, psana.Source('Rayonix'))
  if not rayonix_cfg: return None, None
  return rayonix_cfg.width(), rayonix_cfg.height()

def get_data_from_psana_event(evt, address):
  """ Read the pixel data for a Rayonix image from an event
  @param psana event object
  @param address old style psana detector address
  @return numpy array with raw data"""
  from psana import Source, Camera
  from serialtbx.detector import xtc
  import numpy as np
  address = xtc.old_address_to_new_address(address)
  src=Source(address)
  data = evt.get(Camera.FrameV1,src)
  if data is not None:
    data = data.data16().astype(np.float64)
  return data
