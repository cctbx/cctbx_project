from __future__ import division
# Utility functions for the Rayonix Detector.

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

  @param bin_size rayonix bin size as an integer
  '''
  pixel_size=bin_size*170/3840
  return pixel_size


def get_rayonix_detector_dimensions(bin_size):
  ''' Given the bin size determine the detector dimensions.
      Based on the number of pixels (unbinned) and the bin size calculate
      integer detector dimensions.

      @param bin_size rayonix bin size as an integer
  '''
  assert 3840%bin_size==0
  return 3840//bin_size,3840//bin_size
