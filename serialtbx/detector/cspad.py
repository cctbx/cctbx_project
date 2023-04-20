from __future__ import division

# The CAMP and CSpad counters are both 14 bits wide (Strüder et al
# 2010; Philipp et al., 2007), which means the physical limit is 2**14 - 1.
# However, in practice, when the pixels are in the low gain mode, after
# correcting by a gain value of around 6.87, the pixels tend to saturate
# around 90000. See xpp experiment xppe0314, run 184 as evidence.
cspad_saturated_value = 90000

# The dark average for the CSPAD detector is around 1100-1500. A pixel
# histogram of a minimum projection of an uncorrected (raw) light run shows
# a mostly flat tail up to ~800 ADU with a few bumps in the tail which
# represent true underloads. Assume a dark average of 1200 ADU. After dark
# subtraction, 800 - 1200 gives a minimum trusted value of -400. Reject
# pixels less than this.
cspad_min_trusted_value = -400

# The pixel size in mm.  The pixel size is fixed and square, with side
# length of 110 µm (Philipp et al., 2007).  XXX Should really clarify
# this with Sol and Chris.
#
# XXX Andor: 13.5 µm square, CAMP: 75 µm, square (Strüder et al.,
# 2010)
pixel_size = 110e-3

def dpack(active_areas=None,
          address=None,
          beam_center_x=None,
          beam_center_y=None,
          ccd_image_saturation=None,
          data=None,
          distance=None,
          pixel_size=pixel_size,
          saturated_value=None,
          timestamp=None,
          wavelength=None,
          xtal_target=None,
          min_trusted_value=None):
  """XXX Check completeness.  Should fill in sensible defaults."""

  # Must have data.
  if data is None:
    return None

  # Create a time stamp of the current time if none was supplied.
  if timestamp is None:
    timestamp = evt_timestamp()

  # For unknown historical reasons, the dictionary must contain both
  # CCD_IMAGE_SATURATION and SATURATED_VALUE items.
  if ccd_image_saturation is None:
    if saturated_value is None:
      ccd_image_saturation = cspad_saturated_value
    else:
      ccd_image_saturation = saturated_value
  if saturated_value is None:
    saturated_value = ccd_image_saturation

  # Use a minimum value if provided for the pixel range
  if min_trusted_value is None:
    min_trusted_value = cspad_min_trusted_value

  # By default, the beam center is the center of the image.  The slow
  # (vertical) and fast (horizontal) axes correspond to x and y,
  # respectively.
  if beam_center_x is None:
    beam_center_x = pixel_size * data.focus()[1] / 2
  if beam_center_y is None:
    beam_center_y = pixel_size * data.focus()[0] / 2

  # By default, the entire detector image is an active area.  There is
  # no sensible default for distance nor wavelength.  XXX But setting
  # wavelength to zero may be disastrous?
  if active_areas is None:
    # XXX Verify order with non-square detector
    active_areas = flex.int((0, 0, data.focus()[0], data.focus()[1]))
  if distance is None:
    distance = 0
  if wavelength is None:
    wavelength = 0

  # The size must match the image dimensions.  The length along the
  # slow (vertical) axis is SIZE1, the length along the fast
  # (horizontal) axis is SIZE2.
  return {'ACTIVE_AREAS': active_areas,
          'BEAM_CENTER_X': beam_center_x,
          'BEAM_CENTER_Y': beam_center_y,
          'CCD_IMAGE_SATURATION': ccd_image_saturation,
          'DATA': data,
          'DETECTOR_ADDRESS': address,
          'DISTANCE': distance,
          'PIXEL_SIZE': pixel_size,
          'SATURATED_VALUE': saturated_value,
          'MIN_TRUSTED_VALUE': min_trusted_value,
          'SIZE1': data.focus()[0],
          'SIZE2': data.focus()[1],
          'TIMESTAMP': timestamp,
          'SEQUENCE_NUMBER': 0, # XXX Deprecated
          'WAVELENGTH': wavelength,
          'xtal_target': xtal_target}
