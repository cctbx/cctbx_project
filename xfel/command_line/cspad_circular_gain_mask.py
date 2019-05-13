from __future__ import division
from __future__ import print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cspad.circular_gain_mask
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

"""
This command line function generates a gain ascii file suitable for use
by CXI.  The pixels out to the specified resolution in a circular pattern
will be set to low gain.
"""

import sys, numpy, math
import libtbx.phil
from libtbx.utils import Sorry
from xfel.command_line.make_mask import point_inside_circle

master_phil = libtbx.phil.parse("""
resolution = None
  .type = float
  .help = Low gain pixels will be set out to this resolution. If using an annulus, instead, pixels higher than this resolution will be set to low gain
annulus_inner = None
  .type = float
  .help = Use a low gain annulus instead of masking out all the pixels below the given resolution.
annulus_outer = None
  .type = float
  .help = Use a low gain annulus instead of masking out all the pixels below the given resolution.
distance = None
  .type = float
  .help = Detector distance
wavelength = None
  .type = float
  .help = Beam wavelength
out = circle.gain
  .type = str
  .help = Output file path
optical_metrology_path = None
  .type = str
  .help = Path to slac optical metrology file. If not set, use Run 4 metrology
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    try :
      user_phil.append(libtbx.phil.parse(arg))
    except RuntimeError as e :
      raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  assert params.resolution is not None or (params.annulus_inner is not None and params.annulus_outer is not None)
  assert params.distance is not None
  assert params.wavelength is not None

  annulus = (params.annulus_inner is not None and params.annulus_outer is not None)

  if annulus and params.resolution is not None:
    assert params.resolution < params.annulus_outer

  if annulus:
    if params.resolution is None:
      print("Generating annular gain mask between %f and %f angstroms, assuming a distance %s mm and wavelength %s angstroms" % \
        (params.annulus_inner, params.annulus_outer, params.distance, params.wavelength))
    else:
      print("Generating annular gain mask between %f and %f angstroms, assuming a distance %s mm and wavelength %s angstroms. Also, pixels higher than %f angstroms will be set to low gain." % \
        (params.annulus_inner, params.annulus_outer, params.distance, params.wavelength, params.resolution))
  elif params.resolution is not None:
    print("Generating circular gain mask %s angstroms, assuming a distance %s mm and wavelength %s angstroms" % \
      (params.resolution, params.distance, params.wavelength))

  from xfel.cftbx.detector.cspad_cbf_tbx import read_slac_metrology, get_cspad_cbf_handle
  from dxtbx.format.FormatCBFCspad import FormatCBFCspadInMemory
  from dxtbx.model import BeamFactory
  metro = read_slac_metrology(path = params.optical_metrology_path)
  cbf = get_cspad_cbf_handle(None, metro, 'cbf', None, "test", None, params.distance, verbose = True, header_only = True)
  img = FormatCBFCspadInMemory(cbf)
  beam = BeamFactory.simple(params.wavelength).get_s0()
  beam_center = (0,0)

  data = numpy.zeros((11840,194))

  if annulus:
    inner = params.distance * math.tan(2*math.sinh(params.wavelength/(2*params.annulus_inner)))
    outer = params.distance * math.tan(2*math.sinh(params.wavelength/(2*params.annulus_outer)))
    print("Inner (mm):", inner)
    print("Outer (mm):", outer)
  if params.resolution is not None:
    radius = params.distance * math.tan(2*math.sinh(params.wavelength/(2*params.resolution)))
    print("Radius (mm):", radius)

  print("Panel:", end=' '); sys.stdout.flush()
  for p_id, panel in enumerate(img.get_detector()):
    print(p_id, end=' '); sys.stdout.flush()
    for y in range(185):
      for x in range(194):
        lx, ly, lz = panel.get_pixel_lab_coord((x,y))
        if annulus:
          if not point_inside_circle(lx,ly,beam_center[0],beam_center[1],outer) or point_inside_circle(lx,ly,beam_center[0],beam_center[1],inner):
            data[(p_id*185)+y,x] = 1
        if params.resolution is not None:
          if annulus:
            if not point_inside_circle(lx,ly,beam_center[0],beam_center[1],radius):
              data[(p_id*185)+y,x] = 0
          else:
            if not point_inside_circle(lx,ly,beam_center[0],beam_center[1],radius):
              data[(p_id*185)+y,x] = 1
  print()
  print("Masked out %d pixels out of %d (%.2f%%)"%(data.size-int(data.sum()), data.size, 100*(data.size-int(data.sum()))/data.size))
  numpy.savetxt(params.out, data, fmt="%d")
