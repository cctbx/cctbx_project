from __future__ import division
from __future__ import print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cxi.generate_circular_gain_mask
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
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
from xfel.cxi.cspad_ana.parse_calib import calib2sections
from scitbx.array_family import flex
from iotbx.detectors.npy import NpyImage
from spotfinder.applications.xfel import cxi_phil
from xfel.command_line.make_mask import point_inside_circle

master_phil = libtbx.phil.parse("""
detector_format_version = None
  .type = str
  .help = Detector metrology on which to overlay the gain map
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
      print("Generating annular gain mask using %s metrology between %f and %f angstroms, assuming a distance %s mm and wavelength %s angstroms" % \
        (str(params.detector_format_version), params.annulus_inner, params.annulus_outer, params.distance, params.wavelength))
    else:
      print("Generating annular gain mask using %s metrology between %f and %f angstroms, assuming a distance %s mm and wavelength %s angstroms. Also, pixels higher than %f angstroms will be set to low gain." % \
        (str(params.detector_format_version), params.annulus_inner, params.annulus_outer, params.distance, params.wavelength, params.resolution))
  elif params.resolution is not None:
    print("Generating circular gain mask using %s metrology at %s angstroms, assuming a distance %s mm and wavelength %s angstroms" % \
      (str(params.detector_format_version), params.resolution, params.distance, params.wavelength))

  from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp, cbcaa, pixel_size, CsPadDetector
  from iotbx.detectors.cspad_detector_formats import address_and_timestamp_from_detector_format_version
  from xfel.command_line.convert_gain_map import fake_env, fake_config, fake_evt, fake_cspad_ElementV2
  address, timestamp = address_and_timestamp_from_detector_format_version(params.detector_format_version)
  timestamp = evt_timestamp((timestamp,0))

  raw_data = numpy.zeros((11840,194))
  if params.detector_format_version is not None and "XPP" in params.detector_format_version:
    from xfel.cxi.cspad_ana.cspad_tbx import xpp_active_areas
    active_areas = xpp_active_areas[params.detector_format_version]['active_areas']
    data = flex.int(flex.grid((1765,1765)))
    beam_center = (1765 // 2, 1765 // 2)
  else:
    if params.optical_metrology_path is None:
      calib_dir = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0")
      sections = calib2sections(calib_dir)
    else:
      sections = calib2sections(params.optical_metrology_path)
    asic_start = 0
    data3d = []
    for i_quad in range(4):
      asic_size = 185 * 194
      section_size = asic_size * 4
      quad_start = i_quad * section_size * 4
      quad_asics = []
      for i_2x2 in range(4):
        for i_asic in range(2):
          asic_end = asic_start + 185
          a = raw_data[asic_start:asic_end, :]
          asic_start = asic_end

          asic_end = asic_start + 185
          b = raw_data[asic_start:asic_end, :]
          asic_start = asic_end

          quad_asics.append(numpy.concatenate((a,b),axis=1))
      quad_data = numpy.dstack(quad_asics)
      quad_data = numpy.rollaxis(quad_data, 2,0)
      data3d.append(fake_cspad_ElementV2(quad_data, i_quad))

    env = fake_env(fake_config())
    evt = fake_evt(data3d)
    beam_center, active_areas = cbcaa(fake_config(),sections)
    data = flex.int(CsPadDetector(address, evt, env, sections).astype(numpy.float64))

  img_dict = dpack(
          active_areas=active_areas,
          address=address,
          beam_center_x=beam_center[0]*pixel_size,
          beam_center_y=beam_center[1]*pixel_size,
          data=data,
          distance=params.distance,
          pixel_size=pixel_size,
          timestamp=timestamp,
          wavelength=params.wavelength)

  img = NpyImage("", source_data=img_dict)

  args = ["distl.detector_format_version=%s"%params.detector_format_version]
  horizons_phil = cxi_phil.cxi_versioned_extract(args)

  img.readHeader(horizons_phil)
  img.translate_tiles(horizons_phil)
  tm = img.get_tile_manager(horizons_phil)
  effective_active_areas = tm.effective_tiling_as_flex_int()

  if annulus:
    inner = params.distance * math.tan(2*math.sinh(params.wavelength/(2*params.annulus_inner)))/pixel_size
    outer = params.distance * math.tan(2*math.sinh(params.wavelength/(2*params.annulus_outer)))/pixel_size
    print("Pixel inner:", inner)
    print("Pixel outer:", outer)
  if params.resolution is not None:
    radius = params.distance * math.tan(2*math.sinh(params.wavelength/(2*params.resolution)))/pixel_size
    print("Pixel radius:", radius)

  print("Percent done: 0", end=' '); sys.stdout.flush()
  next_percent = 10
  for y in range(data.focus()[1]):
    if y*100/data.focus()[1] > next_percent:
      print(next_percent, end=' '); sys.stdout.flush()
      next_percent += 10
    for x in range(data.focus()[0]):
      if annulus:
        if not point_inside_circle(x,y,beam_center[0],beam_center[1],outer) or point_inside_circle(x,y,beam_center[0],beam_center[1],inner):
          data[y,x] = 1
      if params.resolution is not None:
        if annulus:
          if not point_inside_circle(x,y,beam_center[0],beam_center[1],radius):
            data[y,x] = 0
        else:
          if not point_inside_circle(x,y,beam_center[0],beam_center[1],radius):
            data[y,x] = 1
  print(100)

  if 'XPP' in params.detector_format_version:
    rotations = xpp_active_areas[params.detector_format_version]['rotations']
    angles = [(x+1) * 90 for x in rotations]
  else:
    angles = []
    for quad in sections:
      for section in quad:
        angles.append(section.angle)
        angles.append(section.angle)
  assert len(angles) == int(len(effective_active_areas)/4)

  raw_data = flex.int(flex.grid((11840,194)))
  for i in range(int(len(effective_active_areas)/4)):
    ul_slow = effective_active_areas[4 * i + 0]
    ul_fast = effective_active_areas[4 * i + 1]
    lr_slow = effective_active_areas[4 * i + 2]
    lr_fast = effective_active_areas[4 * i + 3]

    angle = angles[i]

    block = data.matrix_copy_block(
      i_row=ul_slow,i_column=ul_fast,
      n_rows=lr_slow-ul_slow,n_columns=lr_fast-ul_fast)

    if block.focus()[0] > block.focus()[1]:
      block = block.matrix_rot90(1 + (int(round(angle / 90.0)) % 4))
    else:
      block = block.matrix_rot90(-1 + (int(round(angle / 90.0)) % 4))

    raw_data.matrix_paste_block_in_place(
      block = block,
      i_row = i*block.focus()[0],
      i_column = 0
    )

  numpy.savetxt(params.out, raw_data.as_numpy_array(), fmt="%d")
