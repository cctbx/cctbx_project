def compute_image(work_params):
  dpx,dpy = work_params.detector.pixels
  from scitbx.array_family import flex
  assert work_params.noise.max > 0
  mt = flex.mersenne_twister(seed=work_params.noise.random_seed)
  image = mt.random_size_t(
    size=dpx*dpy,
    modulus=work_params.noise.max).as_int()
  image.reshape(flex.grid(dpx,dpy))
  return image

def process(work_params, image):
  from spotfinder import core_toolbox
  options = ""
  report_overloads = True
  dobj = core_toolbox.w_Distl(options, report_overloads)
  dsx,dsy = work_params.detector.size
  dpx,dpy = work_params.detector.pixels
  pixel_size = dsx / dpx
  assert pixel_size == dsy / dpy
  saturation = work_params.signal_max
  dobj.setspotimg(
    pixel_size = pixel_size,
    distance = work_params.detector.distance,
    wavelength = work_params.wavelength,
    xbeam = dsx/2,
    ybeam = dsy/2,
    rawdata = image,
    peripheral_margin = work_params.detector.peripheral_margin,
    saturation = saturation)
  dobj.set_tiling("")
  dobj.set_scanbox_windows(work_params.detector.scanbox_windows)
  dobj.parameter_guarantees()
  dobj.get_underload()
  dobj.pxlclassify()
  dobj.search_icerings()
  dobj.search_maximas()
  dobj.search_spots()
  dobj.search_overloadpatches()
  dobj.finish_analysis()
  #assert dobj.spots.size() == 0

def run(args):
  from libtbx import phil
  phil_str = """\
wavelength = 1
  .type = float
signal_max = 60000
  .type = int
noise {
  max = 10
    .type = int
  random_seed = 0
    .type = int
}
detector {
  distance = 250
    .type = float
  size = 200 200
    .type = floats(size=2)
    .help = "Detector edge length (x,y)"
  pixels = 1000 1000
    .type = ints(size=2)
    .help = "Number of pixels in each detector dimension (x,y)"
  scanbox_windows = 101 51 51
    .type = ints(size=3, value_min=10)
    .help = "Integer scanbox sizes for calculating background, for cycles 1,2, and 3."
  peripheral_margin = 20
    .type = int(value_min=0)
    .help = "No spot detection inside margin; width in pixels"
}
"""
  import libtbx.phil
  master_phil = libtbx.phil.parse(input_string=phil_str)
  argument_interpreter = master_phil.command_line_argument_interpreter()
  updated_params = master_phil.fetch(sources =
   [argument_interpreter.process(arg=arg) for arg in args])
  work_params = updated_params.extract()
  image = compute_image(work_params)
  process(work_params, image)

def run_scanbox_tests():
  # Large image, normal conditions
  run([])

  # SQUARE vs. CIRCLE image shape detection relies on hard-coded minimum size;
  #  to avoid segmentation fault, return UNKNOWN shape for small images (eg, 45x45).
  run("""detector.scanbox_windows=10,10,10 detector.pixels=45,45""".split(" "))

  # Fix error in iotbx/detectors/scanbox.h; a 90x90 pixel image with 20-pixel margin
  #  should yield exactly one 50x50-pixel scanbox (insert SCITBX_EXAMINE to check)
  run("""detector.scanbox_windows=50,50,50 detector.pixels=90,90""".split(" "))

  # A 89x89 pixel image with 20-pixel margin should yield exactly zero scanboxes;
  #  avoid divide-by-zero error by returning from generate_normal_spacing().
  run("""detector.scanbox_windows=50,50,50 detector.pixels=89,89""".split(" "))

  # Function libdistl.cpp get_underload() makes hard-coded assumptions about
  #  minimum image size; avoid the special heuristics if image is < 100x100 pixels
  #  Note the minimum scanbox_window size is implicity hardcoded in libdistl.cpp,
  #  most likely in the diffimage::search_maximas() procedure.
  run("""detector.scanbox_windows=10,10,10 detector.pixels=25,25 peripheral_margin=0""".split(" "))

  print "OK"

if (__name__ == "__main__"):
  import sys
  if len(sys.argv) > 1:
    run(args=sys.argv[1:])
  else:
    run_scanbox_tests()
