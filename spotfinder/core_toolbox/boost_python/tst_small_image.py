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
    pixel_size,
    work_params.detector.distance,
    work_params.wavelength,
    dsx/2,
    dsy/2,
    image,
    saturation)
  dobj.set_tiling("")
  dobj.parameter_guarantees()
  dobj.get_underload()
  dobj.pxlclassify()
  dobj.search_icerings()
  dobj.search_maximas()
  dobj.search_spots()
  dobj.search_overloadpatches()
  dobj.finish_analysis()
  assert dobj.spots.size() == 0

def run(args):
  assert len(args) == 0
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
  pixels = 52 52
    .type = ints(size=2)
}
"""
  import libtbx.phil
  master_phil = libtbx.phil.parse(input_string=phil_str)
  work_params = master_phil.extract()
  image = compute_image(work_params)
  process(work_params, image)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
