def process(work_params, pixels, show_spots=True):
  from spotfinder import core_toolbox
  import time
  t0 = time.time()
  options = ""
  report_overloads = True
  dobj = core_toolbox.w_Distl(options, report_overloads)
  dsx,dsy = work_params.detector.size
  dpx,dpy = work_params.detector.pixels
  pixel_size = dsx / dpx
  assert pixel_size == dsy / dpy
  dobj.setspotimg(
    pixel_size=pixel_size,
    distance=work_params.detector.distance,
    wavelength=work_params.wavelength,
    xbeam=dsx/2,
    ybeam=dsy/2,
    rawdata=pixels,
    peripheral_margin=work_params.spotfinder.peripheral_margin,
    saturation=work_params.signal_max+2*work_params.noise.max)
  dobj.set_tiling("")
  dobj.set_scanbox_windows(work_params.spotfinder.scanbox_windows)
  dobj.parameter_guarantees()
  dobj.get_underload()
  dobj.pxlclassify()
  dobj.search_icerings()
  dobj.search_maximas()
  dobj.search_spots()
  dobj.search_overloadpatches()
  dobj.finish_analysis()
  print "Time spot finding: %.2f" % (time.time()-t0)
  print "Number of spots:", dobj.spots.size()
  if (show_spots):
    print "        Pixel"
    print "   Center of mass     Weight"
    for spot in dobj.spots:
      print "(%8.3f, %8.3f)  %8.2e" % (
        spot.ctr_mass_x(), spot.ctr_mass_y(), spot.total_mass)
  print
  return dobj.spots

def process_args(args, extra_phil_str=""):
  import libtbx.load_env
  libtbx.env.require_module("spotfinder")
  import spotfinder
  from rstbx.simage import create
  return create.process_args(
    args=args,
    extra_phil_str=spotfinder.phil_str+extra_phil_str)
