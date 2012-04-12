
# TODO display indices of plane on image
#      adjustable minimum spot size

import cctbx.miller.display
from libtbx.utils import Sorry, Usage
import libtbx.phil
from cStringIO import StringIO
import os
import sys

master_phil = libtbx.phil.parse("""
include scope cctbx.miller.display.master_phil
h = None
  .type = int
k = None
  .type = int
l = None
  .type = int
output_file = None
  .type = path
width = 800
  .type = int
height = 800
  .type = int
format = *Auto png eps
  .type = choice
""", process_includes=True)

def run (args) :
  if (len(args) == 0) :
    params_out = StringIO()
    master_phil.show(out=params_out, prefix="  ")
    raise Usage("""\
cctbx.precession_view data.mtz [ [h=X] | [k=X] | [l=X] ] output_file=img.png

Draws a color image of a reciprocal space plane.  For interactive use, run
phenix.data_viewer, which includes both 2D and 3D views (and can also save
images).

Full parameters:
%s""" % params_out.getvalue())
  import iotbx.phil
  from iotbx import file_reader
  import ImageDraw
  import Image
  phil_mods = libtbx.phil.parse("""
expand_to_p1 = True
expand_anomalous = True
slice_mode = True
""")
  master_phil_mod = master_phil.fetch(source=phil_mods)
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil_mod,
    reflection_file_def="data",
    pdb_file_def="symmetry_file")
  params = cmdline.work.extract()
  hkl_in = file_reader.any_file(params.data, force_type="hkl")
  hkl_in.check_file_type("hkl")
  selected_array = None
  obs_arrays = []
  complex_arrays = []
  for array in hkl_in.file_server.miller_arrays :
    if (params.labels is None) :
      if (array.is_xray_amplitude_array() or
          array.is_xray_intensity_array()) :
        obs_arrays.append(array)
      elif (array.is_complex_array()) :
        complex_arrays.append(array)
    elif (array.info().label_string() == params.labels) :
      selected_array = array
      break
  if (selected_array is None) :
    if (len(obs_arrays) == 1) :
      selected_array = obs_arrays[0]
      print "Using %s since no labels were specified." % \
        selected_array.info().label_string()
    elif (len(obs_arrays) > 1) :
      raise Sorry("""Multiple equally suitable arrays of data found:
  %s
Please choose one by specifying the 'labels' parameter.""" %
  "  ".join([ a.info().label_string() for a in obs_arrays ]))
  params.slice_mode = True
  if (params.h is not None) :
    assert (params.k is None) and (params.l is None)
    params.slice_axis = "h"
    params.slice_index = params.h
  elif (params.k is not None) :
    assert (params.h is None) and (params.l is None)
    params.slice_axis = "k"
    params.slice_index = params.k
  elif (params.l is not None) :
    assert (params.h is None) and (params.k is None)
    params.slice_axis = "l"
    params.slice_index = params.l
  scene = cctbx.miller.display.scene(
    miller_array=selected_array,
    settings=params)
  postcript = False
  if (params.output_file is not None) :
    if (params.output_file.endswith("eps")) :
      if (params.format == "png") :
        raise Sorry("Postscript format requested, but output file ends in .png!")
      postcript = True
    elif (not params.output_file.endswith("png")) :
      raise Sorry("Output file extension must be .eps or .png.")
  else :
    prefix = os.path.basename(os.path.splitext(params.data)[0])
    if (params.format == "Auto") or (params.format == "png") :
      ext = "png"
      postscript = False
    else :
      ext = "eps"
      postscript = True
    if (params.slice_axis == "h") :
      params.output_file = "%s_%dkl.%s" % (prefix, params.slice_index, ext)
    elif (params.slice_axis == "k") :
      params.output_file = "%s_h%dl.%s" % (prefix, params.slice_index, ext)
    else :
      params.output_file = "%s_hk%d.%s" % (prefix, params.slice_index, ext)
  if (params.black_background) :
    bg = (0,0,0)
  else :
    bg = (255,255,255)
  if (postscript) :
    render = render_postscript(
      w=params.width,
      h=params.height,
      scene=scene,
      settings=params)
    canvas = open(params.output_file, "w")
    render.paint(canvas)
    canvas.close()
  else :
    render = render_pil(
      w=params.width*8,
      h=params.height*8,
      scene=scene,
      settings=params)
    im = Image.new('RGBA', (params.width*8, params.height*8), bg)
    canvas = ImageDraw.Draw(im)
    render.render(canvas)
    im2 = im.resize((params.width, params.height), Image.ANTIALIAS)
    im2.save(params.output_file)
  print "Wrote %s" % params.output_file

def frac2int (c) :
  return (int(c[0]*255), int(c[1]*255), int(c[2]*255))

class render_pil (cctbx.miller.display.render_2d) :
  def __init__ (self, w, h, *args, **kwds) :
    self._w = w
    self._h = h
    cctbx.miller.display.render_2d.__init__(self, *args, **kwds)

  def GetSize (self) :
    return (self._w, self._h)

  def get_scale_factor (self) : # XXX why is this necessary?
    return 110.

  def draw_line (self, canvas, x1, y1, x2, y2) :
    canvas.line([(x1,y1),(x2,y2)], fill=frac2int(self._foreground), width=2)

  def draw_text (self, canvas, text, x, y) :
    # FIXME this is coming up blank
    canvas.text((x,y), text, fill=frac2int(self._foreground))

  def draw_open_circle (self, canvas, x, y, radius, color=None) :
    radius = max(radius, 4.)
    if (color is None) : color = self._foreground
    canvas.ellipse(
      xy=[(x-radius, y-radius), (x+radius, y+radius)],
      outline=frac2int(color),
      fill=frac2int(self._background))

  def draw_filled_circle (self, canvas, x, y, radius, color) :
    radius = max(radius, 4.)
    canvas.ellipse(
      xy=[(x-radius, y-radius), (x+radius, y+radius)],
      outline=frac2int(color),
      fill=frac2int(color))

class render_postscript (cctbx.miller.display.render_2d) :
  def __init__ (self, w, h, *args, **kwds) :
    self._w = w
    self._h = h
    cctbx.miller.display.render_2d.__init__(self, *args, **kwds)

  def GetSize (self) :
    return (self._w, self._h)

  def get_scale_factor (self) : # XXX why is this necessary?
    return 110.

  def paint (self, canvas) :
    # XXX I have no idea why this is wrong, but actually getting the image
    # dimensions right appears to be far more difficult
    self.setup_colors()
    canvas.write("%!PS-Adobe-3.0 EPSF-3.0\n")
    canvas.write("%%Pages: 1\n")
    canvas.write("%%BeginPreview:\n")
    canvas.write("%%%%BoundingBox: 0 0 %d %d\n" % self.GetSize())
    canvas.write("<< /PageSize [ %d %d ] /ImagingBBox null >> setpagedevice\n"
      % self.GetSize())
    canvas.write("/Helvetica findfont 14 scalefont setfont\n")
    canvas.write("%g %g %g setrgbcolor\n" % self._background)
    canvas.write("0 0 %d %d rectfill\n" % self.GetSize())
    self.render(canvas)

  def draw_line (self, canvas, x1, y1, x2, y2) :
    w, h = self.GetSize()
    c = self._foreground
    canvas.write("""
newpath
%d %d moveto
%d %d lineto
%g %g %g setrgbcolor
stroke
""" % (x1, w-y1, x2, w-y2, c[0], c[1], c[2]))

  def draw_text (self, canvas, text, x, y) :
    w, h = self.GetSize()
    c = self._foreground
    canvas.write("""
newpath
%d %d moveto
(%s) true charpath
%g %g %g setrgbcolor
stroke
""" % (x, w-y, text, c[0], c[1], c[2]))

  def draw_open_circle (self, canvas, x, y, radius, color=None) :
    w, h = self.GetSize()
    if (color is None) : color = self._foreground
    canvas.write("""
newpath
%d %d %d 0 360 arc
%g %g %g setrgbcolor
stroke
""" % (x, w-y, radius, color[0], color[1], color[2]))

  def draw_filled_circle (self, canvas, x, y, radius, color) :
    w, h = self.GetSize()
    canvas.write("""
newpath
%d %d %d 0 360 arc
%g %g %g setrgbcolor
fill
""" % (x, w-y, radius, color[0], color[1], color[2]))

if (__name__ == "__main__") :
  run(sys.argv[1:])
