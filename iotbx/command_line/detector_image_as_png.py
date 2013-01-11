# LIBTBX_SET_DISPATCHER_NAME labelit.png

from __future__ import division
from libtbx.utils import Sorry, Usage
from cStringIO import StringIO
import os
import sys

def run (args) :
  if (len(args) == 0) or (len(args) > 2) :
    raise Usage("""\
labelit.png file_name [output_file]

  Convert an X-ray detector image to PNG format.

  Input file is a detector image (CBF, ADSC, MAR, etc.).  If the name of the
  output file is not specified, it will be identical to the image file but
  with the extension replaced by '.png'.""")
  img_file = args[0]
  output_file = None
  if (len(args) == 2) :
    output_file = args[1]
  convert_image(img_file, output_file)

def convert_image (file_name, output_file=None) :
  from iotbx.detectors import ImageFactory, ImageException
  import Image
  if (output_file is None) :
    output_file = os.path.basename(os.path.splitext(file_name)[0]) + ".png"
  try :
    img = ImageFactory(file_name)
  except ImageException, e :
    raise Sorry(str(e))
  img.read()
  print img.show_header()
  binning = 2
  if (img.vendortype in ["Pilatus-6M","Pilatus-2M","Pilatus-300K"]) :
    binning = 1
  flex_img = img.get_flex_image(binning=binning, brightness=1)
  flex_img.setWindow(0.0, 0.0, 1)
  flex_img.spot_convention(0)
  flex_img.adjust(color_scheme=0)
  flex_img.prep_string()
  data_string = flex_img.export_string
  imageout = Image.fromstring("RGB",
                     (flex_img.ex_size2(), flex_img.ex_size1()),
                     data_string)
  out = StringIO()
  imageout.save(out, "PNG")
  open(output_file, "wb").write(out.getvalue())
  print "Wrote %s" % output_file

if (__name__ == "__main__") :
  run(sys.argv[1:])
