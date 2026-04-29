"""Convert dectector image to png image"""
# LIBTBX_SET_DISPATCHER_NAME labelit.png

from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, Usage
from six import BytesIO
import os
import sys

usage_message ="""\
labelit.png file_name [output_file] [bin=1|2]

  Convert an X-ray detector image to PNG format.

  Input file is a detector image (CBF, ADSC, MAR, etc.).  If the name of the
  output file is not specified, it will be identical to the image file but
  with the extension replaced by '.png'.

  bin argument: 1 ( default, same-sized png  as detector )
                2 ( bin detector pixels 2x2 for display )"""

def run(args):
  if (len(args) == 0) or (len(args) > 3) or "-large" in args or "help" in args:
    raise Usage(usage_message)
  graphics_bin = 1
  if len(args)==3 and args[2].find("bin=")==0 and args[2][4]=="2":
    graphics_bin = 2
  img_file = args[0]
  output_file = None
  if (len(args) >= 2):
    output_file = args[1]
  if (output_file is None):
    output_file = os.path.basename(os.path.splitext(img_file)[0]) + ".png"
  C = convert_image(img_file,graphics_bin)
  C.img.show_header()
  open(output_file, "wb").write(C.output().getvalue())
  print("Wrote %s" % output_file)

class convert_image:

 def __init__(self,file_name,graphics_bin=1):#, output_file=None):
  from rstbx.slip_viewer.slip_viewer_image_factory import SlipViewerImageFactory as ImageFactory
  from iotbx.detectors import ImageException
  try :
    from PIL import Image
  except ImportError :
    import Image
  try :
    img = ImageFactory(file_name)
  except ImageException as e :
    raise Sorry(str(e))
  img.read()
  self.img = img
  if (img.vendortype in ["Pilatus-6M","Pilatus-2M","Pilatus-300K"]):
    graphics_bin = 1
  flex_img = img.get_flex_image(binning=graphics_bin, brightness=1)
  flex_img.setWindow(0.0, 0.0, 1)
  flex_img.spot_convention(0)
  flex_img.adjust(color_scheme=0)
  flex_img.prep_string()
  data_string = flex_img.as_bytes()
  try: # fromstring raises Exception in Pillow >= 3.0.0
    self.imageout = Image.fromstring("RGB",
                         (flex_img.ex_size2(), flex_img.ex_size1()),
                         data_string)
  except Exception:
    self.imageout = Image.frombytes("RGB",
                         (flex_img.ex_size2(), flex_img.ex_size1()),
                         data_string)

 def output(self):
  out = BytesIO()
  self.imageout.save(out, "PNG")
  return out

if (__name__ == "__main__"):
  run(sys.argv[1:])
