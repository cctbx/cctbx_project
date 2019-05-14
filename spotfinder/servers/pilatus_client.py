from __future__ import absolute_import, division, print_function
from six.moves import range
import os
from spotfinder.diffraction.imagefiles import quick_image
from spotfinder.servers.multipart_encoder import post_multipart
from libtbx.development.timers import Timer

def get_spotfinder_url(file_object,host,port):
  testurl = "%s:%d"%(host,port)
  selector = "/spotfinder"
  start_index=0
  stop_index = file_object.linearintdata.size()
  raw_string=file_object.linearintdata.slice_to_byte_str(start_index,stop_index)
  query_object = [
    ("moduleindex",file_object.__dict__.get("sliceindex",-1)),
    ("filename",file_object.filename),
    ("bin",1),
    ("vendortype",file_object.vendortype),
    ("beam_center_reference_frame",file_object.beam_center_reference_frame),
    ("beam_center_convention",file_object.beam_center_convention),
    ("header",file_object.header),
    ("headerlines",""),
  ]
  for item in ['DISTANCE', 'PHI', 'WAVELENGTH',
    'TWOTHETA', 'OSC_RANGE',
    'CCD_IMAGE_SATURATION', 'OSC_START', 'DETECTOR_SN', 'PIXEL_SIZE',
    'SIZE1','SIZE2','BEAM_CENTER_X','BEAM_CENTER_Y'
    ]:
    query_object.append((item,file_object.parameters[item]))

  files = [
    ("adsc_data",file_object.filename,raw_string)
  ]

  print("length of data in ints",stop_index)
  print("length of data in bytes",len(raw_string))
  assert len(raw_string)/4==stop_index

  T = Timer("do_POST")
  Response = post_multipart(host=testurl, selector=selector,
    fields = query_object, files = files)
  del T

  print(Response.getresponse().read())

def get_labelit_image_object(file,convention):
  Q = quick_image(file)
  Q.set_beam_center_convention(convention)
  Q.read()

  return Q

def do_main(filepath, force_binning, convention, host, port):
  absfile = os.path.abspath(filepath)
  Q = get_labelit_image_object(absfile, convention)
  if force_binning:
    Q.setBin(2)
    Q.show_header()
  get_spotfinder_url(Q,host,port)

  for x in range(12):
    file = "file://%s?slice=%d"%(absfile,x)
    Q = get_labelit_image_object(file, convention)
    if force_binning:
      Q.setBin(2)
      Q.show_header()
    get_spotfinder_url(Q,host,port)

if __name__=="__main__":
  import sys
  try:
    filepath, force_binning, convention, host, port = sys.argv[1:6]
    force_binning = bool(force_binning)
    port = int(port)
    convention = int(convention)
  except Exception:
    print("""
Usage:
libtbx.python pilatus_client.py <filepath> <force_binning> <convention> <host> <port>
Four mandatory arguments:
  filepath: absolute or relative path name of the Pilatus test image to be analyzed
  force_binning: True (client-side 2x2-pixel software binning; normally
                        not a good choice for Pilatus) or False
  convention: beam_center_convention as defined on the spotfinder servers wiki ==0
  host: usually "localhost"; in any case, must be machine with same endianness
  port: port number of image analyzer http service
""")
  do_main(filepath, force_binning, convention, host, port)
