from __future__ import absolute_import, division, print_function
import os.path
from libtbx.utils import to_bytes
from six.moves import range

def encode(data):
  edata = ""
  for i in range(len(data)):
    edata += "%.2x" % ord(data[i])
  return edata

def create_encoded(image_file_name):
  import wx
  img = wx.Image(name=image_file_name)
  w,h = img.GetSize()
  name, ext = os.path.splitext(os.path.basename(image_file_name))
  print (
    '%s_img = img_data(width=%d, height=%d, mask=-1, encoded_data = """\\'
     % (name, w, h))
  encoded = encode(img.GetData())
  while (len(encoded) > 0):
    print(encoded[:78]+"\\")
    encoded = encoded[78:]
  print('""")')
  print()

class img_data:

  def __init__(self, width, height, mask, encoded_data):
    self.width = width
    self.height = height
    self.data = self.decode(encoded_data)
    self.mask = mask * 3

  def get_width(self): return self.width
  def get_height(self): return self.height
  def get_size(self): return (self.width, self.height)
  def get_data(self): return self.data
  def get_mask(self): return self.mask

  def decode(self, edata):
    hex_chars = {"0":  0, "1":  1, "2":  2, "3":  3,
                 "4":  4, "5":  5, "6":  6, "7":  7,
                 "8":  8, "9":  9, "a": 10, "b": 11,
                 "c": 12, "d": 13, "e": 14, "f": 15}
    data = ""
    for i in range(0, len(edata), 2):
      data += chr(hex_chars[edata[i]] * 16 + hex_chars[edata[i+1]])
    return data

  def as_wx_Bitmap(self):
    import wx
    w,h = self.get_size()
    data = to_bytes(self.get_data())
    mask = self.get_mask()
    img = wx.Image(w, h)
    img.SetData(data)
    if (mask >= 0):
      img.SetMaskColour(ord(data[mask]), ord(data[mask+1]), ord(data[mask+2]))
      img.SetMask()
    return img.ConvertToBitmap()

if (__name__ == "__main__"):
  import sys
  print('from gltbx.images import img_data\n')
  for arg in sys.argv[1:]:
    create_encoded(arg)
