
from wx.lib.embeddedimage import PyEmbeddedImage
import wx
import libtbx.load_env
import sys, os

blank16 = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAQAAAC1+jfqAAAAAXNSR0IArs4c6QAAAAJiS0dE"
    "AP+Hj8y/AAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH2QsNFBwTBbkZpAAAABl0RVh0"
    "Q29tbWVudABDcmVhdGVkIHdpdGggR0lNUFeBDhcAAAAOSURBVCjPY2AYBaMAAQACEAABFMLA"
    "kgAAAABJRU5ErkJggg==")

icon_lib = libtbx.env.find_in_repositories(
  relative_path="gui_resources/icons",
  test=os.path.isdir)

def find_crystal_icon (icon_class, name, size=32) :
  if icon_lib is not None :
    size_dir = "%dx%d" % (size, size)
    icon_path = os.path.join(icon_lib, "crystal_project", size_dir, icon_class,
      name + ".png")
    if os.path.isfile(icon_path) :
      return icon_path
  return None

image_cache = {}

def fetch_icon_bitmap (*args, **kwds) :
  icon_path = find_crystal_icon(*args, **kwds)
  bmp = image_cache.get(icon_path, None)
  if bmp is None :
    img = wx.Image(icon_path, type=wx.BITMAP_TYPE_PNG, index=-1)
    bmp = img.ConvertToBitmap()
    image_cache[icon_path] = bmp
  return bmp
