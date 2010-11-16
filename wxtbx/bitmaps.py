
from wx.lib.embeddedimage import PyEmbeddedImage
import wx
import libtbx.load_env
import os

blank16 = PyEmbeddedImage(
    "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAQAAAC1+jfqAAAAAXNSR0IArs4c6QAAAAJiS0dE"
    "AP+Hj8y/AAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH2QsNFBwTBbkZpAAAABl0RVh0"
    "Q29tbWVudABDcmVhdGVkIHdpdGggR0lNUFeBDhcAAAAOSURBVCjPY2AYBaMAAQACEAABFMLA"
    "kgAAAABJRU5ErkJggg==")

icon_lib = libtbx.env.find_in_repositories(
  relative_path=os.path.join("gui_resources", "icons"),
  test=os.path.isdir)

image_cache = {}

def load_png_as_bitmap (icon_path) :
  bmp = image_cache.get(icon_path, None)
  if bmp is None :
    img = wx.Image(icon_path, type=wx.BITMAP_TYPE_PNG, index=-1)
    bmp = img.ConvertToBitmap()
    image_cache[icon_path] = bmp
  return bmp

def find_crystal_icon (icon_class, name, size=32) :
  if icon_lib is not None :
    size_dir = "%dx%d" % (size, size)
    icon_path = os.path.join(icon_lib, "crystal_project", size_dir, icon_class,
      name + ".png")
    if os.path.isfile(icon_path) :
      return icon_path
  return None

def fetch_icon_bitmap (*args, **kwds) :
  icon_path = find_crystal_icon(*args, **kwds)
  if icon_path is None :
    return wx.NullBitmap
  return load_png_as_bitmap(icon_path)

def find_custom_icon (name, ext=".png") :
  if icon_lib is not None :
    icon_path = os.path.join(icon_lib, "custom", name + ext)
    if os.path.isfile(icon_path) :
      return icon_path
  return None

def fetch_custom_icon_bitmap (*args, **kwds) :
  icon_path = find_custom_icon(*args, **kwds)
  if icon_path is None :
    return wx.NullBitmap
  return load_png_as_bitmap(icon_path)
