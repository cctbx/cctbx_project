from __future__ import absolute_import, division, print_function

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

def load_png_as_bitmap(icon_path, scale=None):
  bmp = image_cache.get(icon_path, None)
  if bmp is None :
    if (wx.Platform == '__WXMSW__') and (wx.VERSION >= (2,9)):
# libpng purists have now broken backwards compatibility when loading old images.
# This leads to unhelpful error message boxes. Suppress these temporarily with LogNull
      noLog = wx.LogNull()
      img = wx.Image(icon_path, type=wx.BITMAP_TYPE_PNG, index=-1)
    else:
      img = wx.Image(icon_path, type=wx.BITMAP_TYPE_PNG, index=-1)
    if (scale is not None):
      assert isinstance(scale, tuple)
      w, h = scale
      img = img.Scale(w, h, wx.IMAGE_QUALITY_NORMAL)
    bmp = img.ConvertToBitmap()
    image_cache[icon_path] = bmp
  return bmp

def find_crystal_icon(icon_class, name, size=32, scale=None):
  if icon_lib is not None :
    size_dir = "%dx%d" % (size, size)
    icon_path = os.path.join(icon_lib, "crystal_project", size_dir, icon_class,
      name + ".png")
    if os.path.isfile(icon_path):
      return icon_path
  return None

def fetch_icon_bitmap(*args, **kwds):
  icon_path = find_crystal_icon(*args, **kwds)
  if icon_path is None :
    return wx.NullBitmap
  icon_size = kwds.get("icon_size", None)
  if icon_size:
    scale = (icon_size, icon_size)
  else:
    scale = kwds.get("scale", None)
  return load_png_as_bitmap(icon_path, scale=scale)

def find_custom_icon(name, ext=".png", scale=None, icon_size=None):
  if icon_lib is not None :
    icon_path = os.path.join(icon_lib, "custom", name + ext)
    if os.path.isfile(icon_path):
      return icon_path
  return None

def fetch_custom_icon_bitmap(*args, **kwds):
  icon_path = find_custom_icon(*args, **kwds)
  if icon_path is None :
    return wx.NullBitmap
  icon_size = kwds.get("icon_size", None)
  if icon_size:
    scale = (icon_size, icon_size)
  else:
    scale = kwds.get("scale", None)
  return load_png_as_bitmap(icon_path, scale=scale)
