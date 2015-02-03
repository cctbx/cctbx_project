# -*- python -*-
from __future__ import division
import os
import sys
import subprocess
import optparse
import getpass
import urllib
import urlparse

# A reworked base dependency installer.
# This is a work in progress!

class Package(object):
  def __init__(self, python_exe=None):
    self.python_exe = python_exe

  def fetch(self):
    pass

  def unpack(self):
    pass

  def build(self):
    pass

  def install(self):
    pass

##### Python Packages #####

class PythonPackage(Package):
  def build(self):
    pass
    # self.python_exe setup.py build
  def install(self):
    pass
    # self.python_exe setup.py install

class package_numpy(PythonPackage):
  package = "numpy"
  source = "numpy-1.8.1.tar.gz"
  
class package_docutils(PythonPackage):
  package = "docutils"
  source = "docutils-0.12.tar.gz"

class package_biopython(PythonPackage):
  package = "biopython"
  source = "biopython-1.64.tar.gz"
  
class package_imaging(PythonPackage):
  package = "imaging"
  source = "Imaging-1.1.7.tar.gz"
  
class package_scipy(PythonPackage):
  package = "scipy"
  source = "scipy-0.14.0.tar.gz"

class package_py2app(PythonPackage):
  package = "py2app"
  source = "py2app-0.7.3.tar.gz"

class package_reportlab(PythonPackage):
  package = "reportlab"
  source = "reportlab-2.6.tar.gz"
  
class package_sphinx(PythonPackage):
  package = "sphinx"
  pymodule = "sphinx"
  dependencies = ['numpydoc']
  source = "Sphinx-1.2.2.tar.gz"
  
class package_numpydoc(PythonPackage):
  package = "numpydoc"
  pymodule = "numpydoc"
  source = "numpydoc-0.5.tar.gz"
  
class package_ipython(PythonPackage):
  package = "ipython"
  pymodule = "IPython"
  source = "ipython-2.1.0.tar.gz"
  
class package_pyopengl(PythonPackage):
  package = "pyopengl"
  pymodule = "OpenGL"
  source = "PyOpenGL-3.1.0.tar.gz"
  
  
class package_matplotlib(PythonPackage):
  package = "matplotlib"
  source = "matplotlib-1.3.1.tar.gz"

class package_pyrtf(PythonPackage):
  package = "pyrtf"
  pymodule = "PyRTF"
  source = "PyRTF-0.45.tar.gz"
  
class package_send2trash(PythonPackage):
  package = "send2trash"
  pymodule = "send2trash"
  source = "Send2Trash-1.3.0.tar.gz"
  
class package_h5py(PythonPackage):
  package = "h5py"
  pymodule = "h5py"
  source = "h5py-2.0.1-edit.tar.gz"

##### Autoconf-style Packages #####

class AutoconfPackage(Package):
  pass

class package_hdf5(AutoconfPackage):
  package = "hdf5"
  source = "hdf5-1.8.8.tar.bz2"
  
class package_freetype(AutoconfPackage):
  package = "freetype"
  source = "freetype-2.4.2.tar.gz"
  
class package_png(AutoconfPackage):
  package = "png"
  source = "libpng-1.2.52.tar.gz"
  
class package_gettext(AutoconfPackage):
  package = "gettext"
  source = "gettext-0.18.2.tar.gz"
  
class package_glib(AutoconfPackage):
  package = "glib"
  source = "glib-2.12.11.tar.gz"
  
class package_expat(AutoconfPackage):
  package = "expat"
  source = "expat-1.95.8.tar.gz"
  
class package_fontconfig(AutoconfPackage):
  package = "fontconfig"
  source = "fontconfig-2.3.95.tar.gz"
  
class package_render(AutoconfPackage):
  package = "render"
  source = "render-0.8.tar.gz"
  
class package_Xrender(AutoconfPackage):
  package = "xrender"
  source = "xrender-0.8.3.tar.gz"
  
class package_xft(AutoconfPackage):
  package = "xft"
  source = "xft-2.1.2.tar.gz"
  
class package_pixman(AutoconfPackage):
  package = "pixman"
  source = "pixman-0.19.2.tar.gz"
  
class package_cairo(AutoconfPackage):
  package = "cairo"
  source = "cairo-1.8.10.tar.gz"
  
class package_pango(AutoconfPackage):
  package = "pango"
  source = "pango-1.16.1.tar.gz"
  
class package_atk(AutoconfPackage):
  package = "atk"
  source = "atk-1.9.1.tar.gz"

class package_tiff(AutoconfPackage):
  package = "tiff"
  source = "tiff-v3.6.1.tar.gz"
  
class package_gtk(AutoconfPackage):
  package = "gtk"
  source = "gtk+-2.10.11.tar.gz"
  
class package_gtkengine(AutoconfPackage):
  package = "gtkengine"
  source = "clearlooks-0.5.tar.gz"

##### Other packages #####
  
class package_fonts(Package):
  package = "fonts"
  source = "fonts.tar.gz"

class package_wxpython3(Package):
  package = "wxpython3"
  source = "wxPython-src-3.0.1.0.tar.gz"

class package_wxpython2(Package):
  package = "wxpython2"
  source = "wxPython-src-2.8.12.1.tar.gz"
