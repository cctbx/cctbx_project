
"""
Listing of current dependencies for CCTBX and related applications (including
LABELIT, xia2, DIALS, and Phenix with GUI).  Not all of these can be downloaded
via the web (yet).
"""

from __future__ import division
from installer_utils import *
import urllib2
import os.path as op
import os
import sys

BASE_CCI_PKG_URL = "http://cci.lbl.gov/cctbx_dependencies"
BASE_XIA_PKG_URL = "http://www.ccp4.ac.uk/xia"

# from CCI
PYTHON_PKG = "Python-2.7.8_cci.tar.gz"
NUMPY_PKG = "numpy-1.8.1.tar.gz"         # used many places
IMAGING_PKG = "Imaging-1.1.7.tar.gz"     # for labelit, gltbx
REPORTLAB_PKG = "reportlab-2.6.tar.gz"   # for labelit
ZLIB_PKG = "zlib-1.2.7.tar.gz"
SCIPY_PKG = "scipy-0.14.0.tar.gz"        # not used by default
PYRTF_PKG = "PyRTF-0.45.tar.gz"          # for phenix.table_one, etc.
BIOPYTHON_PKG = "biopython-1.64.tar.gz"  # used in iotbx
SPHINX_PKG = "Sphinx-1.2.2.tar.gz"       # for documentation
NUMPYDOC_PKG = "numpydoc-0.5.tar.gz"     # for documentation
IPYTHON_PKG = "ipython-2.1.0.tar.gz"     # IPython

# from xia2 page
HDF5_PKG = "hdf5-1.8.8.tar.bz2"      # dxtbx
H5PY_PKG = "h5py-2.0.1-edit.tar.gz"  # dxtbx

# GUI dependencies
LIBPNG_PKG = "libpng-1.2.32.tar.gz"
FREETYPE_PKG = "freetype-2.4.2.tar.gz"
# Linux-only
# FIXME some of these are getting pretty ancient, time to update?
GETTEXT_PKG = "gettext-0.18.2.tar.gz"
GLIB_PKG = "glib-2.12.11.tar.gz"
EXPAT_PKG = "expat-1.95.8.tar.gz"
FONTCONFIG_PKG = "fontconfig-2.3.95.tar.gz"
RENDER_PKG = "render-0.8.tar.gz"
XRENDER_PKG = "xrender-0.8.3.tar.gz"
XFT_PKG = "xft-2.1.2.tar.gz"
PIXMAN_PKG = "pixman-0.19.2.tar.gz"
CAIRO_PKG = "cairo-1.8.10.tar.gz"
PANGO_PKG = "pango-1.16.1.tar.gz"
ATK_PKG = "atk-1.9.1.tar.gz"
TIFF_PKG = "tiff-v3.6.1.tar.gz"
GTK_PKG = "gtk+-2.10.11.tar.gz"
GTK_ENGINE_PKG = "clearlooks-0.5.tar.gz"
GTK_THEME_PKG = "gtk_themes.tar.gz"
# end Linux-only
FONT_PKG = "fonts.tar.gz"
# FIXME at some point we should switch to using 3.x for all platforms
WXPYTHON_DEV_PKG = "wxPython-src-3.0.1.0.tar.gz"  # Mac 64-bit
WXPYTHON_PKG = "wxPython-src-2.8.12.1.tar.gz"         # Linux, Mac 32-bit
MATPLOTLIB_PKG = "matplotlib-1.3.1.tar.gz"
PY2APP_PKG = "py2app-0.7.3.tar.gz"                    # Mac only
PYOPENGL_PKG = "PyOpenGL-3.1.0.tar.gz"
# https://pypi.python.org/pypi/Send2Trash
SEND2TRASH_PKG = "Send2Trash-1.3.0.tar.gz"

# Various dependencies from external repositories, distributed as static
# tarballs (since they are not under active development by us or our
# collaborators)
dependency_tarballs = {
  "boost":  ("http://cci.lbl.gov/hot", "boost_hot.tar.gz"),
  "scons":  ("http://cci.lbl.gov/hot", "scons_hot.tar.gz"),
  "annlib": ("http://cci.lbl.gov/hot", "annlib_hot.tar.gz"),
}
# External SVN repositories that may be required for certain components of
# CCTBX to work.  This includes forked versions (with minimal changes) of the
# core CCP4 libraries, MUSCLE, and ksDSSP, but also the development branch
# of CBFLIB.
subversion_repositories = {
  "cbflib":"http://svn.code.sf.net/p/cbflib/code-0/trunk/CBFlib_bleeding_edge",
  "ccp4io": "http://cci.lbl.gov/svn/ccp4io/trunk",
  "ccp4io_adaptbx": "http://cci.lbl.gov/svn/ccp4io_adaptbx/trunk",
  "annlib_adaptbx": "http://cci.lbl.gov/svn/annlib_adaptbx/trunk",
  "gui_resources": "http://cci.lbl.gov/svn/gui_resources/trunk",
  "tntbx": "http://cci.lbl.gov/svn/tntbx/trunk",
  "ksdssp": "http://cci.lbl.gov/svn/ksdssp/trunk",
  "muscle": "http://cci.lbl.gov/svn/muscle/trunk",
}

class fetch_packages (object) :
  """
  Download manager for the packages defined by this module - this is used by
  install_base_packages.py but also for setting up installer bundles.
  """
  def __init__ (self, dest_dir, log, pkg_dirs=None, no_download=False,
      copy_files=False) :
    self.dest_dir = dest_dir
    self.log = log
    self.pkg_dirs = pkg_dirs
    self.no_download = no_download
    self.copy_files = copy_files

  def __call__ (self, pkg_name, pkg_url=None, output_file=None) :
    if (pkg_url is None) :
      pkg_url = BASE_CCI_PKG_URL
    if (output_file is None) :
      output_file = pkg_name
    os.chdir(self.dest_dir)
    print >> self.log, "  getting package %s..." % pkg_name
    if (self.pkg_dirs is not None) and (len(self.pkg_dirs) > 0) :
      for pkg_dir in self.pkg_dirs :
        static_file = op.join(pkg_dir, pkg_name)
        if (op.exists(static_file)) :
          print >> self.log, "    using %s" % static_file
          if self.copy_files :
            copy_file(static_file, op.join(self.dest_dir, output_file))
            return op.join(self.dest_dir, output_file)
          else :
            return static_file
    if (op.exists(pkg_name)) :
      print >> self.log, "    using ./%s" % pkg_name
      return op.join(self.dest_dir, pkg_name)
    else :
      if (self.no_download) :
        raise RuntimeError(("Package '%s' not found on local filesystems.  ") %
          pkg_name)
      full_url = "%s/%s" % (pkg_url, pkg_name)
      self.log.write("    downloading from %s : " % pkg_url)
      f = open(output_file, "wb")
      data = urllib2.urlopen(full_url).read()
      assert (len(data) > 0), pkg_name
      self.log.write("%d KB\n" % (len(data) / 1024))
      self.log.flush()
      f.write(data)
      f.close()
      return op.join(self.dest_dir, output_file)

def fetch_all_dependencies (dest_dir,
    log,
    pkg_dirs=None,
    copy_files=True,
    gui_packages=True,
    dials_packages=True) :
  """
  Download or copy all dependencies into a local directory (prepping for
  source installer bundling).
  """
  fetch_package = fetch_packages(
    dest_dir=dest_dir,
    log=log,
    pkg_dirs=pkg_dirs,
    copy_files=copy_files)
  for pkg_name in [
      PYTHON_PKG, NUMPY_PKG, IMAGING_PKG, REPORTLAB_PKG, ZLIB_PKG,
      SCIPY_PKG, PYRTF_PKG, BIOPYTHON_PKG, SPHINX_PKG, NUMPYDOC_PKG,
      IPYTHON_PKG,
    ] :
    fetch_package(pkg_name)
  if (gui_packages) :
    for pkg_name in [
        LIBPNG_PKG, FREETYPE_PKG, GETTEXT_PKG, GLIB_PKG, EXPAT_PKG,
        FONTCONFIG_PKG, RENDER_PKG, XRENDER_PKG, XFT_PKG, PIXMAN_PKG,
        CAIRO_PKG, PANGO_PKG, ATK_PKG, TIFF_PKG, GTK_PKG,
        GTK_ENGINE_PKG, GTK_THEME_PKG, FONT_PKG, WXPYTHON_DEV_PKG, WXPYTHON_PKG,
        MATPLOTLIB_PKG, PY2APP_PKG, SEND2TRASH_PKG,
      ] :
      fetch_package(pkg_name)
  if (dials_packages) :
    for pkg_name in [ HDF5_PKG, H5PY_PKG, PYOPENGL_PKG ] :
      fetch_package(pkg_name, BASE_XIA_PKG_URL)

def fetch_repository (pkg_name, pkg_url=None, working_copy=True,
    delete_if_present=False) :
  """
  Download an SVN repository, with or without metadata required for ongoing
  development.
  """
  if op.exists(pkg_name) :
    if delete_if_present :
      shutil.rmtree(pkg_name)
    else :
      raise OSError("Directory '%s' already exists.")
  if (pkg_url is None) :
    pkg_url = optional_repositories[pkg_name]
  if working_copy :
    call("svn co %s %s" % (pkg_url, pkg_name), sys.stdout)
  else :
    call("svn export %s %s" % (pkg_url, pkg_name), sys.stdout)
  assert op.isdir(pkg_name)

def fetch_remote_package (module_name, log=sys.stdout, working_copy=False) :
  if op.isdir(module_name) :
    shutil.rmtree(module_name)
  if (module_name in dependency_tarballs) :
    pkg_url, pkg_name = dependency_tarballs[module_name]
    tarfile = module_name + ".tar.gz"
    fetch_packages(
      dest_dir=os.getcwd(),
      log=log).__call__(
        pkg_name=pkg_name,
        pkg_url=pkg_url,
        output_file=tarfile)
    untar(tarfile, log)
    os.remove(tarfile)
  elif (module_name in subversion_repositories) :
    pkg_url = subversion_repositories[module_name]
    fetch_repository(
      pkg_name=module_name,
      pkg_url=pkg_url,
      working_copy=working_copy)
