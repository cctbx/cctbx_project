
"""
Listing of current dependencies for CCTBX and related applications (including
LABELIT, xia2, DIALS, and Phenix with GUI).  Not all of these can be downloaded
via the web (yet).
"""

from __future__ import division

BASE_CCI_PKG_URL = "http://cci.lbl.gov/third_party"
BASE_XIA_PKG_URL = "http://www.ccp4.ac.uk/xia"

# from CCI
PYTHON_PKG = "Python-2.7.5_cci.tar.gz"
# XXX we maintain a patched copy to avoid an ICE with gcc 3.4
NUMPY_PKG = "numpy-1.6.2.tar.gz"         # used many places
IMAGING_PKG = "Imaging-1.1.7.tar.gz"     # for labelit, gltbx
REPORTLAB_PKG = "reportlab-2.6.tar.gz"   # for labelit
ZLIB_PKG = "zlib-1.2.7.tar.gz"
SCIPY_PKG = "scipy-0.11.0.tar.gz"        # not used by default
PYRTF_PKG = "PyRTF-0.45.tar.gz"          # for phenix.table_one, etc.
BIOPYTHON_PKG = "biopython-1.58.tar.gz"  # used in iotbx

# from xia2 page
HDF5_PKG = "hdf5-1.8.8.tar.bz2"      # dxtbx
H5PY_PKG = "h5py-2.0.1-edit.tar.gz"  # dxtbx

# GUI dependencies
LIBPNG_PKG = "libpng-1.2.32.tar.gz"
FREETYPE_PKG = "freetype-2.4.2.tar.gz"
# Linux-only
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
WXPYTHON_DEV_PKG = "wxPython-src-2.9.4.1.tar.gz"  # Mac 64-bit
WXPYTHON_PKG = "wxPython-src-2.8.12.1.tar.gz"     # Linux, Mac 32-bit
WEBKIT_PKG = "wxwebkit.tar.gz"                    # not currently used
MATPLOTLIB_PKG = "matplotlib-1.2.1.tar.gz"
PY2APP_PKG = "py2app-0.7.3.tar.gz"                # Mac only
