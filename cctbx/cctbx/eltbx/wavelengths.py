from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.eltbx.wavelengths_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc
