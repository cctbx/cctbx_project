from scitbx.python_utils import misc
ext = misc.import_ext("scitbx_boost.rational_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc
