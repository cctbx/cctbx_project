import scitbx.array_family.flex

from scitbx.python_utils import misc
scitbx_ext = misc.import_ext("scitbx_boost.array_family.flex_scitbx_ext")
cctbx_ext = misc.import_ext("cctbx_boost.array_family.flex_cctbx_ext")
misc.import_regular_symbols(globals(), scitbx_ext.__dict__)
misc.import_regular_symbols(globals(), cctbx_ext.__dict__)
del misc

to_list = scitbx.array_family.flex.to_list
