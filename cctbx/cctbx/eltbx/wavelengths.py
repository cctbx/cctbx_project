from scitbx.python_utils.misc import import_regular_symbols
from cctbx_boost.eltbx import wavelengths_ext as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols
