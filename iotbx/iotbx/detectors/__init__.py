from scitbx.python_utils.misc import import_regular_symbols
from scitbx.array_family import flex
from iotbx_boost import detectors
import_regular_symbols(globals(), detectors.__dict__)
del import_regular_symbols
del detectors
