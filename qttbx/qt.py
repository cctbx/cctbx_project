"""Qt binding shim.

Imports PySide6 if available, otherwise falls back to PySide2. All qttbx code
should import Qt namespaces through this module so it works against either
binding:

    from qttbx.qt import QtCore, QtGui, QtWidgets

The submodule import form is also supported, mirroring PySide directly:

    from qttbx.qt.QtCore import Qt, QSignalBlocker
    from qttbx.qt.QtWidgets import QApplication

QT_API is the string "PySide6" or "PySide2" — useful when behavior must
diverge between the two bindings.
"""

import sys

try:
  from PySide6 import QtCore, QtGui, QtWidgets
  import shiboken6 as shiboken
  QT_API = "PySide6"
except ImportError:
  try:
    from PySide2 import QtCore, QtGui, QtWidgets
    import shiboken2 as shiboken
    QT_API = "PySide2"
  except ImportError as _err:
    raise ImportError(
      "qttbx requires a Qt for Python binding, but neither PySide6 nor "
      "PySide2 could be imported. Install one with:\n"
      "  conda install -c conda-forge pyside6   (recommended)\n"
      "  conda install -c conda-forge pyside2\n"
      "or `pip install PySide6` / `pip install PySide2`."
    ) from _err

# Cross-version aliases for symbols that Qt6 moved out of QtWidgets into
# QtGui. Branch on QT_API rather than try/except: it ties each alias
# directly to the binding it targets, and the cost of the missing-attr
# probe is avoided.
if QT_API == "PySide6":
  QShortcut = QtGui.QShortcut
else:
  QShortcut = QtWidgets.QShortcut

__all__ = ["QtCore", "QtGui", "QtWidgets", "shiboken", "QT_API", "QShortcut"]

# Register the namespace modules under qttbx.qt.* so callers can use the
# submodule import form (e.g. `from qttbx.qt.QtCore import Qt`). PySide
# exposes QtCore/QtGui/QtWidgets as real submodules; we re-bind them here.
sys.modules[__name__ + ".QtCore"] = QtCore
sys.modules[__name__ + ".QtGui"] = QtGui
sys.modules[__name__ + ".QtWidgets"] = QtWidgets
