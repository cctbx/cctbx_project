"""Palette-aware color helpers for validation feedback.

Every color used for validation feedback (invalid background fill, invalid
border, error-icon foreground, low-emphasis side labels) derives from the
active ``QApplication.palette()`` so the widgets adapt to light and dark
themes without any per-widget configuration.

Helpers return ``QColor`` objects; consumers convert to stylesheet strings
or ``QBrush`` as needed.
"""

from PySide2.QtGui import QColor, QPalette


def _active_palette():
  """Return the active QApplication palette, or a default if none exists.

  Returns
  -------
    QPalette
  """
  from PySide2.QtWidgets import QApplication
  app = QApplication.instance()
  if app is not None:
    return app.palette()
  return QPalette()


def is_dark_palette(palette=None):
  """Return True iff ``palette`` (or the active one) is a dark theme.

  The decision is based on the lightness of ``QPalette.Base`` (the
  body color of input widgets): values < 128 indicate a dark theme.

  Parameters
  ----------
    palette: QPalette, optional
      Defaults to the active QApplication palette.

  Returns
  -------
    bool
  """
  if palette is None:
    palette = _active_palette()
  return palette.color(QPalette.Base).lightness() < 128


def invalid_background(palette=None):
  """Return the background ``QColor`` for an invalid input cell or field.

  Light palettes: pale pink (``#fff5f5``).
  Dark palettes:  a desaturated dark red.

  Parameters
  ----------
    palette: QPalette, optional

  Returns
  -------
    QColor
  """
  if is_dark_palette(palette):
    return QColor(64, 28, 28)
  return QColor("#fff5f5")


def invalid_border(palette=None):
  """Return the border ``QColor`` for an invalid input field.

  Parameters
  ----------
    palette: QPalette, optional

  Returns
  -------
    QColor
  """
  if is_dark_palette(palette):
    return QColor(255, 100, 100)
  return QColor(220, 50, 50)


def error_emphasis(palette=None):
  """Return the foreground ``QColor`` for warning icons and strong error text.

  Parameters
  ----------
    palette: QPalette, optional

  Returns
  -------
    QColor
  """
  if is_dark_palette(palette):
    return QColor(255, 130, 130)
  return QColor("#c00")


def secondary_label(palette=None):
  """Return the foreground ``QColor`` for low-emphasis side labels.

  A muted gray suitable for ``SpaceGroupWidget``'s Hall/HM hint, picked
  to remain readable on both light and dark backgrounds.

  Parameters
  ----------
    palette: QPalette, optional

  Returns
  -------
    QColor
  """
  if is_dark_palette(palette):
    return QColor(170, 170, 170)
  return QColor(85, 85, 85)
