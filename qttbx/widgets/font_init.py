"""Suppress Qt's 'Sans Serif alias' startup warning.

On platforms where the "Sans Serif" font family isn't a real concrete
family (notably the ``offscreen`` QPA plugin used by the chat tests,
but also some bare Linux installs), Qt's font system has to populate
an alias map the first time a widget asks for a font. It does this
exactly once per process, caches the result, and prints

    qt.qpa.fonts: Populating font family aliases took N ms.
    Replace uses of missing font family "Sans Serif" with one that
    exists to avoid this cost.

This is purely informational; the alias scan happens regardless of
whether we suppress the message. There is no way to set the default
font early enough to avoid the scan (the typical "set
``QApplication.setFont`` to a real family" advice doesn't work --
``QFontDatabase.systemFont(GeneralFont)`` returns the same alias on
the affected platforms, and even probing for a real family via
``QFontInfo`` triggers the scan). The only reliable suppression is a
Qt message handler installed BEFORE ``QApplication`` construction.

Call :func:`init_default_app_font` once at the very top of any entry
point that builds a ``QApplication`` (launcher, test fixture).
Idempotent: subsequent calls reinstall the same handler.
"""

import sys

from qttbx.qt import QtCore


_FILTERED_NEEDLE = "Populating font family aliases"


def _msg_handler(msg_type, ctx, msg):
  # Drop the one-shot alias-population notice; forward everything else
  # to stderr the way Qt's default handler does.
  if _FILTERED_NEEDLE in msg:
    return
  sys.stderr.write(msg + "\n")


def init_default_app_font(app=None):
  """Install a Qt message handler that drops the 'Sans Serif alias'
  startup notice.

  Parameters
  ----------
  app : QtWidgets.QApplication, optional
      Accepted for call-site symmetry with the typical
      ``init_X(app)`` pattern; unused -- Qt's message handler is
      process-wide.

  Returns
  -------
  QtWidgets.QApplication or None
      The ``app`` argument unchanged (so callers can chain
      ``init_default_app_font(QtWidgets.QApplication(...))``).
  """
  QtCore.qInstallMessageHandler(_msg_handler)
  return app
