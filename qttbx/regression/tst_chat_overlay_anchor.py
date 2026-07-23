"""OverlayAnchor tests: half-width sizing with the sizeHint floor and
narrow-host cap, scrollbar clearance with the deferred reposition on
rangeChanged, no auto-reposition while the overlay is hidden, and the
destroyed-hook teardown.

The host is a plain QScrollArea -- the anchor must work without any
chat widget, which is the point of the extraction."""

import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times

try:
  from qttbx.qt import QtCore, QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)


def _app():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  return app


def _pump(app, ms=100):
  timer = QtCore.QElapsedTimer()
  timer.start()
  while timer.elapsed() < ms:
    app.processEvents(QtCore.QEventLoop.AllEvents, 20)


class _Overlay(QtWidgets.QFrame):
  """Overlay with a fixed, known sizeHint."""

  def sizeHint(self):
    return QtCore.QSize(120, 24)


def _anchored_host(app, width=400, height=300, content_height=10):
  """Shown QScrollArea + overlay + anchor; content height is settable."""
  from qttbx.widgets.chat.overlay_anchor import OverlayAnchor
  host = QtWidgets.QScrollArea()
  host.setWidgetResizable(True)
  inner = QtWidgets.QWidget()
  inner.setMinimumHeight(content_height)
  host.setWidget(inner)
  host.resize(width, height)
  host.show()
  overlay = _Overlay(host)
  anchor = OverlayAnchor(host, overlay)
  _pump(app, 30)
  return host, inner, overlay, anchor


def exercise_half_width_floor_and_cap():
  """Sizing: half the host width; never narrower than the overlay's
  sizeHint; capped to the host minus margins for narrow hosts. The
  host Resize repositions a VISIBLE overlay automatically."""
  app = _app()
  host, _inner, overlay, anchor = _anchored_host(app, width=400)
  anchor.reposition()
  overlay.show()
  _pump(app, 30)
  # Half of 400 beats the 120 hint: width 200, floated at y=8 with an
  # 8 px right margin (no scrollbar for a 10 px content).
  assert not host.verticalScrollBar().isVisible()
  assert overlay.width() == 200, overlay.width()
  assert overlay.y() == 8, overlay.y()
  assert overlay.x() == 400 - 200 - 8, overlay.x()
  # Shrink: half (100) drops under the 120 hint -> floor at the hint.
  host.resize(200, 300)
  _pump(app, 30)
  assert overlay.width() == 120, overlay.width()
  # Very narrow: the host cap (width - 16) beats even the hint.
  host.resize(100, 300)
  _pump(app, 30)
  assert overlay.width() == 100 - 16, overlay.width()
  assert overlay.x() >= 0


def exercise_scrollbar_clearance_deferred_on_range_change():
  """Content growth summons the vertical scrollbar WITHOUT any host
  Resize; the deferred rangeChanged reposition must shift the overlay
  clear of it."""
  app = _app()
  host, inner, overlay, anchor = _anchored_host(app, width=400)
  anchor.reposition()
  overlay.show()
  _pump(app, 30)
  x_without_sb = overlay.x()
  sb = host.verticalScrollBar()
  assert not sb.isVisible()
  inner.setMinimumHeight(2000)                  # summon the scrollbar
  _pump(app)
  assert sb.isVisible()
  assert overlay.x() < x_without_sb, (overlay.x(), x_without_sb)
  assert overlay.x() + overlay.width() <= host.width() - sb.width()


def exercise_hidden_overlay_is_left_alone():
  """Host resizes and range changes must not move a HIDDEN overlay
  (the search bar parks closed over the view); an explicit
  reposition() still places it -- open() calls that before show()."""
  app = _app()
  host, inner, overlay, _anchor = _anchored_host(app, width=400)
  overlay.move(3, 3)
  host.resize(500, 300)
  inner.setMinimumHeight(2000)
  _pump(app)
  assert not overlay.isVisible()
  assert (overlay.x(), overlay.y()) == (3, 3), (overlay.x(), overlay.y())
  _anchor.reposition()
  assert overlay.y() == 8
  assert overlay.width() == 250, overlay.width()


def exercise_wide_scrollbar_narrow_host_cap():
  """A style-scaled scrollbar wider than the 16 px margin (fontless
  Windows CI scales everything up): the width cap must subtract the
  scrollbar too, or a hint-floored overlay lands on top of it. The
  stylesheet is applied before show so the scrollbar is wide from the
  first layout pass."""
  app = _app()
  from qttbx.widgets.chat.overlay_anchor import OverlayAnchor
  host = QtWidgets.QScrollArea()
  host.setStyleSheet("QScrollBar:vertical { width: 30px; }")
  host.setWidgetResizable(True)
  inner = QtWidgets.QWidget()
  inner.setMinimumHeight(2000)
  host.setWidget(inner)
  host.resize(140, 300)
  host.show()
  overlay = _Overlay(host)
  anchor = OverlayAnchor(host, overlay)
  _pump(app, 30)
  anchor.reposition()
  overlay.show()
  _pump(app, 30)
  sb = host.verticalScrollBar()
  assert sb.isVisible() and sb.width() == 30, sb.width()
  assert overlay.width() == 140 - 30 - 16, overlay.width()
  assert overlay.x() + overlay.width() <= host.width() - sb.width(), \
      (overlay.x(), overlay.width(), host.width(), sb.width())


def exercise_destroyed_host_drops_references():
  """Killing the host trips the destroyed-hook: the anchor drops its
  references, and late calls (a queued reposition, a stale rangeChanged
  slot) are silent no-ops instead of touching dying widgets.

  shiboken.delete destroys the C++ host synchronously -- deleteLater
  would never fire here, since DeferredDelete events are only handled
  by a running exec loop, not by processEvents."""
  from qttbx.qt import shiboken
  app = _app()
  host, _inner, overlay, anchor = _anchored_host(app)
  assert anchor._host is host and anchor._overlay is overlay
  shiboken.delete(host)                         # overlay dies with it
  assert anchor._host is None, anchor._host
  assert anchor._overlay is None, anchor._overlay
  anchor.reposition()                           # all must no-op quietly
  anchor._on_scroll_range(0, 10)
  anchor._reposition_if_visible()
  _pump(app, 20)


def exercise():
  exercise_half_width_floor_and_cap()
  exercise_scrollbar_clearance_deferred_on_range_change()
  exercise_hidden_overlay_is_left_alone()
  exercise_wide_scrollbar_narrow_host_cap()
  exercise_destroyed_host_drops_references()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
