"""Floating-overlay geometry helper.

Keeps one overlay widget floated over a host QAbstractScrollArea's
top-right corner, clear of its vertical scrollbar. Nothing here is
search-specific; any floating toolbar over a scroll area can reuse it.
"""

from qttbx.qt import QtCore


class OverlayAnchor(QtCore.QObject):
  """Float ``overlay`` at ``host``'s top-right corner.

  The overlay is sized to about half the host's width -- never
  narrower than its sizeHint, capped to the host for narrow windows;
  the extra width goes to whatever its layout stretches. Repositions
  when the host resizes, when the scroll range changes
  (content growth can summon the vertical scrollbar WITHOUT any host
  Resize, shifting the usable right edge), and on demand via
  ``reposition()`` -- callers show/hide the overlay themselves and call
  that right before ``show()``.

  Parameters
  ----------
  host : QtWidgets.QAbstractScrollArea
      The widget the overlay floats over.
  overlay : QtWidgets.QWidget
      The floating widget (a child of ``host``).
  parent : QtCore.QObject, optional
      Qt parent; defaults to ``overlay``.
  """

  def __init__(self, host, overlay, parent=None):
    super().__init__(parent or overlay)
    self._host = host
    self._overlay = overlay
    host.installEventFilter(self)
    host.verticalScrollBar().rangeChanged.connect(self._on_scroll_range)
    # Lifetime pin: when the host's C++ object goes down, drop our
    # references so anything still queued cannot touch a dying widget.
    # (The installed filter and the scrollbar connection die with the
    # host; this hook covers callbacks already in flight.)
    host.destroyed.connect(self._on_host_destroyed)

  def _on_host_destroyed(self, *_args):
    self._host = None
    self._overlay = None

  def reposition(self):
    """Size and place the overlay now (cheap no-op when nothing moved)."""
    host = getattr(self, "_host", None)
    overlay = getattr(self, "_overlay", None)
    if host is None or overlay is None:
      return
    sb = host.verticalScrollBar()
    sb_w = sb.width() if sb.isVisible() else 0
    hint = overlay.sizeHint()
    width = max(hint.width(), host.width() // 2)
    # The cap subtracts the scrollbar too: a style-scaled scrollbar can
    # be wider than the 16 px margin, and a hint-floored overlay capped
    # only to the host would land on top of it.
    width = min(width, max(0, host.width() - sb_w - 16))
    overlay.resize(width, hint.height())
    x = max(0, host.width() - width - sb_w - 8)
    overlay.move(x, 8)

  # ---- change tracking -----------------------------------------------------

  def eventFilter(self, obj, event):
    # The host <-> anchor reference cycle (Qt parentage + installed
    # event filter) can be reclaimed by cyclic GC at teardown, clearing
    # this object's Python attributes while its C++ half is still
    # installed as the host's filter; a late event must not read the
    # cleared attrs. The destroyed-hook narrows this window but cannot
    # close it (GC can run first), so getattr stays.
    host = getattr(self, "_host", None)
    overlay = getattr(self, "_overlay", None)
    if (host is not None and obj is host
        and event.type() == QtCore.QEvent.Resize):
      try:
        visible = overlay is not None and overlay.isVisible()
      except RuntimeError:              # C++ overlay already deleted
        visible = False
      if visible:
        self.reposition()
    return super().eventFilter(obj, event)

  def _on_scroll_range(self, _min, _max):
    # Same teardown posture as eventFilter: the connected scrollbar can
    # outlive this object's Python attributes (cyclic GC) or the C++
    # widgets can be gone (hard teardown) when a late signal lands.
    overlay = getattr(self, "_overlay", None)
    if overlay is None:
      return
    try:
      visible = overlay.isVisible()
    except RuntimeError:
      return
    if visible:
      # rangeChanged fires BEFORE Qt flips the scrollbar's visibility
      # on the next layout pass; defer so reposition reads the settled
      # state (a synchronous call would still see the old sb width).
      QtCore.QTimer.singleShot(0, self._reposition_if_visible)

  def _reposition_if_visible(self):
    overlay = getattr(self, "_overlay", None)
    if overlay is None:
      return
    try:
      if overlay.isVisible():
        self.reposition()
    except RuntimeError:
      return
