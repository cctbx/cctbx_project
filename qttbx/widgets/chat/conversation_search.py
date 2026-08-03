"""Find-in-conversation controller.

Binds one SearchBar to one ConversationView: scans the view's
searchable cells for a literal, case-insensitive query; paints every
match as a yellow extra selection (current match orange); keeps the
bar's count label fresh. Overlay geometry (floating the bar at the
view's top-right) is OverlayAnchor's job.

Highlights are a paint-time overlay (``setExtraSelections``) -- the
documents are never modified. The one actively streaming markdown cell
re-sets its whole document per delta, which drops that cell's
selections and re-flows its plain-text offsets; re-scans therefore
happen on each interaction (query edit, scope toggle, navigation), not
per delta.
"""

import re

from qttbx.qt import QtCore, QtGui, QtWidgets

# All-matches / current-match highlight colors. Explicit background AND
# foreground so contrast holds on light and dark themes.
_MATCH_BG = "#FFF176"
_CURRENT_BG = "#FF9800"
_MATCH_FG = "#000000"

# Scan cap: painting thousands of extra selections (a one-letter query
# against 50K of tool output) stalls the UI. The label shows "n/500+".
MAX_MATCHES = 500


class ConversationSearch(QtCore.QObject):
  """Search lifecycle + scanning + highlighting for one view.

  Parameters
  ----------
  view : ConversationView
      The conversation to search. The SearchBar is created as a child
      of this view and floated over its top-right corner.
  parent : QtCore.QObject, optional
      Qt parent; defaults to ``view``.
  """

  closed = QtCore.Signal()

  def __init__(self, view, parent=None):
    super().__init__(parent or view)
    from qttbx.widgets.chat.overlay_anchor import OverlayAnchor
    from qttbx.widgets.chat.search_bar import SearchBar
    self._view = view
    self.bar = SearchBar(view)
    self.bar.hide()
    self._anchor = OverlayAnchor(view, self.bar)
    self._matches = []        # (widget, start, length), display order
    self._current = -1        # index into _matches; -1 = none
    self._capped = False
    self._highlighted = []    # widgets currently carrying selections
    self.bar.query_changed.connect(self._on_query_changed)
    self.bar.scope_changed.connect(self._on_scope_changed)
    self.bar.close_requested.connect(self.close)
    self.bar.next_requested.connect(lambda: self._go(+1))
    self.bar.prev_requested.connect(lambda: self._go(-1))

  # ---- lifecycle -----------------------------------------------------------

  def is_open(self):
    return self.bar.isVisible()

  def toggle(self):
    if self.is_open():
      self.close()
    else:
      self.open()

  def open(self):
    self._anchor.reposition()
    self.bar.show()
    self.bar.raise_()
    self.bar.focus_query()
    self._rescan(reset_current=True)
    self._apply()

  def close(self):
    if not self.is_open():
      return
    self._clear_highlights()
    self._matches = []
    self._current = -1
    # Park focus on the view BEFORE hiding the bar. Hiding the focus
    # widget makes Qt walk the tab chain (focusNextPrevChild), and
    # QScrollArea scrolls whatever cell receives focus into view --
    # yanking the viewport back to the first bubble. Focus on the view
    # itself is not a child of the scroll widget, so nothing scrolls;
    # window-level wiring (closed -> composer focus) may then refine it.
    fw = QtWidgets.QApplication.focusWidget()
    if fw is not None and (fw is self.bar or self.bar.isAncestorOf(fw)):
      self._view.setFocus()
    self.bar.hide()
    self.release_hold()
    self.closed.emit()

  def release_hold(self):
    """Release the viewport hold so add_message follows again (send /
    close / query-goes-empty)."""
    self._view.set_autofollow_held(False)

  # ---- scanning ------------------------------------------------------------

  def _scoped_cells(self):
    allowed = self.bar.included_kinds()
    return [w for kind, w in self._view.searchable_cells()
            if kind in allowed]

  def _rescan(self, reset_current=False):
    prev = None
    if not reset_current and 0 <= self._current < len(self._matches):
      prev = self._matches[self._current]
    self._matches = []
    self._capped = False
    query = self.bar.query()
    if query:
      # re.escape + IGNORECASE: literal, case-insensitive matching with
      # offsets that refer to the original string (str.lower()/casefold()
      # can change string LENGTH; regex matching cannot).
      pattern = re.compile(re.escape(query), re.IGNORECASE)
      for w in self._scoped_cells():
        try:
          text = w.toPlainText()
        except RuntimeError:            # C++ widget already deleted
          continue
        for m in pattern.finditer(text):
          self._matches.append((w, m.start(), m.end() - m.start()))
          if len(self._matches) >= MAX_MATCHES:
            self._capped = True
            break
        if self._capped:
          break
    self._current = self._relocate(prev)
    self._update_count()

  def _relocate(self, prev):
    """Index of the previous current match in the fresh match list.

    Exact ``(widget, start)`` identity first; else the nearest-offset
    match in the same widget (the expected path when the match sat in
    the streaming markdown cell, whose plain-text offsets re-flow on
    every delta); else the first match. -1 when there are no matches.
    """
    if not self._matches:
      return -1
    if prev is None:
      return 0
    prev_w, prev_start, _length = prev
    best = None                          # (distance, index)
    for i, (w, start, _n) in enumerate(self._matches):
      if w is not prev_w:
        continue
      if start == prev_start:
        return i
      d = abs(start - prev_start)
      if best is None or d < best[0]:
        best = (d, i)
    return best[1] if best is not None else 0

  # ---- navigation ----------------------------------------------------------

  def find_next(self):
    """Step to the next match; open the bar first when it is closed.

    The window-level Find Next shortcut lands here. From a closed bar
    the retained query's CURRENT match is revealed without stepping
    (reopening must not skip match 1); with it open this is exactly
    the bar's Enter.
    """
    self._find_step(+1)

  def find_previous(self):
    """Step to the previous match; open the bar first when closed.

    Mirror of :meth:`find_next` (the bar's Shift+Enter).
    """
    self._find_step(-1)

  def _find_step(self, delta):
    if not self.is_open():
      self.open()
      if self._matches:
        self._reveal()
      return
    self._go(delta)

  def _go(self, delta):
    """Step to the next/previous match (wrapping) and reveal it.

    The lazy re-scan first picks up anything streamed since the last
    interaction. When a current match existed before the re-scan, step
    ``delta`` from its relocated position; otherwise land ON the first
    match rather than stepping past it.
    """
    had_current = 0 <= self._current < len(self._matches)
    self._rescan(reset_current=False)
    if not self._matches:
      self._apply()
      return
    if had_current:
      self._current = (self._current + delta) % len(self._matches)
      self._update_count()
    self._apply()
    self._reveal()

  def _reveal(self):
    # Engage the hold NOW, in the navigation's own tick -- before the
    # reveal below expands content (rangeChanged would snap to the
    # bottom while follow is still engaged) and before any streamed
    # message can slip in ahead of the deferred scroll and schedule its
    # own snap-to-bottom. ensure_child_rect_visible re-asserts both
    # flags later, but by then the race is already lost without this.
    self._view.set_autofollow_held(True)
    w, start, _length = self._matches[self._current]
    # A match can sit inside collapsed content (e.g. a tool
    # disclosure's body). Ask every ancestor that knows how to reveal
    # hidden searchable content to do so; the duck-typed protocol keeps
    # this controller ignorant of the concrete cell classes.
    p = w.parent()
    while p is not None:
      reveal = getattr(p, "ensure_revealed", None)
      if reveal is not None:
        reveal()
      p = p.parent()
    # Geometry from a just-revealed body (and any auto-height growth)
    # commits on the next event-loop pass; scroll after it settles.
    QtCore.QTimer.singleShot(0, lambda: self._scroll_to(w, start))

  def _scroll_to(self, w, start):
    # A singleShot(0) scheduled by _reveal can land after close() while
    # this object is still alive; a stale deferred scroll must not
    # re-engage the hold (via ensure_child_rect_visible) for a search
    # the user already closed. The WHOLE body sits under the
    # RuntimeError guard: during hard window teardown the bar / cell /
    # view C++ objects may already be gone when the deferred call
    # fires, and even is_open() would raise.
    try:
      if not self.is_open():
        return
      cursor = QtGui.QTextCursor(w.document())
      cursor.setPosition(start)
      rect = w.cursorRect(cursor)
      self._view.ensure_child_rect_visible(w, rect)
    except RuntimeError:                 # widget deleted mid-defer
      return

  # ---- highlighting --------------------------------------------------------

  def _apply(self):
    per_widget = {}                      # id(w) -> (w, [ExtraSelection])
    for i, (w, start, length) in enumerate(self._matches):
      try:
        doc = w.document()
      except RuntimeError:
        continue
      if id(w) not in per_widget:
        per_widget[id(w)] = (w, [])
      sel = QtWidgets.QTextEdit.ExtraSelection()
      cursor = QtGui.QTextCursor(doc)
      cursor.setPosition(start)
      cursor.setPosition(start + length, QtGui.QTextCursor.KeepAnchor)
      sel.cursor = cursor
      bg = _CURRENT_BG if i == self._current else _MATCH_BG
      sel.format.setBackground(QtGui.QColor(bg))
      sel.format.setForeground(QtGui.QColor(_MATCH_FG))
      per_widget[id(w)][1].append(sel)
    stale = [w for w in self._highlighted if id(w) not in per_widget]
    self._highlighted = []
    for w, sels in per_widget.values():
      try:
        w.setExtraSelections(sels)
        self._highlighted.append(w)
      except RuntimeError:
        pass
    for w in stale:
      try:
        w.setExtraSelections([])
      except RuntimeError:
        pass

  def _clear_highlights(self):
    for w in self._highlighted:
      try:
        w.setExtraSelections([])
      except RuntimeError:
        pass
    self._highlighted = []

  def _update_count(self):
    if not self.bar.query():
      self.bar.clear_count()
      return
    current = self._current + 1 if self._current >= 0 else 0
    self.bar.set_count(current, len(self._matches), capped=self._capped)

  # ---- signal handlers -----------------------------------------------------

  def _on_query_changed(self, _text):
    # An emptied query ends the navigation the hold was protecting; let
    # the view follow the bottom again (spec: the hold releases on
    # close, send, and query-goes-empty).
    if not self.bar.query():
      self.release_hold()
    self._rescan(reset_current=True)
    self._apply()

  def _on_scope_changed(self):
    self._rescan(reset_current=True)
    self._apply()
