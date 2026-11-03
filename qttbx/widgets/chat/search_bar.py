"""Floating find-in-conversation bar.

Pure UI: a query line edit, a match-count label, previous/next arrows,
a "Tools" scope checkbox (default OFF), and a close button. The widget owns no search logic; ConversationSearch drives it
through signals and the small query/scope/count API.

Key handling while the bar has focus: Enter -> next, Shift+Enter ->
previous, Escape -> close. Escape is claimed through the
ShortcutOverride event so it outranks the window-level Escape=stop
shortcut ONLY while focus is inside the bar; Escape anywhere else still
stops the turn.
"""

from qttbx.qt import QtCore, QtWidgets

# Kinds searched unconditionally: regular text and thinking cells both
# read as the conversation itself, so neither is toggleable.
ALWAYS_SEARCHED = frozenset(("text", "thinking"))

# Optional search scopes: (kind, checkbox label, tooltip, default_on).
# Each entry gets a checkbox on the bar, keyed by the same kind string
# the cell classes report from their searchable_cells(). Adding a cell
# kind = one row here (or in ALWAYS_SEARCHED) plus the cell class
# reporting it. Tool payloads are bulky, so they stay opt-in.
OPTIONAL_SCOPES = (
  ("tool", "Tools", "Also search tool-call args and results", False),
)


class SearchBar(QtWidgets.QFrame):
  """Find bar floated over the conversation by ConversationSearch."""

  query_changed = QtCore.Signal(str)
  next_requested = QtCore.Signal()
  prev_requested = QtCore.Signal()
  scope_changed = QtCore.Signal()
  close_requested = QtCore.Signal()

  def __init__(self, parent=None):
    super().__init__(parent)
    # Opaque panel: the bar floats over conversation text.
    self.setAutoFillBackground(True)
    self.setFrameStyle(QtWidgets.QFrame.StyledPanel)
    layout = QtWidgets.QHBoxLayout(self)
    layout.setContentsMargins(6, 4, 6, 4)
    layout.setSpacing(4)
    self._edit = QtWidgets.QLineEdit(self)
    self._edit.setPlaceholderText("Find in conversation")
    self._edit.setClearButtonEnabled(True)
    self._edit.textChanged.connect(self.query_changed)
    layout.addWidget(self._edit, 1)
    self._count = QtWidgets.QLabel("", self)
    self._count.setStyleSheet("color: palette(mid);")
    layout.addWidget(self._count)
    self._prev_btn = QtWidgets.QToolButton(self)
    self._prev_btn.setText("▲")
    self._prev_btn.setToolTip("Previous match (Shift+Enter)")
    self._prev_btn.setAutoRaise(True)
    self._prev_btn.clicked.connect(
      lambda _checked=False: self.prev_requested.emit())
    layout.addWidget(self._prev_btn)
    self._next_btn = QtWidgets.QToolButton(self)
    self._next_btn.setText("▼")
    self._next_btn.setToolTip("Next match (Enter)")
    self._next_btn.setAutoRaise(True)
    self._next_btn.clicked.connect(
      lambda _checked=False: self.next_requested.emit())
    layout.addWidget(self._next_btn)
    self._scope_boxes = {}
    for kind, label, tip, default_on in OPTIONAL_SCOPES:
      box = QtWidgets.QCheckBox(label, self)
      box.setToolTip(tip)
      # Default BEFORE connecting: the initial check must not fire
      # scope_changed into a half-built controller.
      box.setChecked(default_on)
      box.toggled.connect(lambda _on: self.scope_changed.emit())
      layout.addWidget(box)
      # Breathing room: the checkbox label would otherwise sit only
      # the default 4 px from the next control.
      layout.addSpacing(8)
      self._scope_boxes[kind] = box
    self._close_btn = QtWidgets.QToolButton(self)
    self._close_btn.setText("✕")
    self._close_btn.setToolTip("Close search (Esc)")
    self._close_btn.setAutoRaise(True)
    self._close_btn.clicked.connect(
      lambda _checked=False: self.close_requested.emit())
    layout.addWidget(self._close_btn)
    # Escape/Enter must work with focus ANYWHERE inside the bar, not
    # just the line edit: a checkbox takes click focus, and an Escape
    # there would otherwise fall through to the window-level
    # Escape=stop shortcut and kill the running turn. One shared filter
    # on every interactive child; the handler treats them alike (Space
    # and Tab are not intercepted, so checkbox toggling and focus
    # traversal stay native).
    children = [self._edit, self._prev_btn, self._next_btn]
    children.extend(self._scope_boxes.values())
    children.append(self._close_btn)
    for w in children:
      w.installEventFilter(self)

  # ---- controller API ------------------------------------------------------

  def query(self):
    return self._edit.text()

  def included_kinds(self):
    """Kinds currently in scope.

    Returns
    -------
    set of str
        The always-searched kinds (text, thinking) plus the kind of
        each checked optional-scope checkbox.
    """
    kinds = set(ALWAYS_SEARCHED)
    for kind, box in self._scope_boxes.items():
      if box.isChecked():
        kinds.add(kind)
    return kinds

  def set_count(self, current, total, capped=False):
    """Show ``current/total`` (1-based; ``0/0`` = no matches).

    Parameters
    ----------
    current : int
        1-based index of the current match; 0 when there is none.
    total : int
        Number of collected matches.
    capped : bool, optional
        True when scanning stopped at the match cap; renders ``total+``.
    """
    total_text = "%d+" % total if capped else "%d" % total
    self._count.setText("%d/%s" % (current, total_text))

  def clear_count(self):
    """Blank the count label (empty query)."""
    self._count.setText("")

  def focus_query(self):
    """Focus the line edit with the retained query pre-selected."""
    self._edit.setFocus()
    self._edit.selectAll()

  # ---- key handling --------------------------------------------------------

  def eventFilter(self, obj, event):
    if self.isAncestorOf(obj):
      t = event.type()
      if t == QtCore.QEvent.ShortcutOverride:
        # Claim Escape while the bar has focus so the window-level
        # Escape=stop shortcut doesn't fire; the KeyPress then arrives
        # below and closes the bar.
        if event.key() == QtCore.Qt.Key_Escape:
          event.accept()
          return True
      elif t == QtCore.QEvent.KeyPress:
        key = event.key()
        if key == QtCore.Qt.Key_Escape:
          self.close_requested.emit()
          return True
        if key in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
          if event.modifiers() & QtCore.Qt.ShiftModifier:
            self.prev_requested.emit()
          else:
            self.next_requested.emit()
          return True
    return super().eventFilter(obj, event)
