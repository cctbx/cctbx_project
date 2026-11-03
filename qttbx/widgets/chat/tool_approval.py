"""Inline card requesting approval for one or more tool calls.

For a single ``ToolApprovalRequest`` the card shows the tool, three
buttons -- Approve / Deny / Stop -- and an "Always allow this tool"
checkbox. For a batch (multiple requests sharing ``batch_id``) the card
collapses to Approve all / Deny all / Stop. Each response carries its
``request_id``; the ``ApprovalCoordinator`` pairs it with the matching
pending request, so responses need not arrive in dispatch order.
"""

import json

from qttbx.qt import QtCore, QtWidgets

from qttbx.widgets.chat.agent.tools import ToolApprovalResponse
from qttbx.widgets.chat.eliding import ElidingCheckBox, ElidingLabel

# Per-risk text styling for the title (single-card) / per-line
# (batched). 'write' (the most common risk) uses no color override so
# it inherits the active theme's text color -- a hardcoded amber like
# '#c80' renders as low-contrast orange on dark themes. Read uses a
# palette-aware dim grey; destructive uses a red that stays readable
# on light + dark.
_RISK_STYLES = {
  "read":        "color: palette(mid);",
  "write":       "",
  "destructive": "color: #c0392b;",
}

# Characters of the "<source> requests:" lead that survive any squeeze:
# the whole of the longest kind, "builtin", or "mcp:" and three of the
# server. The source is a trust signal -- the Phenix agent labels an
# MCP tool mcp:<server> so it does not read as a built-in -- and a long
# tool name beside it took it to a lone ellipsis.
_SOURCE_KEEP = len("builtin")


class ToolApprovalCard(QtWidgets.QFrame):
  """Single or batched approval card.

  The card emits one ``decided(list)`` per click; the list always
  contains one ``ToolApprovalResponse`` per ``ToolApprovalRequest``, in
  dispatch order.
  """

  decided = QtCore.Signal(list)                  # list[ToolApprovalResponse]

  def __init__(self, parent=None):
    super().__init__(parent)
    self.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self._requests = []
    self._decided = False
    self._remember_tool = False
    self._layout = QtWidgets.QVBoxLayout(self)
    self._layout.setContentsMargins(8, 6, 8, 6)
    self._content_widget = None
    self._buttons_widget = None

  # ---- public API ----------------------------------------------------------

  def set_requests(self, requests):
    self._requests = list(requests)
    self._rebuild()

  def set_remember_tool(self, value):
    self._remember_tool = bool(value)

  def is_decided(self):
    """Report whether a decision has already been emitted.

    Once decided, the card is hidden and its buttons disabled.
    ``ConversationView`` consults this so it never appends a newly
    arrived same-batch request to an already-answered card.

    Returns
    -------
    bool
        ``True`` once a decision has been emitted.
    """
    return self._decided

  def click_approve_all(self):
    self._emit("approve")

  def click_deny_all(self):
    self._emit("deny")

  def click_stop(self):
    self._emit("deny_and_stop")

  def finalize(self):
    """Disable this card WITHOUT emitting a decision.

    Called when the card's turn ends while still undecided (the user stopped
    the turn) so a later click can't emit a stale response into a subsequent
    turn -- the approval-misroute guard. Marks the card decided (so ``_emit``
    won't fire), disables its buttons, and hides it, mirroring the
    post-decision state. A no-op once a real decision has been emitted.
    """
    if self._decided:
      return
    self._decided = True
    if self._buttons_widget is not None:
      self._buttons_widget.setEnabled(False)
    self.hide()

  # ---- UI build ------------------------------------------------------------

  def _rebuild(self):
    if self._content_widget is not None:
      self._content_widget.setParent(None)
    if self._buttons_widget is not None:
      self._buttons_widget.setParent(None)
    # The flag mirrors the checkbox, and the rebuilt one starts unticked
    # -- or is not there at all once a same-batch request has turned the
    # card batched, which must not approve the whole batch for good.
    self._remember_tool = False
    self._content_widget = self._build_content()
    self._buttons_widget = self._build_buttons()
    self._layout.addWidget(self._content_widget)
    self._layout.addWidget(self._buttons_widget)

  def _build_content(self):
    box = QtWidgets.QWidget(self)
    layout = QtWidgets.QVBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
    # No free text in a plain QLabel: a QLabel reports its full text width
    # as its minimum (word wrap does not take it below the widest
    # unbreakable token), so a long MCP tool name or argument floored the
    # ConversationView's container and put a horizontal scrollbar under
    # every narrower window while the tool awaited approval. Names elide
    # in the middle, the whole name on the tooltip; arguments go in a
    # read-only text view that wraps anywhere and grows with its content,
    # so what _short_json leaves of them is all there to read.
    if len(self._requests) == 1:
      r = self._requests[0]
      # The tool's own name beside the source, so what an elided title
      # keeps is the verb: elided as one string it kept the server twice
      # and lost "delete". The registered name stays on the tooltip.
      row, labels = self._name_row(
        box, r, _display_name(r), lead="%s requests:" % r.tool_source)
      for label in labels:
        f = label.font(); f.setBold(True); label.setFont(f)
      layout.addLayout(row)
      layout.addWidget(_args_view(box, _short_json(r.input, 200)))
    else:
      head = QtWidgets.QLabel(
        "Batched tool requests (%d)" % len(self._requests), box)
      f = head.font(); f.setBold(True); head.setFont(f)
      layout.addWidget(head)
      # The registered name on every line: with their servers stripped,
      # mcp__coot__delete_all and mcp__phenix__delete_all read as two
      # identical "delete_all" lines above "Approve all". Indented
      # rather than bulleted: a bullet in the name's text put it on the
      # tooltip, and one in a label of its own elided to an ellipsis.
      # The preview is indented under its name, one block per request.
      # It keeps the theme's text color. The one-label line it replaced
      # painted the arguments in the name's risk color too, but a color
      # stylesheet on a text view is the trap set_transparent_background
      # describes, and the arguments carry no risk the name does not
      # already show.
      for r in self._requests:
        block = QtWidgets.QVBoxLayout()
        block.setContentsMargins(12, 0, 0, 0)
        row, _ = self._name_row(box, r, r.tool_name)
        block.addLayout(row)
        block.addWidget(_args_view(box, _short_json(r.input, 120)))
        layout.addLayout(block)
    return box

  def _name_row(self, parent, r, name, lead=None):
    """A row reading ``lead``, if any, then ``name``, the tool's name.

    Eliding labels both, so each gives way on its own and neither
    floors the card: a plain label for the source floored it at a long
    user-configured server name, with the name squeezed to a lone
    ellipsis. Tool name and source are externally supplied: rendered
    literally, on the tooltip too. The registered name goes on the
    tooltip when the label shows a shortened one.
    """
    row = QtWidgets.QHBoxLayout()
    row.setContentsMargins(0, 0, 0, 0)
    # The lead elides at its end with its first characters kept, so the
    # source's kind outlasts the squeeze: elided in the middle it lost
    # the server first ("mc…s:"), and with nothing kept a long tool
    # name took it to a lone ellipsis. The rest is on its tooltip. The
    # name elides in the middle: its verb is at one end and its object
    # at the other.
    specs = [(QtCore.Qt.ElideMiddle, 0, name,
              None if name == r.tool_name else r.tool_name)]
    if lead is not None:
      specs.insert(0, (QtCore.Qt.ElideRight, _SOURCE_KEEP, lead, None))
    labels = []
    for mode, keep, text, tooltip in specs:
      label = ElidingLabel(parent, mode=mode, keep=keep)
      label.setStyleSheet(_RISK_STYLES.get(r.risk, ""))
      label.set_full_text(text, tooltip=tooltip)
      labels.append(label)
    # The name is what the user decides on, so it is the last to give
    # way: a layout squeezes a stretched item first (its hint counts as
    # its minimum) and takes the rest equally from the others, so the
    # stretch goes on the lead. A Maximum policy keeps the lead from
    # growing past its text, and the trailing spacer carries no stretch
    # of its own, or the two would split the room and the lead elide
    # with half the row empty; expanding, it still takes what the lead
    # cannot use.
    if lead is not None:
      policy = labels[0].sizePolicy()
      policy.setHorizontalPolicy(QtWidgets.QSizePolicy.Maximum)
      labels[0].setSizePolicy(policy)
      row.addWidget(labels[0], 1)
    row.addWidget(labels[-1])
    row.addStretch()
    return row, labels

  def _build_buttons(self):
    box = QtWidgets.QWidget(self)
    rows = QtWidgets.QVBoxLayout(box)
    rows.setContentsMargins(0, 0, 0, 0)
    layout = QtWidgets.QHBoxLayout()
    layout.setContentsMargins(0, 0, 0, 0)
    rows.addLayout(layout)
    if len(self._requests) <= 1:
      # Offer "Always allow this tool" only when the request permits it; a
      # request with allow_remember=False (e.g. a destructive recovery tool)
      # must not grant standing per-tool auto-approval.
      allow_remember = (getattr(self._requests[0], "allow_remember", True)
                        if self._requests else True)
      approve = QtWidgets.QPushButton("Approve", box)
      deny = QtWidgets.QPushButton("Deny", box)
      stop = QtWidgets.QPushButton("Stop", box)
      approve.clicked.connect(self.click_approve_all)
      deny.clicked.connect(self.click_deny_all)
      stop.clicked.connect(self.click_stop)
      layout.addWidget(approve)
      layout.addWidget(deny)
      layout.addStretch(1)
      layout.addWidget(stop)
      if allow_remember:
        # Its own row under the buttons: in their row it was the only
        # shrinkable item, so it absorbed the whole shortfall and read
        # "Always all..." at the ChatWindow's narrowest column. Eliding
        # still, as the labels above: a plain QCheckBox's minimum is its
        # full text width, and under a font wider than the buttons' that
        # would floor the card.
        always = ElidingCheckBox(box)
        always.set_full_text("Always allow this tool")
        always.toggled.connect(self.set_remember_tool)
        rows.addWidget(always)
    else:
      approve_all = QtWidgets.QPushButton("Approve all", box)
      deny_all = QtWidgets.QPushButton("Deny all", box)
      stop = QtWidgets.QPushButton("Stop", box)
      approve_all.clicked.connect(self.click_approve_all)
      deny_all.clicked.connect(self.click_deny_all)
      stop.clicked.connect(self.click_stop)
      layout.addWidget(approve_all)
      layout.addWidget(deny_all)
      layout.addStretch(1)
      layout.addWidget(stop)
    return box

  # ---- emit ----------------------------------------------------------------

  def _emit(self, decision):
    if self._decided:
      return               # already decided or finalized -- never emit twice
    remember = "none"
    if decision == "approve" and self._remember_tool:
      remember = "tool"
    responses = [
      ToolApprovalResponse(request_id=r.request_id,
                           decision=decision, remember=remember)
      for r in self._requests
    ]
    # Card has served its purpose -- disable so a stray re-click can't
    # emit again, then hide so the conversation view doesn't keep a
    # stale prompt around. Qt's QVBoxLayout skips hidden widgets so
    # the visible layout collapses naturally.
    self._decided = True
    if self._buttons_widget is not None:
      self._buttons_widget.setEnabled(False)
    self.hide()
    self.decided.emit(responses)


def _display_name(r):
  """``r``'s tool name without an ``mcp__<server>__`` prefix that only
  repeats its source (``mcp:<server>``), so what an elided title keeps is
  the tool's own name -- the verb -- and not the server twice. Not the
  ``<server>__`` of a tool the registry renamed for a collision: that is
  what tells the two apart, and what "Always allow" is keyed by."""
  if r.tool_source.startswith("mcp:"):
    prefix = "mcp__%s__" % r.tool_source[len("mcp:"):]
    if r.tool_name.startswith(prefix) and len(r.tool_name) > len(prefix):
      return r.tool_name[len(prefix):]
  return r.tool_name


def _args_view(parent, text):
  """Read-only argument preview that wraps anywhere and grows with its
  content: the tool disclosure's argument view, in miniature. Unlike a
  label it neither floors the card at its widest token nor clips one.
  Plain text, so untrusted tool arguments can't inject rich text."""
  from qttbx.widgets.chat.auto_height import plain_text_view
  view = plain_text_view(parent, fixed_pitch=True)
  # Not a Tab stop: it sat, invisibly, ahead of Approve, and while it
  # held the focus Space and PageDown went to it and not the view.
  view.setFocusPolicy(QtCore.Qt.NoFocus)
  view.setPlainText(text)
  return view


def _short_json(d, limit=80):
  try:
    s = json.dumps(d, default=str)
  except Exception:
    s = repr(d)
  if len(s) > limit:
    s = s[:limit - 1] + "..."
  return s
