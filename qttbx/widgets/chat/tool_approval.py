"""Inline card requesting approval for one or more tool calls.

For a single ``ToolApprovalRequest`` the card shows the tool and four
buttons: Approve / Deny / Always allow / Stop. For a batch (multiple
requests sharing ``batch_id``) the card collapses to Approve all /
Deny all / Stop. Each response carries its ``request_id``; the
``ApprovalCoordinator`` pairs it with the matching pending request, so
responses need not arrive in dispatch order.
"""

import json

from qttbx.qt import QtCore, QtWidgets

from qttbx.widgets.chat.agent.tools import ToolApprovalResponse

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
    self._content_widget = self._build_content()
    self._buttons_widget = self._build_buttons()
    self._layout.addWidget(self._content_widget)
    self._layout.addWidget(self._buttons_widget)

  def _build_content(self):
    box = QtWidgets.QWidget(self)
    layout = QtWidgets.QVBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
    if len(self._requests) == 1:
      r = self._requests[0]
      title = "%s requests: %s" % (r.tool_source, r.tool_name)
      head = QtWidgets.QLabel(title, box)
      # tool name/source are externally supplied: render literally.
      head.setTextFormat(QtCore.Qt.PlainText)
      f = head.font(); f.setBold(True); head.setFont(f)
      head.setStyleSheet(_RISK_STYLES.get(r.risk, ""))
      layout.addWidget(head)
      # The input preview is untrusted tool-argument data shown on the
      # security-decision surface: PlainText so it can't inject rich text.
      preview = QtWidgets.QLabel(_short_json(r.input, 200), box)
      preview.setTextFormat(QtCore.Qt.PlainText)
      preview.setWordWrap(True)
      layout.addWidget(preview)
    else:
      head = QtWidgets.QLabel(
        "Batched tool requests (%d)" % len(self._requests), box)
      f = head.font(); f.setBold(True); head.setFont(f)
      layout.addWidget(head)
      for r in self._requests:
        line = QtWidgets.QLabel(
          "  - %s - %s" % (r.tool_name, _short_json(r.input, 120)), box)
        line.setTextFormat(QtCore.Qt.PlainText)
        line.setWordWrap(True)
        line.setStyleSheet(_RISK_STYLES.get(r.risk, ""))
        layout.addWidget(line)
    return box

  def _build_buttons(self):
    box = QtWidgets.QWidget(self)
    layout = QtWidgets.QHBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
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
      if allow_remember:
        always = QtWidgets.QCheckBox("Always allow this tool", box)
        always.toggled.connect(self.set_remember_tool)
        layout.addWidget(always)
      layout.addWidget(deny)
      layout.addStretch(1)
      layout.addWidget(stop)
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


def _short_json(d, limit=80):
  try:
    s = json.dumps(d, default=str)
  except Exception:
    s = repr(d)
  if len(s) > limit:
    s = s[:limit - 1] + "..."
  return s
