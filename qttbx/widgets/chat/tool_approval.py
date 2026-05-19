"""``ToolApprovalCard`` — inline card requesting approval for one or
more tool calls.

For a single ``ToolApprovalRequest`` the card shows the tool and four
buttons: Approve / Deny / Always allow / Stop. For a batch (multiple
requests sharing ``batch_id``) the card collapses to Approve all /
Deny all / Stop. The worker receives batched responses in dispatch
order so the session's ``_dispatch_and_build_results`` loop can pair
each response with its request.
"""

import json

from qttbx.qt import QtCore, QtWidgets

from qttbx.widgets.chat.agent.tools import ToolApprovalResponse

_RISK_COLORS = {
  "read": "#888",
  "write": "#c80",
  "destructive": "#c33",
}


class ToolApprovalCard(QtWidgets.QFrame):
  """Single or batched approval card. The card emits one decided(list) per
  click; the list always contains one ToolApprovalResponse per
  ToolApprovalRequest, in dispatch order."""

  decided = QtCore.Signal(list)                  # list[ToolApprovalResponse]

  def __init__(self, parent=None):
    super().__init__(parent)
    self.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self._requests = []
    self._remember_tool = False
    self._remember_server = False
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

  def set_remember_server(self, value):
    self._remember_server = bool(value)

  def click_approve_all(self):
    self._emit("approve")

  def click_deny_all(self):
    self._emit("deny")

  def click_stop(self):
    self._emit("deny_and_stop")

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
      f = head.font(); f.setBold(True); head.setFont(f)
      head.setStyleSheet(
        "color: %s;" % _RISK_COLORS.get(r.risk, "#888"))
      layout.addWidget(head)
      preview = QtWidgets.QLabel(_short_json(r.input, 200), box)
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
        line.setWordWrap(True)
        line.setStyleSheet(
          "color: %s;" % _RISK_COLORS.get(r.risk, "#888"))
        layout.addWidget(line)
    return box

  def _build_buttons(self):
    box = QtWidgets.QWidget(self)
    layout = QtWidgets.QHBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
    if len(self._requests) <= 1:
      approve = QtWidgets.QPushButton("Approve", box)
      always = QtWidgets.QCheckBox("Always allow this tool", box)
      always.toggled.connect(self.set_remember_tool)
      deny = QtWidgets.QPushButton("Deny", box)
      stop = QtWidgets.QPushButton("Stop", box)
      approve.clicked.connect(self.click_approve_all)
      deny.clicked.connect(self.click_deny_all)
      stop.clicked.connect(self.click_stop)
      layout.addWidget(approve)
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
    remember = "none"
    if decision == "approve":
      if self._remember_server:
        remember = "server"
      elif self._remember_tool:
        remember = "tool"
    responses = [
      ToolApprovalResponse(request_id=r.request_id,
                           decision=decision, remember=remember)
      for r in self._requests
    ]
    self.decided.emit(responses)


def _short_json(d, limit=80):
  try:
    s = json.dumps(d, default=str)
  except Exception:
    s = repr(d)
  if len(s) > limit:
    s = s[:limit - 1] + "..."
  return s
