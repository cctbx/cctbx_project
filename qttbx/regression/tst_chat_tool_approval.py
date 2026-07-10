"""Tool approval card tests. Covers single-request card, batched card,
and the deny/approve/deny-and-stop response wiring."""

import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times

try:
  from qttbx.qt import QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)

from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
from qttbx.widgets.font_init import init_default_app_font


def _req(request_id="r1", tool_name="phenix_start_job",
         source="mcp:phenix", risk="write", batch_id=None):
  return ToolApprovalRequest(
    request_id=request_id, tool_name=tool_name, tool_source=source,
    input={"program": "phenix.refine"}, risk=risk, batch_id=batch_id)


def exercise_single_card_approve_emits_response():
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req()])
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_approve_all()
  assert len(decisions) == 1
  assert decisions[0][0].decision == "approve"
  assert decisions[0][0].request_id == "r1"


def exercise_single_card_deny_and_stop():
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req()])
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_stop()
  assert decisions[0][0].decision == "deny_and_stop"


def exercise_batched_card_approve_all_emits_n_responses_in_order():
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([
    _req(request_id="a", batch_id="B"),
    _req(request_id="b", batch_id="B"),
    _req(request_id="c", batch_id="B"),
  ])
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_approve_all()
  assert len(decisions) == 1
  ids = [r.request_id for r in decisions[0]]
  assert ids == ["a", "b", "c"], ids
  assert all(r.decision == "approve" for r in decisions[0])


def exercise_batched_deny_all():
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req(request_id="a", batch_id="B"),
                     _req(request_id="b", batch_id="B")])
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_deny_all()
  assert [r.decision for r in decisions[0]] == ["deny", "deny"]


def exercise_remember_tool_checkbox_sets_remember_field():
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req()])
  card.set_remember_tool(True)
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_approve_all()
  assert decisions[0][0].remember == "tool"


def exercise_card_hides_and_disables_buttons_after_click():
  """After a decision is emitted, the card must (a) hide so the
  conversation view doesn't carry a stale prompt forward, and (b)
  disable its button row so a double-click before Qt removes the
  widget can't fire the decision a second time. Applies to every
  button -- approve, deny, and stop."""
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  for click_method in ("click_approve_all", "click_deny_all", "click_stop"):
    card = ToolApprovalCard()
    card.set_requests([_req()])
    card.show()
    decisions = []
    card.decided.connect(lambda r: decisions.append(r))
    assert not card.isHidden(), click_method
    getattr(card, click_method)()
    assert len(decisions) == 1, (click_method, decisions)
    assert card.isHidden(), click_method
    assert not card._buttons_widget.isEnabled(), click_method


def exercise_is_decided_reflects_decision_state():
  """is_decided() is False until the card emits a decision, then True.
  ConversationView consults it to tell an already-answered batch card
  from a live one when a serially-dispatched same-batch request arrives,
  so it never appends a new request to a hidden, decided card."""
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req()])
  assert not card.is_decided()
  card.click_approve_all()
  assert card.is_decided()


def exercise_untrusted_input_rendered_as_plain_text():
  """The approval card shows untrusted tool input (and tool name/source)
  -- the very surface the user reads to decide. Its labels must use
  PlainText so a tool argument containing HTML can't inject rich text
  (e.g. a clickable file:// link) into the decision prompt."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  evil = "<a href='file:///etc/passwd'>ok</a>"
  # Single-request card: the input preview label.
  card = ToolApprovalCard()
  card.set_requests([ToolApprovalRequest(
    request_id="r1", tool_name="phenix_start_job", tool_source="mcp:phenix",
    input={"path": evil}, risk="write", batch_id=None)])
  hits = [l for l in card.findChildren(QtWidgets.QLabel)
          if "file:///" in l.text()]
  assert hits, [l.text() for l in card.findChildren(QtWidgets.QLabel)]
  for l in hits:
    assert l.textFormat() == QtCore.Qt.PlainText, l.text()
  # Batched card: the per-request line labels.
  card2 = ToolApprovalCard()
  card2.set_requests([
    ToolApprovalRequest(request_id="a", tool_name="t1", tool_source="s",
                        input={"p": evil}, risk="write", batch_id="B"),
    ToolApprovalRequest(request_id="b", tool_name="t2", tool_source="s",
                        input={"p": "ok"}, risk="write", batch_id="B")])
  hits2 = [l for l in card2.findChildren(QtWidgets.QLabel)
           if "file:///" in l.text()]
  assert hits2, [l.text() for l in card2.findChildren(QtWidgets.QLabel)]
  for l in hits2:
    assert l.textFormat() == QtCore.Qt.PlainText, l.text()


def exercise():
  exercise_single_card_approve_emits_response()
  exercise_single_card_deny_and_stop()
  exercise_batched_card_approve_all_emits_n_responses_in_order()
  exercise_batched_deny_all()
  exercise_remember_tool_checkbox_sets_remember_field()
  exercise_card_hides_and_disables_buttons_after_click()
  exercise_is_decided_reflects_decision_state()
  exercise_untrusted_input_rendered_as_plain_text()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
