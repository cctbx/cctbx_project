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


def _req(request_id="r1", tool_name="phenix_start_job",
         source="mcp:phenix", risk="write", batch_id=None):
  return ToolApprovalRequest(
    request_id=request_id, tool_name=tool_name, tool_source=source,
    input={"program": "phenix.refine"}, risk=risk, summary=None,
    batch_id=batch_id)


def exercise_single_card_approve_emits_response():
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
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
  card = ToolApprovalCard()
  card.set_requests([_req()])
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_stop()
  assert decisions[0][0].decision == "deny_and_stop"


def exercise_batched_card_approve_all_emits_n_responses_in_order():
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
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
  card = ToolApprovalCard()
  card.set_requests([_req()])
  card.set_remember_tool(True)
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_approve_all()
  assert decisions[0][0].remember == "tool"


def exercise():
  exercise_single_card_approve_emits_response()
  exercise_single_card_deny_and_stop()
  exercise_batched_card_approve_all_emits_n_responses_in_order()
  exercise_batched_deny_all()
  exercise_remember_tool_checkbox_sets_remember_field()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
