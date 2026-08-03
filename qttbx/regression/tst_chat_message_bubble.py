"""MessageBubble widget test. Covers: text rendering, tool_use cell,
thinking block, streaming append, and image placeholder rendering."""

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

from qttbx.widgets.chat.agent.conversation import ContentBlock, Message, now
from qttbx.widgets.font_init import init_default_app_font


def _qapp():
  return QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)


# The shape of error a force-killed Coot bridge produces: long, single-line,
# and routed to the failing tool's disclosure cell. Kept verbatim in
# tst_chat_conversation_widgets.py, which exercises the same failure one
# layer out (at the ConversationView) -- change both together.
LONG_TOOL_ERROR = (
  "MCP error -32000: Connection closed. Failed to connect to the Coot RPC "
  "bridge at 127.0.0.1:44100: [Errno 61] Connection refused. The Coot "
  "process may have exited or been killed; restart Coot and retry.")


def exercise_renders_user_text():
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  m = Message(role="user", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "hello world"})])
  b = MessageBubble(m)
  assert "hello world" in b.combined_text()


def exercise_renders_tool_use_cell():
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "calling tool"}),
    ContentBlock(type="tool_use", data={
      "id": "t1", "name": "phenix_start_job",
      "input": {"program": "phenix.refine"}}),
  ])
  b = MessageBubble(m)
  assert "phenix_start_job" in b.combined_text()


def exercise_renders_thinking_block():
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="thinking", data={"text": "I should think hard"}),
    ContentBlock(type="text", data={"text": "answer"}),
  ])
  b = MessageBubble(m)
  text = b.combined_text()
  assert "answer" in text
  assert "think hard" in text


def exercise_image_block_renders_as_placeholder():
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="image", data={
      "attachment_sha256": "abcdef0123456789",
      "mime": "image/png", "caption": None}),
  ])
  b = MessageBubble(m)
  text = b.combined_text()
  assert "image/png" in text
  assert "abcdef01" in text


def exercise_streaming_append():
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": ""})])
  b = MessageBubble(m)
  b.append_text_delta("hel")
  b.append_text_delta("lo")
  assert "hello" in b.combined_text()


def exercise_thinking_delta_after_text_starts_new_block():
  """Streaming thinking AFTER text should append a new thinking block, not
  retroactively extend an earlier thinking block. Mirrors the strict
  last-block semantics of append_text_delta."""
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="thinking", data={"text": "T1"}),
    ContentBlock(type="text", data={"text": "middle"}),
  ])
  b = MessageBubble(m)
  b.append_thinking_delta("T2")
  # Original thinking block must be unchanged.
  assert m.content[0].data["text"] == "T1", m.content[0].data["text"]
  # A new thinking block must have been appended.
  assert m.content[-1].type == "thinking"
  assert m.content[-1].data["text"] == "T2"


def exercise_image_cell_renders_real_image():
  import shutil, tempfile
  from pathlib import Path
  from qttbx.widgets.chat.agent.storage import ConversationStorage
  from qttbx.widgets.chat.message_bubble import MessageBubble
  from libtbx.utils import null_out

  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    png = _png_bytes(40, 30)
    att = storage.store_attachment("c1", png, "image/png")
    m = Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="image", data={
        "attachment_sha256": att.sha256,
        "mime": "image/png", "caption": "hi"})])
    bub = MessageBubble(m, storage=storage, conv_id="c1")
    # Walk children to find the _ImageCell. Its inner pixmap label should
    # carry a non-null QPixmap.
    pixmap_labels = [
      w for w in bub.findChildren(QtWidgets.QLabel)
      if w.pixmap() is not None and not w.pixmap().isNull()]
    assert pixmap_labels, "no rendered image label found"
  finally:
    shutil.rmtree(tmp)


def exercise_image_cell_click_emits_signal():
  import shutil, tempfile
  from pathlib import Path
  from qttbx.widgets.chat.agent.storage import ConversationStorage
  from qttbx.widgets.chat.message_bubble import MessageBubble
  from libtbx.utils import null_out

  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    att = storage.store_attachment("c1", _png_bytes(), "image/png")
    m = Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="image", data={
        "attachment_sha256": att.sha256, "mime": "image/png"})])
    bub = MessageBubble(m, storage=storage, conv_id="c1")
    fired = []
    bub.image_clicked.connect(lambda c, s: fired.append((c, s)))
    cells = [w for w in bub.findChildren(QtWidgets.QFrame)
             if hasattr(w, "click") and hasattr(w, "sha256")]
    assert cells, "no _ImageCell with .click()"
    cells[0].click()
    assert fired == [("c1", att.sha256)], fired
  finally:
    shutil.rmtree(tmp)


def exercise_bubble_has_no_frame_border():
  """Flattened bubble: QFrame style is NoFrame and no fillable
  background. Spec section 3 -- spacing between turns, not borders."""
  from qttbx.widgets.chat.message_bubble import MessageBubble
  _qapp()
  bubble = MessageBubble(role="user")
  assert bubble.frameStyle() == QtWidgets.QFrame.NoFrame, \
    bubble.frameStyle()


def exercise_role_renders_as_bold_prefix_on_first_text_line():
  """The role word is bold-prefixed onto the first text cell as
  'You: ...' or 'Claude: ...'."""
  from qttbx.widgets.chat.message_bubble import MessageBubble
  _qapp()
  bubble = MessageBubble(role="user")
  bubble.append_text_delta("refine 1yjp")
  text = bubble.first_text_cell_html()
  assert "You" in text and "refine 1yjp" in text, text


def exercise_assistant_label_from_backend_stamp():
  """An assistant message stamped with a backend renders the matching display
  name (openai -> 'GPT'), not a hard-coded 'Claude'."""
  from qttbx.widgets.chat.agent.conversation import Message, now
  from qttbx.widgets.chat.message_bubble import MessageBubble
  _qapp()
  msg = Message(role="assistant", timestamp=now(), content=[], backend="openai")
  bubble = MessageBubble(msg)
  bubble.append_text_delta("I will run phenix.refine.")
  text = bubble.first_text_cell_html()
  assert "GPT" in text and "phenix.refine" in text, text
  assert "Claude" not in text, text


def exercise_assistant_label_falls_back_to_passed_name_then_assistant():
  """With no backend stamp the bubble uses the assistant_label passed in (the
  current backend's name); with neither it is the generic 'Assistant'."""
  from qttbx.widgets.chat.message_bubble import MessageBubble
  _qapp()
  b1 = MessageBubble(role="assistant", assistant_label="Gemini")
  b1.append_text_delta("hi")
  assert "Gemini" in b1.first_text_cell_html()
  b2 = MessageBubble(role="assistant")
  b2.append_text_delta("hi")
  t2 = b2.first_text_cell_html()
  assert "Assistant" in t2 and "Claude" not in t2, t2


def exercise_tool_use_cell_uses_disclosure_widget():
  """The tool-use cell must be a ToolCallDisclosure -- collapsed by
  default, expands on click, status flips from 'running' to a
  terminal state when the bubble is told the tool finished."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.widgets.chat.message_bubble import MessageBubble
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  bubble = MessageBubble(role="assistant")
  cell = bubble.add_tool_use_cell(
    tool_id="t1", name="phenix_start_job",
    args={"program": "phenix.refine", "args": ["model.pdb"]})
  assert isinstance(cell, ToolCallDisclosure), type(cell)
  assert cell.body.isHidden(), "tool cell starts collapsed"
  bubble.set_tool_use_finished(
    tool_id="t1", elapsed="12s",
    result="R-free: 0.32 -> 0.27")
  assert "finished, 12s" in cell.header_button.text()
  cell.header_button.click()
  assert "R-free" in cell.result_view.toPlainText()


def exercise_text_after_tool_cell_does_not_merge_into_prior_text_view():
  """After a tool cell is inserted, the next text block must land in a
  fresh MarkdownView (not the one used before the tool). Otherwise the
  visual ordering text->tool->text collapses into text+text (with the
  tool below them, out of order)."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.widgets.chat.message_bubble import MessageBubble
  bubble = MessageBubble(role="assistant")
  bubble.append_text_delta("I'll run the job.")
  first_view = bubble._text_view
  assert first_view is not None
  bubble.add_tool_use_cell(tool_id="t1", name="phenix_start_job", args={})
  # After a non-text widget, the bubble must drop its reference to
  # _text_view so the NEXT add_text() creates a new MarkdownView.
  assert bubble._text_view is None, \
    "tool cell insertion must reset _text_view"
  bubble.append_text_delta("Now waiting for it to finish.")
  second_view = bubble._text_view
  assert second_view is not None
  assert second_view is not first_view, \
    "text after a tool cell must use a fresh MarkdownView"


def exercise_image_cell_renders_inline_thumbnail_with_click_to_lightbox():
  """Image cell: bounded-height thumbnail in the bubble; click opens
  ImageLightbox. Pins the inline-rendering contract so a regression
  that routes images back to a side panel surfaces here."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.image_lightbox import ImageLightbox
  from qttbx.widgets.chat.message_bubble import MessageBubble
  pm = QtGui.QPixmap(800, 600)
  pm.fill(QtGui.QColor(QtCore.Qt.red))
  bubble = MessageBubble(role="assistant")
  cell = bubble.add_image_cell(pixmap=pm, caption="model.pdb electron density")
  # Thumbnail height capped at MAX_HEIGHT (240), preserving aspect.
  assert cell.thumbnail.pixmap().height() <= 240, \
    cell.thumbnail.pixmap().height()
  # Caption rendered below the image (if present).
  assert "electron density" in cell.caption_label.text()
  # Open the lightbox -- cell exposes the dialog instance for test
  # inspection (set in _open_lightbox before exec).
  cell._open_lightbox()
  assert isinstance(cell._last_lightbox, ImageLightbox)
  assert cell._last_lightbox.pixmap is pm


def exercise_image_cell_without_caption_hides_caption_label():
  """No caption => caption_label is hidden, not just empty."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.message_bubble import MessageBubble
  pm = QtGui.QPixmap(100, 100)
  pm.fill(QtGui.QColor(QtCore.Qt.green))
  bubble = MessageBubble(role="assistant")
  cell = bubble.add_image_cell(pixmap=pm, caption=None)
  assert cell.caption_label.isHidden()


def exercise_server_tool_use_then_result_folds_into_same_disclosure():
  """A server_tool_use block creates a ToolCallDisclosure cell; the
  matching server_tool_result folds into it as the result body (same
  pattern the regular tool_use / tool_result chain uses). The
  display name carries a 'server:' prefix so the user can tell
  API-executed tools apart from client-dispatched ones."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.widgets.chat.agent.conversation import ContentBlock
  from qttbx.widgets.chat.message_bubble import MessageBubble
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  bubble = MessageBubble(role="assistant")
  bubble.append_block(ContentBlock(type="server_tool_use", data={
    "id": "srvtu_1", "name": "web_search",
    "input": {"query": "phenix.refine"}}))
  cell = bubble._tool_cells_by_id.get("srvtu_1")
  assert isinstance(cell, ToolCallDisclosure), cell
  assert "server:" in cell.header_button.text()
  assert "web_search" in cell.header_button.text()
  bubble.append_block(ContentBlock(type="server_tool_result", data={
    "tool_use_id": "srvtu_1",
    "content": {"type": "web_search_tool_result", "content": [
      {"type": "web_search_result",
       "title": "phenix.refine — Phenix",
       "url": "https://phenix-online.org/x.html"}]}}))
  cell.header_button.click()
  body = cell.result_view.toPlainText()
  assert "phenix.refine" in body, body
  assert "phenix-online.org" in body, body


def exercise_orphan_tool_result_renders_collapsed_disclosure():
  """A tool_result whose matching tool_use is NOT in this bubble (an
  'orphan' -- the normal shape when a stored conversation is reloaded,
  since the tool_use lives in the assistant message and the tool_result
  in the following user message) must render as a collapsed
  ToolCallDisclosure, not a full-text label. Large Phenix tool output
  (phenix_get_phil can exceed 50K chars) is hidden behind the disclosure
  and shown only when the user clicks to expand."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.widgets.chat.message_bubble import MessageBubble
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  marker = "PHIL_SCOPE_MARKER_" + ("x" * 8000)        # stand-in for bulk output
  m = Message(role="user", timestamp=now(), content=[
    ContentBlock(type="tool_result", data={
      "tool_use_id": "toolu_orphan", "is_error": False,
      "content": [ContentBlock(type="text", data={"text": marker})]})])
  b = MessageBubble(m)
  discs = b.findChildren(ToolCallDisclosure)
  assert len(discs) == 1, "orphan tool_result must render one disclosure"
  disc = discs[0]
  # Collapsed by default: the bulk output is NOT visible until clicked.
  assert disc.body.isHidden(), "orphan result disclosure must start collapsed"
  # The bulk text must NOT be dumped into a wrapped QLabel (the old
  # _ToolResultCell behavior that showed 50K chars in full).
  assert not any(marker in lbl.text()
                 for lbl in b.findChildren(QtWidgets.QLabel)), \
    "bulk tool output must not be shown in a full-text label"
  # The text is held in the disclosure result view, revealed on expand.
  disc.header_button.click()
  assert not disc.body.isHidden(), "click expands the disclosure"
  assert marker in disc.result_view.toPlainText()


# ---- ToolCallDisclosure (the widget bubbles use for tool cells) ---------
# Tests live here because every tool cell rendered inside a MessageBubble
# is a ToolCallDisclosure; pinning the disclosure contract alongside the
# bubble tests keeps the assertions on the same surface together.


def _qapp():
  from qttbx.qt import QtWidgets
  return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def exercise_disclosure_starts_collapsed_with_running_status():
  _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  w = ToolCallDisclosure(name="phenix_start_job", status="running")
  assert w.header_button.text().startswith("▸ ")
  assert "phenix_start_job" in w.header_button.text()
  assert "running" in w.header_button.text()
  assert w.body.isHidden()


def exercise_disclosure_click_expands_and_re_click_collapses():
  _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  w = ToolCallDisclosure(name="phenix_get_log_summary", status="running")
  w.header_button.click()
  assert not w.body.isHidden()
  assert w.header_button.text().startswith("▾ ")
  w.header_button.click()
  assert w.body.isHidden()
  assert w.header_button.text().startswith("▸ ")


def exercise_disclosure_set_status_updates_header_text_and_color():
  _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  w = ToolCallDisclosure(name="phenix_refine", status="running")
  w.set_status("finished, 12s")
  assert "finished, 12s" in w.header_button.text()
  ss_running = w.header_button.styleSheet()
  w.set_status("failed: bad input", color="error")
  ss_failed = w.header_button.styleSheet()
  # Color must change between running / finished / failed so the user
  # sees state transitions even when the row stays collapsed.
  assert ss_failed != ss_running, (ss_failed, ss_running)


def exercise_disclosure_set_args_and_result_populates_body():
  _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  w = ToolCallDisclosure(name="phenix_start_job", status="finished")
  w.set_args({"program": "phenix.refine", "args": ["model.pdb"]})
  w.set_result("R-free: 0.32 -> 0.27 over 5 macrocycles.")
  # Expand to inspect body text.
  w.header_button.click()
  assert "phenix.refine" in w.args_view.toPlainText()
  assert "R-free" in w.result_view.toPlainText()


def exercise_disclosure_expanded_result_view_fits_all_content():
  """Per-bubble scrollbars are disabled (the outer ConversationView is
  the sole scroller), so when the user clicks to expand the body the
  result_view MUST be tall enough to render every line of output -- a
  capped height clips content the user can never reveal otherwise.

  Regression: QPlainTextDocumentLayout's documentSize().height()
  returned the BLOCK count instead of pixels, so auto_height capped
  the widget at ~one line regardless of how many lines of bash output
  it held."""
  from qttbx.qt import QtWidgets
  app = _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  host = QtWidgets.QWidget()
  host.resize(600, 600)
  v = QtWidgets.QVBoxLayout(host)
  w = ToolCallDisclosure(name="Bash", status="finished", parent=host)
  v.addWidget(w)
  v.addStretch(1)
  host.show()
  n_lines = 30
  long_result = "\n".join(
    "line %02d of bash output" % i for i in range(n_lines))
  w.set_args({"command": "ls -la", "timeout": 30})
  w.set_result(long_result)
  w.header_button.click()
  for _ in range(5):
    app.processEvents()
  line_h = w.result_view.fontMetrics().lineSpacing()
  # Allow ~2 lines of slack for frame chrome / margins; the height
  # MUST cover at least N - 2 lines or the user can't read the tail.
  needed = (n_lines - 2) * line_h
  assert w.result_view.height() >= needed, (
    w.result_view.height(), needed, line_h, n_lines)


def exercise_disclosure_no_hardcoded_color_on_args_view():
  """Args view must NOT carry a hardcoded text color. Theme palettes
  (dark mode) make a fixed mid-grey unreadable; the body should
  inherit the platform's text color instead."""
  _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  w = ToolCallDisclosure(name="ToolSearch", status="finished")
  ss = w.args_view.styleSheet() or ""
  # The earlier hardcoded color was '#555'; any explicit 'color:' rule
  # on the args view would break dark-theme readability the same way.
  assert "color:" not in ss.replace(" ", "").lower(), ss


def exercise_disclosure_is_running_reflects_status():
  """``is_running`` is the predicate the cancel sweep uses to find tool cells
  that never reached a terminal state -- True for the initial 'running' status,
  False once a finished/failed/cancelled status is set."""
  _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  w = ToolCallDisclosure(name="phenix_start_job", status="running")
  assert w.is_running()
  w.set_status("finished, 3s", color="default")
  assert not w.is_running()


def exercise_disclosure_cancelled_state_distinct_from_finished_and_failed():
  """The cancelled terminal state must be visually distinct from BOTH the
  success (finished) and failure (failed) states, so a tool aborted when the
  user hits Stop reads differently from one that completed or errored."""
  _qapp()
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  finished = ToolCallDisclosure(name="t", status="running")
  finished.set_status("finished", color="default")
  failed = ToolCallDisclosure(name="t", status="running")
  failed.set_status("failed: boom", color="error")
  cancelled = ToolCallDisclosure(name="t", status="running")
  cancelled.set_status("cancelled", color="cancelled")
  assert "cancelled" in cancelled.header_button.text()
  ss_cancelled = cancelled.header_button.styleSheet()
  assert ss_cancelled != finished.header_button.styleSheet(), ss_cancelled
  assert ss_cancelled != failed.header_button.styleSheet(), ss_cancelled


def exercise_tool_cell_cancelled_marks_terminal_state():
  """``set_tool_use_finished(cancelled=True)`` transitions a running tool cell
  to the cancelled terminal state: no longer running, header shows
  'cancelled' (distinct from the finished/failed wording)."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.widgets.chat.message_bubble import MessageBubble
  bubble = MessageBubble(role="assistant")
  cell = bubble.add_tool_use_cell(
    tool_id="t1", name="phenix_start_job", args={})
  assert cell.is_running()
  bubble.set_tool_use_finished(tool_id="t1", cancelled=True)
  assert not cell.is_running()
  assert "cancelled" in cell.header_button.text()


def exercise_cancel_running_tools_sweeps_only_running_cells():
  """``cancel_running_tools()`` marks every still-running tool cell cancelled
  and leaves already-finished cells untouched -- the finalize-on-Stop sweep."""
  os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
  _qapp()
  from qttbx.widgets.chat.message_bubble import MessageBubble
  bubble = MessageBubble(role="assistant")
  done = bubble.add_tool_use_cell(tool_id="done", name="a", args={})
  running = bubble.add_tool_use_cell(tool_id="run", name="b", args={})
  bubble.set_tool_use_finished(tool_id="done", result="ok")
  assert not done.is_running()
  assert running.is_running()
  bubble.cancel_running_tools()
  assert not running.is_running()
  assert "cancelled" in running.header_button.text()
  # The already-finished cell is left exactly as it was.
  assert "finished" in done.header_button.text()
  assert "cancelled" not in done.header_button.text()


def exercise_untrusted_text_labels_use_plain_text_format():
  """Thinking text, image captions, and tool-result text are
  model/tool-controlled, so an embedded ``<img src="file://...">`` must
  be shown literally rather than rendered as rich text (which would load
  the local file at paint time). The caption renders in a PlainText
  QLabel; thinking and tool-result text render in QPlainTextEdits, which
  have no rich-text path at all."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.message_bubble import MessageBubble
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  evil = '<img src="file:///etc/passwd">'
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="thinking", data={"text": evil}),
    ContentBlock(type="tool_result", data={
      "tool_use_id": "", "is_error": False,
      "content": [ContentBlock(type="text", data={"text": evil})]}),
    ContentBlock(type="image", data={
      "attachment_sha256": "", "mime": "image/png", "caption": evil})])
  b = MessageBubble(m)
  # Image caption carries the literal text in a PlainText QLabel.
  hits = [lbl for lbl in b.findChildren(QtWidgets.QLabel)
          if "<img" in lbl.text()]
  assert len(hits) == 1, [lbl.text()
                          for lbl in b.findChildren(QtWidgets.QLabel)]
  for lbl in hits:
    assert lbl.textFormat() == QtCore.Qt.PlainText, lbl.text()
  # Thinking text is carried by the thinking cell's read-only
  # QPlainTextEdit -- plain text only, no rich-text path.
  assert b._thinking_cell is not None
  assert evil in b._thinking_cell.view.toPlainText()
  assert b._thinking_cell.view.isReadOnly()
  # Tool-result text is likewise a QPlainTextEdit.
  discs = b.findChildren(ToolCallDisclosure)
  assert len(discs) == 1, discs
  assert evil in discs[0].result_view.toPlainText()


def exercise_thinking_cell_append_streams_and_survives_highlights():
  """append() extends the plain text in place, and appending while extra
  selections (search highlights) are applied neither corrupts the text
  nor raises -- selections are a paint overlay, not document content."""
  from qttbx.qt import QtGui
  from qttbx.widgets.chat.message_bubble import _ThinkingCell
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  cell = _ThinkingCell("alpha")
  assert cell.view.toPlainText() == "[thinking] alpha"
  cell.append(" beta")
  assert cell.view.toPlainText() == "[thinking] alpha beta"
  # Visual-parity pins: frameless, no doc margin, italic font,
  # transparent background, auto-height applied (no inner scrollbars,
  # refresh hook installed), and mouse-selectable text.
  from qttbx.qt import QtCore
  assert cell.view.frameShape() == QtWidgets.QFrame.NoFrame
  assert cell.view.document().documentMargin() == 0
  assert cell.view.font().italic()
  assert "background: transparent" in cell.view.styleSheet()
  assert cell.view.verticalScrollBarPolicy() == QtCore.Qt.ScrollBarAlwaysOff
  assert hasattr(cell.view, "_auto_height_refresh")
  assert cell.view.textInteractionFlags() & QtCore.Qt.TextSelectableByMouse
  # Apply a highlight-style extra selection, then keep streaming.
  sel = QtWidgets.QTextEdit.ExtraSelection()
  cursor = QtGui.QTextCursor(cell.view.document())
  cursor.setPosition(0)
  cursor.setPosition(5, QtGui.QTextCursor.KeepAnchor)
  sel.cursor = cursor
  sel.format.setBackground(QtGui.QColor("#FFF176"))
  cell.view.setExtraSelections([sel])
  cell.append(" gamma")
  assert cell.view.toPlainText() == "[thinking] alpha beta gamma"


def exercise_short_json_shared_from_tool_approval_not_redefined():
  """_short_json has one definition (in tool_approval); message_bubble shares
  it rather than carrying a byte-identical copy. Pin that no module-level copy
  comes back AND that combined_text still truncates tool_use input via the
  shared helper."""
  import qttbx.widgets.chat.message_bubble as mb
  _qapp()
  assert not hasattr(mb, "_short_json"), \
    "message_bubble must not re-introduce a local _short_json copy"
  m = Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="tool_use", data={
      "id": "t1", "name": "phenix_start_job",
      "input": {"k": "v" * 200}})])
  text = mb.MessageBubble(m).combined_text()
  assert "phenix_start_job" in text
  assert "..." in text                  # long input truncated by _short_json


def exercise():
  exercise_renders_user_text()
  exercise_renders_tool_use_cell()
  exercise_renders_thinking_block()
  exercise_image_block_renders_as_placeholder()
  exercise_untrusted_text_labels_use_plain_text_format()
  exercise_thinking_cell_append_streams_and_survives_highlights()
  exercise_streaming_append()
  exercise_thinking_delta_after_text_starts_new_block()
  exercise_image_cell_renders_real_image()
  exercise_image_cell_click_emits_signal()
  exercise_bubble_has_no_frame_border()
  exercise_role_renders_as_bold_prefix_on_first_text_line()
  exercise_assistant_label_from_backend_stamp()
  exercise_assistant_label_falls_back_to_passed_name_then_assistant()
  exercise_tool_use_cell_uses_disclosure_widget()
  exercise_text_after_tool_cell_does_not_merge_into_prior_text_view()
  exercise_image_cell_renders_inline_thumbnail_with_click_to_lightbox()
  exercise_image_cell_without_caption_hides_caption_label()
  exercise_server_tool_use_then_result_folds_into_same_disclosure()
  exercise_orphan_tool_result_renders_collapsed_disclosure()
  exercise_disclosure_starts_collapsed_with_running_status()
  exercise_disclosure_click_expands_and_re_click_collapses()
  exercise_disclosure_set_status_updates_header_text_and_color()
  exercise_disclosure_set_args_and_result_populates_body()
  exercise_disclosure_expanded_result_view_fits_all_content()
  exercise_disclosure_no_hardcoded_color_on_args_view()
  exercise_disclosure_is_running_reflects_status()
  exercise_disclosure_cancelled_state_distinct_from_finished_and_failed()
  exercise_tool_cell_cancelled_marks_terminal_state()
  exercise_cancel_running_tools_sweeps_only_running_cells()
  exercise_short_json_shared_from_tool_approval_not_redefined()
  exercise_failed_tool_error_goes_to_the_body_not_the_header()
  exercise_failed_tool_header_status_is_clamped_to_one_short_line()
  exercise_long_tool_name_alone_does_not_floor_the_header_width()
  exercise_tool_header_hugs_its_text_instead_of_spanning_the_bubble()
  exercise_failed_tool_shows_the_error_even_when_a_result_is_present()
  exercise_non_str_tool_error_renders_instead_of_crashing()
  exercise_non_str_status_renders_instead_of_crashing()


def exercise_failed_tool_error_goes_to_the_body_not_the_header():
  """A failed tool's error text is shown in the disclosure body, not inlined
  into the header.

  The header is a QToolButton: it never elides and its minimumSizeHint is its
  full text width, so an inlined error floors the bubble's -- and hence the
  whole ConversationView's -- minimum width. The body's result view already
  wraps and auto-heights, and is where bulk tool output belongs. The header
  still says 'failed' so the state stays readable while collapsed.
  """
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  b = MessageBubble(role="assistant")
  b.add_tool_use_cell(tool_id="t1", name="mcp__coot__run_python", args={})
  b.set_tool_use_finished(tool_id="t1", error=LONG_TOOL_ERROR)
  cell = b._tool_cells_by_id["t1"]
  header = cell.header_button.text()
  assert "failed" in header, header
  assert LONG_TOOL_ERROR not in header, "full error inlined into the header"
  assert LONG_TOOL_ERROR in cell.result_view.toPlainText(), \
    "error text was dropped instead of being shown in the body"
  # The floor this bug was made of: the header must stay narrow enough that a
  # bubble can still shrink with the window.
  assert cell.header_button.minimumSizeHint().width() < 400, \
    cell.header_button.minimumSizeHint().width()


def exercise_failed_tool_header_status_is_clamped_to_one_short_line():
  """Defense in depth for the header floor: whatever a caller passes as a
  status, the rendered header stays a short single line (full text on the
  tooltip). Without this, any future long status would re-floor the layout the
  way the inlined error did."""
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  cell = ToolCallDisclosure(name="Bash", status="running")
  cell.set_status("failed: line one\nline two\n" + "x" * 500, color="error")
  header = cell.header_button.text()
  assert "\n" not in header, "header must stay a single line"
  assert len(header) < 120, len(header)
  assert cell.header_button.minimumSizeHint().width() < 400, \
    cell.header_button.minimumSizeHint().width()
  assert "line two" in cell.header_button.toolTip(), \
    "full status should remain available on the tooltip"


def exercise_long_tool_name_alone_does_not_floor_the_header_width():
  """The tool name reaches the header with no clamp of any kind.

  Nothing bounds a server-registered tool name, and add_tool_use_cell renders
  it straight into the QToolButton -- so a long one floors the view on its own,
  before any status is considered, and does so again on every reload.
  """
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  cell = ToolCallDisclosure(
    name="mcp__structure_tools__" + "compute_real_space_correlation" * 3,
    status="running")
  assert cell.header_button.minimumSizeHint().width() < 400, (
    "header minimum width %d is floored by the tool name alone"
    % cell.header_button.minimumSizeHint().width())


def exercise_tool_header_hugs_its_text_instead_of_spanning_the_bubble():
  """The disclosure header is a button, not a banner.

  QToolButton defaults to a Fixed horizontal policy exactly so that it hugs
  its text: with autoRaise on, the hover highlight then covers the row you can
  actually click and nothing more. Making the header shrinkable (so a long
  tool name cannot floor the view) must not also make it growable -- that
  stretches it to the full bubble width and lights the whole bubble up on
  hover.
  """
  from qttbx.widgets.chat.message_bubble import MessageBubble
  from qttbx.widgets.chat.conversation_view import ConversationView
  from qttbx.widgets.chat.agent.conversation import Message
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  view = ConversationView()
  bub = view.add_message(Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "Fitting the ligand now."})]))
  bub.add_tool_use_cell(tool_id="t1", name="coot_ping", args={})
  bub.set_tool_use_finished(tool_id="t1", result="pong", elapsed="0.1s")
  view.resize(900, 600)
  view.show()
  app.processEvents()
  header = bub._tool_cells_by_id["t1"].header_button
  assert header.width() <= header.sizeHint().width() + 2, (
    "header grew to %d px past its %d px of text"
    % (header.width(), header.sizeHint().width()))
  assert header.width() < bub.width(), (
    "header spans the whole %d px bubble instead of hugging its text"
    % bub.width())
  del MessageBubble


def exercise_failed_tool_shows_the_error_even_when_a_result_is_present():
  """Passing both result and error must not silently discard the error.

  Both are documented parameters of set_tool_use_finished, and the old header
  showed 'failed: <error>' regardless of result. Routing the error into the
  body only when result is None loses it entirely for a caller that has both:
  the header says a bare 'failed' and the error appears nowhere at all.
  """
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  b = MessageBubble(role="assistant")
  b.add_tool_use_cell(tool_id="t1", name="mcp__coot__run_python", args={})
  b.set_tool_use_finished(
    tool_id="t1", result="partial output before the bridge died",
    error=LONG_TOOL_ERROR)
  cell = b._tool_cells_by_id["t1"]
  body = cell.result_view.toPlainText()
  assert "partial output before the bridge died" in body, body
  assert LONG_TOOL_ERROR in body, (
    "the error was discarded because a result was also passed")


def exercise_non_str_tool_error_renders_instead_of_crashing():
  """A non-str error must render, not raise.

  The docstring asks only for an 'error message', and the old header coerced
  whatever it got via '%s'. Routing the value into set_result drops that
  coercion, so an exception object reaches QPlainTextEdit.setPlainText and
  raises TypeError inside the handler processing the tool result -- the cell
  never reaches a terminal state and the traceback surfaces mid-turn.
  """
  from qttbx.widgets.chat.message_bubble import MessageBubble
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  b = MessageBubble(role="assistant")
  b.add_tool_use_cell(tool_id="t1", name="mcp__coot__run_python", args={})
  b.set_tool_use_finished(tool_id="t1", error=RuntimeError("bridge is gone"))
  cell = b._tool_cells_by_id["t1"]
  assert "bridge is gone" in cell.result_view.toPlainText(), \
    cell.result_view.toPlainText()
  assert "failed" in cell.header_button.text(), cell.header_button.text()


def exercise_non_str_status_renders_instead_of_crashing():
  """Same coercion drop, second site: set_status on a non-str value.

  The old header build coerced any object via '%s'; splitting the raw status
  to clamp it calls .split() on whatever it is, so an exception object raises
  AttributeError out of set_status instead of rendering.
  """
  from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  cell = ToolCallDisclosure(name="Bash", status="running")
  cell.set_status(ValueError("boom"), color="error")
  assert "boom" in cell.header_button.text(), cell.header_button.text()


def _png_bytes(width=10, height=10):
  from qttbx.qt import QtCore, QtGui
  img = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
  img.fill(QtGui.QColor(0, 128, 255))
  buf = QtCore.QBuffer()
  buf.open(QtCore.QBuffer.WriteOnly)
  img.save(buf, "PNG")
  return bytes(buf.data())


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
