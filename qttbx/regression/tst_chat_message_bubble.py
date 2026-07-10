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
  bubble.add_text("refine 1yjp")
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
  bubble.add_text("I will run phenix.refine.")
  text = bubble.first_text_cell_html()
  assert "GPT" in text and "phenix.refine" in text, text
  assert "Claude" not in text, text


def exercise_assistant_label_falls_back_to_passed_name_then_assistant():
  """With no backend stamp the bubble uses the assistant_label passed in (the
  current backend's name); with neither it is the generic 'Assistant'."""
  from qttbx.widgets.chat.message_bubble import MessageBubble
  _qapp()
  b1 = MessageBubble(role="assistant", assistant_label="Gemini")
  b1.add_text("hi")
  assert "Gemini" in b1.first_text_cell_html()
  b2 = MessageBubble(role="assistant")
  b2.add_text("hi")
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
  bubble.add_text("I'll run the job.")
  first_view = bubble._text_view
  assert first_view is not None
  bubble.add_tool_use_cell(tool_id="t1", name="phenix_start_job", args={})
  # After a non-text widget, the bubble must drop its reference to
  # _text_view so the NEXT add_text() creates a new MarkdownView.
  assert bubble._text_view is None, \
    "tool cell insertion must reset _text_view"
  bubble.add_text("Now waiting for it to finish.")
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
  the local file at paint time). Thinking and caption render in PlainText
  QLabels; tool-result text renders in the disclosure's QPlainTextEdit,
  which has no rich-text path at all."""
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
  # Thinking cell and image caption carry the literal text in PlainText
  # QLabels.
  hits = [lbl for lbl in b.findChildren(QtWidgets.QLabel)
          if "<img" in lbl.text()]
  assert len(hits) >= 2, [lbl.text() for lbl in b.findChildren(QtWidgets.QLabel)]
  for lbl in hits:
    assert lbl.textFormat() == QtCore.Qt.PlainText, lbl.text()
  # Tool-result text is carried by the disclosure's QPlainTextEdit, which
  # renders only plain text -- there is no HTML/rich-text path through
  # which the embedded file:// reference could be loaded.
  discs = b.findChildren(ToolCallDisclosure)
  assert len(discs) == 1, discs
  assert evil in discs[0].result_view.toPlainText()


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
