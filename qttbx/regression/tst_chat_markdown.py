"""Markdown rendering tests: MarkdownView (the read-only widget that
streams assistant text into a QTextBrowser) and conversation_to_markdown
(the 'Save chat' export). Both are about turning conversation content
into markdown surfaces, just at different ends of the pipeline."""

import os
import shutil
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times, null_out

try:
  from qttbx.qt import QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)


# ---- MarkdownView widget -------------------------------------------------


def exercise_set_and_append_markdown():
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  v = MarkdownView()
  v.set_markdown("**hello**")
  assert "hello" in v.toPlainText()
  v.append_markdown("\n\nworld")
  text = v.toPlainText()
  assert "hello" in text and "world" in text


def exercise_clear():
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  v = MarkdownView()
  v.set_markdown("x")
  v.clear()
  assert v.toPlainText().strip() == ""


def exercise_auto_height_no_scrollbar():
  """After the chat UI redesign, MarkdownView has no scrollbar of its
  own (auto_height). The outer ConversationView is the sole scroller.
  Pins the contract so a regression that re-enables the widget's own
  scrollbar surfaces here instead of as a stacked-scrolls UX bug."""
  from qttbx.qt import QtCore
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  v = MarkdownView()
  assert v.verticalScrollBarPolicy() == QtCore.Qt.ScrollBarAlwaysOff
  assert v.horizontalScrollBarPolicy() == QtCore.Qt.ScrollBarAlwaysOff


# ---- conversation_to_markdown export -------------------------------------


def _make_conv(title="My chat"):
  from qttbx.widgets.chat.agent.conversation import (
    Conversation, Message, ContentBlock, now)
  conv = Conversation.new(profile_name="default", model="claude-opus-4-7",
                          title=title)
  # Pin created_at so the header asserts deterministically.
  from datetime import datetime, timezone
  conv.meta.created_at = datetime(2026, 5, 18, 14, 30, tzinfo=timezone.utc)
  return conv, Message, ContentBlock, now


def exercise_header_renders_title_meta_and_separator():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, _, _, _ = _make_conv("Refining 1yjp")
  md = conversation_to_markdown(conv)
  assert md.startswith("# Refining 1yjp\n"), md
  assert "2026-05-18 14:30 UTC" in md, md
  assert "model: claude-opus-4-7" in md, md
  assert "profile: default" in md, md
  assert "\n---\n" in md, md


def exercise_user_and_assistant_text_blocks_alternate():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="user", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "hello claude"})]))
  conv.append(Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "hi! how can i help?"})]))
  md = conversation_to_markdown(conv)
  assert "## You\n\nhello claude\n" in md, md
  assert "## Claude\n\nhi! how can i help?\n" in md, md


def exercise_thinking_blocks_are_skipped():
  """Extended-thinking is internal -- exporting it would leak noise
  into the user's archived chat."""
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="thinking", data={
      "text": "let me reason carefully", "signature": "sig"}),
    ContentBlock(type="text", data={"text": "the answer is 42"})]))
  md = conversation_to_markdown(conv)
  assert "let me reason" not in md, md
  assert "the answer is 42" in md, md


def exercise_tool_use_renders_as_fenced_block():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="tool_use", data={
      "id": "tu_1", "name": "phenix_get_status",
      "input": {"job_id": "abc-123"}})]))
  md = conversation_to_markdown(conv)
  assert "```tool-use phenix_get_status" in md, md
  assert "abc-123" in md, md


def exercise_tool_result_renders_text_content_in_fence():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="user", timestamp=now(), content=[
    ContentBlock(type="tool_result", data={
      "tool_use_id": "tu_1",
      "content": [ContentBlock(type="text", data={
        "text": "status: running\nlast: 3 min ago"})],
      "is_error": False})]))
  md = conversation_to_markdown(conv)
  assert "```tool-result" in md, md
  assert "status: running" in md, md
  assert "last: 3 min ago" in md, md


def exercise_image_renders_link_to_attachment_when_storage_present():
  from qttbx.widgets.chat.agent.storage import ConversationStorage
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  from qttbx.qt import QtCore, QtGui
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=tmp, log=null_out())
    conv, Message, ContentBlock, now = _make_conv()
    storage.save(conv)
    img = QtGui.QImage(2, 2, QtGui.QImage.Format_ARGB32)
    img.fill(QtGui.QColor(0, 255, 0))
    buf = QtCore.QBuffer()
    buf.open(QtCore.QBuffer.WriteOnly)
    img.save(buf, "PNG")
    att = storage.store_attachment(
      conv.meta.id, bytes(buf.data()), "image/png")
    conv.append(Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="image", data={
        "attachment_sha256": att.sha256, "mime": "image/png",
        "caption": "a green square"})]))
    md = conversation_to_markdown(conv, storage=storage)
    assert "![a green square](" in md, md
    needle_start = "![a green square]("
    start = md.index(needle_start) + len(needle_start)
    end = md.index(")", start)
    img_path = md[start:end]
    assert os.path.exists(img_path), (img_path, md)
    assert img_path.endswith(".png"), img_path
  finally:
    shutil.rmtree(tmp)


def exercise_image_without_storage_falls_back_to_placeholder():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="image", data={
      "attachment_sha256": "deadbeef" * 8, "mime": "image/png",
      "caption": None})]))
  md = conversation_to_markdown(conv)
  assert "[image deadbeefdead image/png]" in md, md


def exercise_malformed_block_does_not_crash():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="assistant", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "before"}),
    ContentBlock(type="surprise", data={"who": "knows"}),
    ContentBlock(type="text", data={"text": "after"})]))
  md = conversation_to_markdown(conv)
  assert "before" in md, md
  assert "after" in md, md
  assert "[unknown block: surprise]" in md, md


def exercise_ends_with_newline():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="user", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "hi"})]))
  md = conversation_to_markdown(conv)
  assert md.endswith("\n"), repr(md[-10:])
  assert not md.endswith("\n\n\n"), repr(md[-10:])


def exercise():
  exercise_set_and_append_markdown()
  exercise_clear()
  exercise_auto_height_no_scrollbar()
  exercise_header_renders_title_meta_and_separator()
  exercise_user_and_assistant_text_blocks_alternate()
  exercise_thinking_blocks_are_skipped()
  exercise_tool_use_renders_as_fenced_block()
  exercise_tool_result_renders_text_content_in_fence()
  exercise_image_renders_link_to_attachment_when_storage_present()
  exercise_image_without_storage_falls_back_to_placeholder()
  exercise_malformed_block_does_not_crash()
  exercise_ends_with_newline()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
