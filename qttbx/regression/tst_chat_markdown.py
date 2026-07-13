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


def exercise_append_markdown():
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  v = MarkdownView()
  v.append_markdown("**hello**")
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
  v.append_markdown("x")
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


def exercise_raw_html_in_markdown_is_not_rendered_as_rich_text():
  """Assistant/tool text is untrusted. Embedded raw HTML (e.g.
  ``<img src="file:///etc/passwd">``) must NOT be parsed into live rich
  text -- otherwise QTextBrowser loads the local resource at paint time.
  With HTML disabled the tag survives as literal text instead."""
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  v = MarkdownView()
  v.append_markdown('look <img src="file:///etc/passwd"> here')
  assert "<img" in v.toPlainText(), repr(v.toPlainText())


def exercise_loadresource_refuses_local_file():
  """Even a markdown image (``![](file:///...)``) -- which is plain
  markdown, not raw HTML -- must not read a local file. loadResource is
  overridden to refuse every external resource; inline chat images go
  through ImageCell, not this view."""
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.markdown_view import MarkdownView
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  v = MarkdownView()
  d = tempfile.mkdtemp()
  try:
    secret = os.path.join(d, "secret.txt")
    with open(secret, "w") as fh:
      fh.write("TOP-SECRET")
    res = v.loadResource(int(QtGui.QTextDocument.ImageResource),
                         QtCore.QUrl.fromLocalFile(secret))
    assert not res, res
    data = bytes(res) if res else b""
    assert b"TOP-SECRET" not in data, data
  finally:
    shutil.rmtree(d)


def exercise_escape_approx_tildes_is_narrow():
  """The tilde filter is deliberately narrow: only a lone '~' immediately
  before a digit (the '~2.5' "approximately" shorthand) is escaped. A real
  '~~word~~' / '~word~', a '~/path', and a bare '~' are all left untouched,
  so deliberate strikethrough and ordinary tildes render unchanged."""
  from qttbx.widgets.chat.markdown_view import _escape_approx_tildes as esc
  assert esc('~0.02 and ~0.001') == r'\~0.02 and \~0.001'
  assert esc('resolution ~2.5 A') == r'resolution \~2.5 A'
  assert esc('drop ~~this~~ now') == 'drop ~~this~~ now'    # double tilde
  assert esc('a ~word~ here') == 'a ~word~ here'            # '~' before letter
  assert esc('path ~/data/x.pdb') == 'path ~/data/x.pdb'   # '~' before slash
  assert esc('no tildes here') == 'no tildes here'


def exercise_escape_approx_tildes_skips_code():
  """The tilde filter must leave code alone: a '~<digit>' inside an inline
  `code` span or a fenced code block is literal code (a version like ~1.2.0,
  a size like ~5GB), where GFM strikethrough never applies and a backslash
  would show up verbatim. Only prose '~<digit>' is still escaped."""
  from qttbx.widgets.chat.markdown_view import _escape_approx_tildes as esc
  # Inline code span: the tilde stays literal (no spurious backslash).
  assert esc('pin `~1.2.0` please') == 'pin `~1.2.0` please'
  # Fenced code block: the tilde stays literal inside the fence.
  fenced = 'run:\n```\nsize = ~5GB\n```\ndone'
  assert esc(fenced) == fenced
  # Prose is still escaped -- the strikethrough guard is not regressed.
  assert esc('about ~5 apples') == r'about \~5 apples'
  assert esc('costs ~0.02 and ~0.001') == r'costs \~0.02 and \~0.001'


def exercise_approx_tilde_not_rendered_as_strikethrough():
  """End-to-end: two '~<number>' approximate values in one block (the second
  closable because it follows a quote) would be paired by Qt's GFM parser and
  struck. append_markdown escapes them first, so nothing is struck and the
  tildes survive as literal text. The precondition check pins that the raw
  input genuinely strikes, so the test fails if the filter is removed."""
  from qttbx.qt import QtGui
  from qttbx.widgets.chat.markdown_view import MarkdownView, _MD_FEATURES
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  from qttbx.widgets.font_init import init_default_app_font
  init_default_app_font(app)
  raw = 'costs ~0.02 rather than the "~0.001" I estimated'
  probe = QtGui.QTextDocument()
  probe.setMarkdown(raw, _MD_FEATURES)
  assert "line-through" in probe.toHtml().lower(), "input no longer strikes"
  v = MarkdownView()
  v.append_markdown(raw)
  assert "line-through" not in v.document().toHtml().lower()
  plain = v.toPlainText()
  assert "~0.02" in plain and "~0.001" in plain, repr(plain)


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


def exercise_assistant_label_reflects_backend_stamp():
  """The assistant section is labelled by the backend that produced the turn
  (google -> '## Gemini'); an unstamped message falls back to '## Assistant'."""
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="assistant", timestamp=now(), backend="google",
    content=[ContentBlock(type="text", data={"text": "rendered a plot"})]))
  conv.append(Message(role="assistant", timestamp=now(),
    content=[ContentBlock(type="text", data={"text": "no stamp here"})]))
  md = conversation_to_markdown(conv)
  assert "## Gemini\n\nrendered a plot\n" in md, md
  assert "## Assistant\n\nno stamp here\n" in md, md


def exercise_user_and_assistant_text_blocks_alternate():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="user", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "hello claude"})]))
  conv.append(Message(role="assistant", timestamp=now(), backend="anthropic",
    content=[ContentBlock(type="text", data={"text": "hi! how can i help?"})]))
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


def exercise_tool_result_fence_outlasts_inner_backticks():
  """A tool_result whose own text contains a ``` code fence must not break out
  of the export's code fence. Per CommonMark the exporter opens with a fence
  one backtick longer than the longest inner backtick run (never fewer than 3),
  so the inner ``` stays inside the block and the prose after it is not
  promoted to live markdown."""
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  inner = "log start\n```\ndef f(): return 1\n```\nAFTER_FENCE: done"
  conv.append(Message(role="user", timestamp=now(), content=[
    ContentBlock(type="tool_result", data={
      "tool_use_id": "tu_1",
      "content": [ContentBlock(type="text", data={"text": inner})],
      "is_error": False})]))
  md = conversation_to_markdown(conv)
  # The longest inner run is 3 backticks, so the fence is 4; the whole tool
  # output (its own ``` fences and the line after them) is preserved verbatim
  # between the 4-backtick open/close rather than breaking out of the fence.
  assert "````tool-result\n%s\n````" % inner in md, md


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


def exercise_image_resolves_via_public_conv_dir_accessor():
  """The export must resolve attachment paths through storage's PUBLIC
  conv_dir() accessor, not the private _conv_dir(). A storage-like object that
  exposes only the public method still resolves the image link; reverting to
  the private name AttributeErrors -> _resolve_attachment_path returns None ->
  the export falls back to the placeholder, failing the link assert below."""
  from pathlib import Path
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  tmp = tempfile.mkdtemp()
  try:
    sha = "abc123" + "0" * 58            # plain hex; only used as a filename token
    att_dir = os.path.join(tmp, "conv", "attachments")
    os.makedirs(att_dir)
    img_path = os.path.join(att_dir, "sha256-%s.png" % sha)
    with open(img_path, "wb") as fh:
      fh.write(b"PNG")

    class PublicOnlyStorage(object):
      # Exposes ONLY the public conv_dir(); it has no _conv_dir, so reaching
      # into the private name would AttributeError.
      def conv_dir(self, conv_id):
        return Path(tmp) / "conv"

    conv.append(Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="image", data={
        "attachment_sha256": sha, "mime": "image/png",
        "caption": "green square"})]))
    md = conversation_to_markdown(conv, storage=PublicOnlyStorage())
    assert "![green square](%s)" % img_path in md, md
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


def exercise_tool_result_non_iterable_content_does_not_sink_export():
  """The module docstring promises a malformed block falls back to a
  placeholder so a single bad block doesn't sink the export. A tool_result
  whose `content` is a truthy NON-iterable (e.g. a JSON number) makes
  `for inner in content` raise TypeError in _render_tool_result; without a
  per-block guard in conversation_to_markdown that TypeError aborts the whole
  'Save chat'. The export must catch it, emit a placeholder for the bad block,
  and still render the good blocks on either side."""
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="user", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "before bad block"}),
    ContentBlock(type="tool_result", data={
      "tool_use_id": "tu_1",
      "content": 42,                 # truthy non-iterable -> TypeError in loop
      "is_error": False}),
    ContentBlock(type="text", data={"text": "after bad block"})]))
  md = conversation_to_markdown(conv)              # must NOT raise
  assert "before bad block" in md, md
  assert "after bad block" in md, md
  assert "[unknown block]" in md, md


def exercise_ends_with_newline():
  from qttbx.widgets.chat.markdown_export import conversation_to_markdown
  conv, Message, ContentBlock, now = _make_conv()
  conv.append(Message(role="user", timestamp=now(), content=[
    ContentBlock(type="text", data={"text": "hi"})]))
  md = conversation_to_markdown(conv)
  assert md.endswith("\n"), repr(md[-10:])
  assert not md.endswith("\n\n\n"), repr(md[-10:])


def exercise():
  exercise_append_markdown()
  exercise_clear()
  exercise_auto_height_no_scrollbar()
  exercise_raw_html_in_markdown_is_not_rendered_as_rich_text()
  exercise_loadresource_refuses_local_file()
  exercise_escape_approx_tildes_is_narrow()
  exercise_escape_approx_tildes_skips_code()
  exercise_approx_tilde_not_rendered_as_strikethrough()
  exercise_header_renders_title_meta_and_separator()
  exercise_assistant_label_reflects_backend_stamp()
  exercise_user_and_assistant_text_blocks_alternate()
  exercise_thinking_blocks_are_skipped()
  exercise_tool_use_renders_as_fenced_block()
  exercise_tool_result_renders_text_content_in_fence()
  exercise_tool_result_fence_outlasts_inner_backticks()
  exercise_image_renders_link_to_attachment_when_storage_present()
  exercise_image_resolves_via_public_conv_dir_accessor()
  exercise_image_without_storage_falls_back_to_placeholder()
  exercise_malformed_block_does_not_crash()
  exercise_tool_result_non_iterable_content_does_not_sink_export()
  exercise_ends_with_newline()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
