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


def _settle(app, card):
  """Pump events until a pass changes neither the card's size nor any
  preview's height: a preview that grows on rewrap resizes the card."""
  last = None
  for _ in range(50):
    app.processEvents()
    now = (card.size(), tuple(
      v.height() for v in card.findChildren(QtWidgets.QPlainTextEdit)))
    if now == last:
      return
    last = now
  raise AssertionError("the card did not settle: %r" % (last,))


def _title(card):
  """The single card's title labels: the "<source> requests:" lead and
  the tool name."""
  from qttbx.widgets.chat.eliding import ElidingLabel
  lead, name = card._content_widget.findChildren(ElidingLabel)
  return lead, name


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
  """Ticking "Always allow this tool" is what puts remember='tool' on the
  responses: the box's toggled signal drives the field, so click the box
  rather than the setter it is wired to."""
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req()])
  box, = card.findChildren(QtWidgets.QCheckBox)
  box.click()
  assert box.isChecked()
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_approve_all()
  assert decisions[0][0].remember == "tool"


def exercise_title_drops_the_mcp_prefix_that_repeats_the_source():
  """An MCP tool registered as mcp__<server>__<name> under source
  mcp:<server> is titled by <name> beside a "<source> requests:" lead,
  so what an elided title keeps is the tool's own name -- the verb --
  and not the server twice. The registered name stays on the tooltip.
  Any other name is shown as it is. A batched card's lines keep the
  registered name: stripped, mcp__coot__delete_all and
  mcp__phenix__delete_all were two identical lines above "Approve all".
  """
  from qttbx.widgets.chat.eliding import ElidingLabel
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  cards = []                               # keep the labels' parents alive

  def title(tool_name, source):
    card = ToolApprovalCard()
    card.set_requests([_req(tool_name=tool_name, source=source)])
    cards.append(card)
    lead, name = _title(card)
    return lead.full_text(), name.full_text(), name.toolTip()

  assert title("mcp__coot__delete_residue_range_and_refine", "mcp:coot") == (
    "mcp:coot requests:", "delete_residue_range_and_refine",
    "mcp__coot__delete_residue_range_and_refine")
  assert title("phenix_start_job", "mcp:phenix") == (
    "mcp:phenix requests:", "phenix_start_job", "")
  assert title("mcp__coot__x", "builtin") == (
    "builtin requests:", "mcp__coot__x", "")
  # A name that is nothing but the prefix is shown whole, not as "".
  assert title("mcp__coot__", "mcp:coot") == (
    "mcp:coot requests:", "mcp__coot__", "")
  # A tool the registry renamed <server>__<name> for a collision keeps
  # the prefix: it is what tells the two tools apart, and "Always allow"
  # is keyed by the registered name.
  assert title("coot__delete_all", "mcp:coot") == (
    "mcp:coot requests:", "coot__delete_all", ""), (
      title("coot__delete_all", "mcp:coot"))

  card = ToolApprovalCard()
  card.set_requests([
    _req("a", tool_name="mcp__coot__delete_all", source="mcp:coot",
         batch_id="B"),
    _req("b", tool_name="mcp__phenix__delete_all", source="mcp:phenix",
         batch_id="B")])
  lines = [w.full_text()
           for w in card._content_widget.findChildren(ElidingLabel)]
  assert lines == ["mcp__coot__delete_all", "mcp__phenix__delete_all"], lines


def exercise_long_source_gives_way_before_the_tool_name():
  """A user-configured server name floors nothing and yields first.

  The "<source> requests:" lead was a plain label, so a long source
  floored the card at its full width and left the name -- the only
  shrinkable thing in the row -- at 0 px, a lone ellipsis above
  Approve. Now the lead elides, and before the name does: the name is
  what the user decides on, and the whole source is on the lead's
  tooltip. What the lead gives up, the name gets: a row that let the
  lead elide with room still empty beside the name had the stretch on
  the wrong item. And with room for both, the name sits beside the
  lead's text, not pushed to the far end of a lead grown past it.
  """
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  long_source = "mcp:a_rather_long_user_configured_server_name"

  def card_for(source):
    card = ToolApprovalCard()
    card.set_requests([_req(tool_name="phenix_start_job", source=source)])
    return card

  def empty_room(card):
    """Width of the title row neither label is using."""
    lead, name = _title(card)
    row = card._content_widget.width()
    return row - lead.width() - name.width()

  short, long = card_for("mcp:phenix"), card_for(long_source)
  assert short.minimumSizeHint().width() == long.minimumSizeHint().width(), (
    short.minimumSizeHint(), long.minimumSizeHint())
  margins = long._layout.contentsMargins()
  chrome = margins.left() + margins.right()
  # Room for the whole name beside half the source -- from the hints,
  # not a pixel literal: the lead has to elide and the name must not,
  # on this machine's font and on the 17 px glyphs of a runner without.
  lead, name = _title(long)
  long.resize(max(long.minimumSizeHint().width(), chrome
                  + lead.sizeHint().width() // 2 + name.sizeHint().width()),
              300)
  long.show()
  _settle(app, long)
  lead, name = _title(long)
  assert name.text() == "phenix_start_job", name.text()
  assert lead.text() != lead.full_text(), lead.text()
  assert lead.toolTip() == long_source + " requests:", lead.toolTip()
  # The lead elided, so nothing but the layout's spacing is spare.
  spare = empty_room(long)
  assert spare <= 2 * long._layout.spacing(), (
    "the lead is elided with %d px of its row empty" % spare)
  long.resize(long.minimumSizeHint().width(), 300)
  _settle(app, long)
  lead, name = _title(long)
  assert name.text() != "\u2026", "the name is a lone ellipsis"
  # Room for both: the name follows the lead's text, within a spacing.
  lead, name = _title(short)
  short.resize(max(short.minimumSizeHint().width(), chrome
                   + lead.sizeHint().width() + name.sizeHint().width()) + 100,
               300)
  short.show()
  _settle(app, short)
  lead, name = _title(short)
  assert lead.text() == lead.full_text() == "mcp:phenix requests:", (
    lead.text())
  assert name.x() - lead.x() <= (
    lead.sizeHint().width() + 2 * short._layout.spacing()), (
    "the name is %d px from a lead whose text takes %d px"
    % (name.x() - lead.x(), lead.sizeHint().width()))


def exercise_the_source_kind_outlasts_any_squeeze():
  """The first characters of the source are on the card at every width.

  The source is a trust signal: the Phenix agent labels an MCP tool
  mcp:<server> so that it does not read as a trusted built-in. With the
  name the last to give way, a long tool name took the lead to a lone
  ellipsis, and elided in the middle the lead dropped the server first
  ("mc…s:"). So it elides at its end and keeps enough characters for
  "builtin" whole or "mcp:" and three of the server, hoverable for the
  rest; the name, elided in the middle, takes what is left. The card's
  own minimum width -- its buttons -- is not raised by that.
  """
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  plain = ToolApprovalCard()
  plain.set_requests([_req()])
  for source, kept in (("mcp:phenix", "mcp:phe"), ("builtin", "builtin"),
                       ("skill", "skill r")):
    card = ToolApprovalCard()
    card.set_requests([_req(tool_name="mcp__coot__" + "x" * 80,
                            source=source)])
    assert card.minimumSizeHint().width() == plain.minimumSizeHint().width()
    card.resize(card.minimumSizeHint().width(), 300)
    card.show()
    _settle(app, card)
    lead, name = _title(card)
    assert lead.text().startswith(kept), (source, lead.text())
    assert lead.text().endswith("\u2026"), (source, lead.text())
    assert lead.toolTip() == source + " requests:", lead.toolTip()
    shown = name.fontMetrics().horizontalAdvance(name.text())
    assert "\u2026" in name.text() and shown <= name.width(), (
      name.text(), shown, name.width())
    card.close()


def exercise_long_text_elides_or_wraps_instead_of_clipping():
  """At the card's own minimum width every character is still reachable.

  The argument preview wraps an unbreakable token (a sequence, a path)
  onto further lines and grows to show them all: a label that merely
  stopped flooring the card clipped the token instead, with no tooltip
  and no scrollbar. A batched card's tool names elide, with the whole
  name on the tooltip, rather than clip -- "Approve all" must not hide
  which tool runs.
  """
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.eliding import ElidingLabel
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  token = "ACGT" * 35
  card = ToolApprovalCard()
  card.set_requests([ToolApprovalRequest(
    request_id="r1", tool_name="phenix_start_job", tool_source="mcp:phenix",
    input={"sequence": token}, risk="write", batch_id=None)])
  card.resize(card.minimumSizeHint().width(), 600)
  card.show()
  _settle(app, card)
  view, = card._content_widget.findChildren(QtWidgets.QPlainTextEdit)
  assert token in view.toPlainText()
  cursor = view.textCursor()
  cursor.movePosition(QtGui.QTextCursor.End)
  end = view.cursorRect(cursor)
  assert view.document().lineCount() > 1, "the token did not wrap"
  # bottom() is the last pixel row; row height() is the first outside.
  assert end.bottom() < view.viewport().height(), (
    "the preview's last line is clipped", end, view.viewport().height())
  # And no taller than its text: the view grows with its content, it
  # does not take the card's spare height for a block of blank.
  blank = view.height() - end.bottom()
  assert blank < view.fontMetrics().lineSpacing(), (
    "%d px of blank under the preview's last line" % blank)
  # Nothing left for the preview's own (hidden) scrollbars: with word
  # wrap only, the token ran off to the right; a view a line short of
  # its text kept the last line for a wheel notch to reveal.
  _assert_nothing_to_scroll(view)
  # Nor is the preview a Tab stop, sitting invisibly ahead of Approve.
  assert view.focusPolicy() == QtCore.Qt.NoFocus, view.focusPolicy()

  card = ToolApprovalCard()
  card.set_requests([ToolApprovalRequest(
    request_id="r%d" % i, tool_name="mcp__coot__" + "x" * 80,
    tool_source="mcp:coot", input={"sequence": token}, risk="write",
    batch_id="B") for i in (1, 2)])
  card.resize(card.minimumSizeHint().width(), 600)
  card.show()
  _settle(app, card)
  names = card._content_widget.findChildren(ElidingLabel)
  assert len(names) == 2, names
  for name in names:
    shown = name.fontMetrics().horizontalAdvance(name.text())
    assert shown <= name.width(), (
      "a batched tool name is clipped: %d px of text in %d px"
      % (shown, name.width()))
    assert name.toolTip() == "mcp__coot__" + "x" * 80, name.toolTip()
  views = card._content_widget.findChildren(QtWidgets.QPlainTextEdit)
  assert len(views) == 2, views
  for name, view in zip(names, views):
    _assert_nothing_to_scroll(view)
    # One indented block per request: the preview under its name.
    assert view.x() == name.x() > 0, (view.x(), name.x())


def exercise_zoomed_preview_grows_with_its_font():
  """Ctrl+wheel over the preview zooms it, as over any read-only
  QPlainTextEdit; the height has to follow, or the grown lines sit
  behind an inner scroll range that swallows the next plain notches."""
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([ToolApprovalRequest(
    request_id="r1", tool_name="phenix_start_job", tool_source="mcp:phenix",
    input={"sequence": "ACGT" * 35}, risk="write", batch_id=None)])
  card.resize(card.minimumSizeHint().width(), 600)
  card.show()
  _settle(app, card)
  view, = card._content_widget.findChildren(QtWidgets.QPlainTextEdit)
  before = (view.font().pointSize(), view.height())
  for _ in range(3):
    app.sendEvent(view.viewport(), QtGui.QWheelEvent(
      QtCore.QPointF(5, 5), QtCore.QPointF(5, 5), QtCore.QPoint(0, 0),
      QtCore.QPoint(0, 120), QtCore.Qt.NoButton, QtCore.Qt.ControlModifier,
      QtCore.Qt.ScrollUpdate, False))
  _settle(app, card)
  assert view.font().pointSize() > before[0], (
    "the view did not zoom", view.font().pointSize(), before)
  assert view.height() > before[1], (view.height(), before)
  _assert_nothing_to_scroll(view)
  cursor = view.textCursor()
  cursor.movePosition(QtGui.QTextCursor.End)
  assert view.cursorRect(cursor).bottom() < view.viewport().height(), (
    "the zoomed preview's last line is clipped",
    view.cursorRect(cursor).bottom(), view.viewport().height())


def _assert_nothing_to_scroll(view):
  ranges = (view.verticalScrollBar().maximum(),
            view.horizontalScrollBar().maximum())
  assert ranges == (0, 0), (
    "the preview has %d px / %d lines of its own to scroll"
    % (ranges[1], ranges[0]))


def exercise_remember_checkbox_sits_under_the_buttons():
  """The checkbox has a row of its own: in the button row it was the only
  shrinkable item, so it absorbed the whole shortfall and read
  "Always all..." at the ChatWindow's narrowest column."""
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req()])
  card.resize(card.minimumSizeHint().width(), 300)
  card.show()
  _settle(app, card)
  box, = card.findChildren(QtWidgets.QCheckBox)
  buttons = card._buttons_widget.findChildren(QtWidgets.QPushButton)
  assert box.y() >= max(b.geometry().bottom() for b in buttons), (
    box.geometry(), [b.geometry() for b in buttons])
  # And it floors nothing. The view's floor guards cannot see this: they
  # compare against a view whose card carries the same checkbox.
  assert box.minimumSizeHint().width() == 0, box.minimumSizeHint()


def exercise_batch_growth_forgets_a_ticked_remember():
  """A same-batch request arriving after "Always allow this tool" was
  ticked rebuilds the card batched, without the checkbox; the tick
  must go with it, or Approve all remembers every tool in the batch."""
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req(batch_id="B")])
  box, = card.findChildren(QtWidgets.QCheckBox)
  box.click()
  card.set_requests([_req(batch_id="B"), _req("r2", batch_id="B")])
  assert not card.findChildren(QtWidgets.QCheckBox)
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  card.click_approve_all()
  assert [r.remember for r in decisions[0]] == ["none", "none"], decisions


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
  -- the very surface the user reads to decide. The argument preview is
  a plain-text view, which cannot render markup at all, and the name
  labels use PlainText, so a tool argument or name containing HTML can't
  inject rich text (e.g. a clickable file:// link) into the prompt."""
  from qttbx.qt import QtCore, QtGui
  from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  evil = "<a href='file:///etc/passwd'>ok</a>"

  def check(card, previews, labels=1):
    views = [v for v in card.findChildren(QtWidgets.QPlainTextEdit)
             if "file:///" in v.toPlainText()]
    assert len(views) == previews, [
      v.toPlainText() for v in card.findChildren(QtWidgets.QPlainTextEdit)]
    names = [l for l in card.findChildren(QtWidgets.QLabel)
             if "file:///" in l.text()]
    assert len(names) == labels, [
      l.text() for l in card.findChildren(QtWidgets.QLabel)]
    for l in names:
      assert l.textFormat() == QtCore.Qt.PlainText, l.text()

  # Single-request card: the source and name labels and the preview.
  # The source is the user's server name, from the profile: no more
  # trusted than the tool's, on the surface that says who is asking.
  card = ToolApprovalCard()
  card.set_requests([ToolApprovalRequest(
    request_id="r1", tool_name=evil, tool_source="mcp:" + evil,
    input={"path": evil}, risk="write", batch_id=None)])
  check(card, previews=1, labels=2)
  # The tooltip too: a tooltip auto-detects rich text, and a name shown
  # shortened carries the registered one there.
  registered = "mcp__phenix__" + evil
  card3 = ToolApprovalCard()
  card3.set_requests([ToolApprovalRequest(
    request_id="r1", tool_name=registered, tool_source="mcp:phenix",
    input={}, risk="write", batch_id=None)])
  lead, name = _title(card3)
  assert name.full_text() == evil, name.full_text()
  shown = QtGui.QTextDocumentFragment.fromHtml(name.toolTip()).toPlainText()
  assert shown == registered, (name.toolTip(), shown)
  # Batched card: the per-request name labels and previews.
  card2 = ToolApprovalCard()
  card2.set_requests([
    ToolApprovalRequest(request_id="a", tool_name=evil, tool_source="s",
                        input={"p": evil}, risk="write", batch_id="B"),
    ToolApprovalRequest(request_id="b", tool_name="t2", tool_source="s",
                        input={"p": "ok"}, risk="write", batch_id="B")])
  check(card2, previews=1)


def exercise_allow_remember_hides_checkbox():
  """A request with allow_remember=False yields a card without the
  'Always allow this tool' checkbox; the default (True) keeps it."""
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)

  def _checkbox_labels(card):
    # full_text(): the rendered text elides once the card is squeezed.
    return [w.full_text() for w in card.findChildren(QtWidgets.QCheckBox)]

  # Default (allow_remember=True) keeps the "Always allow this tool" checkbox.
  card_default = ToolApprovalCard()
  card_default.set_requests([ToolApprovalRequest(
    request_id="r1", tool_name="phenix_start_job", tool_source="mcp:phenix",
    input={}, risk="write", batch_id=None)])
  assert "Always allow this tool" in _checkbox_labels(card_default)

  # allow_remember=False suppresses it (a destructive tool opts out of the
  # standing per-tool auto-approval the checkbox grants).
  card_no_remember = ToolApprovalCard()
  card_no_remember.set_requests([ToolApprovalRequest(
    request_id="r2", tool_name="coot_close_unresponsive",
    tool_source="mcp:phenix_agent", input={}, risk="destructive",
    batch_id=None, allow_remember=False)])
  assert "Always allow this tool" not in _checkbox_labels(card_no_remember)


def exercise_finalize_disables_card_without_emitting():
  """finalize() disables an undecided card WITHOUT emitting a decision, so a
  later stray/queued click can't route a stale response into a subsequent turn
  (the approval-misroute guard). A no-op once a real decision has emitted."""
  from qttbx.widgets.chat.tool_approval import ToolApprovalCard
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  card = ToolApprovalCard()
  card.set_requests([_req()])
  decisions = []
  card.decided.connect(lambda resps: decisions.append(resps))
  assert not card.is_decided()
  card.finalize()
  assert card.is_decided()                       # marked decided (finalized)
  assert card.isHidden()                         # collapsed like a decided card
  assert not card._buttons_widget.isEnabled()    # buttons disabled
  # A later click on the abandoned card must NOT emit a stale response.
  card.click_stop()
  card.click_approve_all()
  assert decisions == [], decisions
  # finalize() after a real decision is a harmless no-op (no re-emit).
  card2 = ToolApprovalCard()
  card2.set_requests([_req()])
  d2 = []
  card2.decided.connect(lambda resps: d2.append(resps))
  card2.click_approve_all()
  assert len(d2) == 1
  card2.finalize()
  assert len(d2) == 1, d2


_LAYOUT_EXERCISES = (
  exercise_long_source_gives_way_before_the_tool_name,
  exercise_the_source_kind_outlasts_any_squeeze,
  exercise_long_text_elides_or_wraps_instead_of_clipping,
  exercise_remember_checkbox_sits_under_the_buttons,
)


def exercise_layout_again_under_a_font_this_machine_lacks():
  """The layout checks again under a font wide enough that any width
  still written as a pixel literal cannot be enough: the Windows CI
  runner has no font files and measures every glyph at 17 px (see
  tst_chat_eliding)."""
  from qttbx.qt import QtGui
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  original = QtGui.QFont(app.font())
  wide = QtGui.QFont(original)
  wide.setPointSize(max(24, original.pointSize() * 2))
  app.setFont(wide)
  try:
    for fn in _LAYOUT_EXERCISES:
      fn()
  finally:
    app.setFont(original)


def exercise():
  exercise_single_card_approve_emits_response()
  exercise_single_card_deny_and_stop()
  exercise_batched_card_approve_all_emits_n_responses_in_order()
  exercise_batched_deny_all()
  exercise_remember_tool_checkbox_sets_remember_field()
  exercise_title_drops_the_mcp_prefix_that_repeats_the_source()
  exercise_long_source_gives_way_before_the_tool_name()
  exercise_the_source_kind_outlasts_any_squeeze()
  exercise_long_text_elides_or_wraps_instead_of_clipping()
  exercise_zoomed_preview_grows_with_its_font()
  exercise_remember_checkbox_sits_under_the_buttons()
  exercise_batch_growth_forgets_a_ticked_remember()
  exercise_card_hides_and_disables_buttons_after_click()
  exercise_is_decided_reflects_decision_state()
  exercise_finalize_disables_card_without_emitting()
  exercise_untrusted_input_rendered_as_plain_text()
  exercise_allow_remember_hides_checkbox()
  exercise_layout_again_under_a_font_this_machine_lacks()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
