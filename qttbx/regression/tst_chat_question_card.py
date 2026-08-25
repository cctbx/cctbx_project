"""QuestionCard widget test. Covers the radio-button (single-select)
and checkbox (multi-select) modes, the 'Other' free-form override, the
emitted answer payload, and the disable-on-submit behaviour."""

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


def _qapp():
  from qttbx.widgets.font_init import init_default_app_font
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  return app


def _single_select_question(label_choice):
  return [{
    "question": "Which backend?",
    "header": "Backend",
    "multiSelect": False,
    "options": [
      {"label": "anthropic", "description": "Use the Anthropic API"},
      {"label": "claude_code", "description": "Use Claude Code OAuth"},
    ],
  }]


def exercise_single_select_emits_chosen_label():
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  card.set_request("q_1", _single_select_question("claude_code"))
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  # Pick the second radio.
  state = card._question_state[0]
  state["option_buttons"][1][0].setChecked(True)
  card.click_submit()
  assert fired == [("q_1", {"Which backend?": "claude_code"})], fired


def exercise_single_select_other_overrides_radio():
  """When the user types in 'Other', that beats any radio selection.
  Intent: whatever they last typed is what they meant."""
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  card.set_request("q_2", _single_select_question("claude_code"))
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  state = card._question_state[0]
  state["option_buttons"][0][0].setChecked(True)
  state["other_edit"].setText("custom_backend")
  card.click_submit()
  assert fired == [("q_2", {"Which backend?": "custom_backend"})], fired


def exercise_multi_select_collects_all_checked():
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  card.set_request("q_3", [{
    "question": "Which programs?",
    "multiSelect": True,
    "options": [
      {"label": "phenix.refine", "description": "Refinement"},
      {"label": "phenix.xtriage", "description": "Data quality"},
      {"label": "phenix.molprobity", "description": "Validation"},
    ],
  }])
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  state = card._question_state[0]
  state["option_buttons"][0][0].setChecked(True)
  state["option_buttons"][2][0].setChecked(True)
  card.click_submit()
  assert fired == [("q_3", {"Which programs?": [
    "phenix.refine", "phenix.molprobity"]})], fired


def exercise_multi_select_other_is_appended():
  """Multi-select 'Other' adds the typed value to the list rather than
  replacing the checked boxes."""
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  card.set_request("q_4", [{
    "question": "Which?",
    "multiSelect": True,
    "options": [
      {"label": "a"}, {"label": "b"},
    ],
  }])
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  state = card._question_state[0]
  state["option_buttons"][0][0].setChecked(True)
  state["other_edit"].setText("c")
  card.click_submit()
  assert fired == [("q_4", {"Which?": ["a", "c"]})], fired


def exercise_multiple_questions_in_one_card():
  """The model groups its questions into a single tool call -- the
  card MUST render them all and emit a dict keyed by each question."""
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  card.set_request("q_5", [
    {"question": "A?", "options": [{"label": "a1"}, {"label": "a2"}]},
    {"question": "B?", "options": [{"label": "b1"}, {"label": "b2"}]},
  ])
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  card._question_state[0]["option_buttons"][0][0].setChecked(True)
  card._question_state[1]["option_buttons"][1][0].setChecked(True)
  card.click_submit()
  assert fired == [("q_5", {"A?": "a1", "B?": "b2"})], fired


def exercise_duplicate_question_text_keeps_both_answers():
  """Two questions may carry the SAME text (or empty/missing text). Keying
  the answers dict by raw question text drops one -- the second write
  clobbers the first (answers[text]=picked). Each of the N questions must
  still yield its own entry, so no user selection is silently lost."""
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  card.set_request("q_dup", [
    {"question": "Same?", "options": [{"label": "a1"}, {"label": "a2"}]},
    {"question": "Same?", "options": [{"label": "b1"}, {"label": "b2"}]},
  ])
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  card._question_state[0]["option_buttons"][0][0].setChecked(True)
  card._question_state[1]["option_buttons"][1][0].setChecked(True)
  card.click_submit()
  assert fired, fired
  rid, ans = fired[0]
  assert rid == "q_dup", rid
  # Two questions -> two entries; neither answer overwrites the other.
  assert len(ans) == 2, ans
  picked = list(ans.values())
  assert "a1" in picked and "b2" in picked, ans


def exercise_header_uses_assistant_label():
  """The card header names the active assistant ('GPT needs an answer:'),
  not a hard-coded 'Claude', so it matches the chosen backend."""
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard(assistant_label="GPT")
  card.set_request("q_h", _single_select_question("anthropic"))
  labels = [w.text() for w in card.findChildren(QtWidgets.QLabel)]
  assert any("GPT needs an answer" in t for t in labels), labels
  assert not any("Claude needs an answer" in t for t in labels), labels


def exercise_submit_disables_buttons_and_hides_card():
  """After submit, the buttons are disabled and the card is hidden so
  a stray re-click can't re-emit. Same UX contract as
  ToolApprovalCard."""
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  card.set_request("q_6", _single_select_question("anthropic"))
  card._question_state[0]["option_buttons"][0][0].setChecked(True)
  card.click_submit()
  assert card.isHidden()
  assert not card._buttons_widget.isEnabled()


def exercise_malformed_question_payload_unblocks_worker():
  """A malformed ask_user_question payload (not a list of question dicts)
  must not crash the card. The card is built on the GUI thread while the
  agent worker is parked on question_queue.get(); a crash here strands that
  worker forever. set_request must stay quiet AND release the waiting path
  by emitting `answered` so the worker unblocks with a safe error answer."""
  from qttbx.widgets.chat.question_card import QuestionCard
  app = _qapp()
  card = QuestionCard()
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  # A bare string is the classic shape: list("...") explodes into single
  # characters, each a non-dict that crashes q.get(...) when the card
  # builds. set_request must not raise.
  card.set_request("q_bad", "totally malformed")
  # The unblock is emitted deferred (the real caller connects `answered`
  # only after set_request returns), so let the event loop run it.
  app.processEvents()
  assert fired, "answered not emitted; the parked worker would hang"
  rid, ans = fired[0]
  assert rid == "q_bad", rid
  assert isinstance(ans, dict), ans


def exercise_non_dict_options_do_not_crash():
  """An option that isn't a {label, description} dict (a bare string, or
  junk) must not crash card construction at opt.get('label'). The card
  still renders and a normal Submit emits an answer for the question."""
  from qttbx.widgets.chat.question_card import QuestionCard
  _qapp()
  card = QuestionCard()
  # Mixed junk: a plain string, a proper dict, and entries that are neither.
  card.set_request("q_opt", [{
    "question": "Pick?",
    "options": ["yes", {"label": "no"}, 5, None],
  }])
  fired = []
  card.answered.connect(lambda rid, ans: fired.append((rid, ans)))
  card.click_submit()
  assert fired and fired[0][0] == "q_opt", fired


def _card_min_width(questions):
  """Minimum width a fresh card reports for ``questions``."""
  from qttbx.widgets.chat.question_card import QuestionCard
  card = QuestionCard()
  card.set_request("r1", questions)
  return card.minimumSizeHint().width()


def exercise_unbroken_question_header_does_not_floor_the_card_width():
  """A long header chip must not widen the card, even as one unbroken token.

  'header' is model-controlled: meant to be a ~12-char chip, but nothing
  enforces that, and an over-long one floors the card's -- and hence the whole
  ConversationView's -- minimum width, stopping the bubbles tracking the
  window. The chip elides rather than wraps: word wrap would not have unfloored
  it anyway, since a wrapped QLabel reports its widest unbreakable token as its
  minimumSizeHint width, so a space-free header (a snake_case identifier, a
  path, a URL) would floor the card exactly as the whole string would.

  Measured against the same card with a conventional short header, not against
  an absolute pixel budget: the card's fixed chrome (the '<assistant> needs an
  answer:' head, the Other row, Submit) scales with the platform's font
  metrics, and a font-less Qt -- the Windows CI conda build ships no
  Library/lib/fonts, so every glyph measures as a fallback box the size of the
  font -- pushes that chrome past any absolute budget while the chip behaves
  perfectly. Only the header differs between the two payloads, so the delta is
  exactly what the chip adds.
  """
  _qapp()
  base_width = _card_min_width([{
    "question": "Which backend?",
    "header": "Backend",
    "multiSelect": False,
    "options": [{"label": "a", "description": "d"},
                {"label": "b", "description": "d"}],
  }])
  wide_width = _card_min_width([{
    "question": "Which backend?",
    "header": "backend_selection_strategy_for_this_conversation_including_"
              "every_fallback_we_might_conceivably_need",
    "multiSelect": False,
    "options": [{"label": "a", "description": "d"},
                {"label": "b", "description": "d"}],
  }])
  assert wide_width <= base_width, (
    "question card minimum width grew from %d to %d when the header became "
    "one unbroken token" % (base_width, wide_width))


def exercise_long_option_description_does_not_floor_the_card_width():
  """Option rows must not widen the card either.

  QCheckBox / QRadioButton neither wrap nor elide, so folding the
  model-controlled 'label -- description' text into the button text reports
  its full width as the row's minimumSizeHint. Option descriptions are
  sentence-length by design, so this is ordinary output rather than an edge
  case -- and it re-floors the whole ConversationView the way the inlined Coot
  error did. Same short-payload baseline as the header test, for the same
  font-metrics reason; only the descriptions differ between the two payloads.
  """
  _qapp()
  base_width = _card_min_width([{
    "question": "Which backend?",
    "header": "Backend",
    "multiSelect": False,
    "options": [{"label": "claude_code", "description": "d"},
                {"label": "anthropic", "description": "d"}],
  }])
  wide_width = _card_min_width([{
    "question": "Which backend?",
    "header": "Backend",
    "multiSelect": False,
    "options": [
      {"label": "claude_code",
       "description": "Reuse the existing Claude Code OAuth login so no API "
                      "key is needed; this is the default and works for "
                      "every profile that does not pin its own backend."},
      {"label": "anthropic",
       "description": "Call the Anthropic API directly with "
                      "PHENIX_ANTHROPIC_API_KEY, falling back to "
                      "ANTHROPIC_API_KEY when the prefixed one is unset."},
    ],
  }])
  assert wide_width <= base_width, (
    "question card minimum width grew from %d to %d under sentence-length "
    "option descriptions" % (base_width, wide_width))


def exercise_long_option_description_stays_readable():
  """Unflooring an option row must not cost the user the description text.

  The description is what the user reads to choose between options, so it has
  to stay fully rendered at a narrow width rather than being elided away to an
  ellipsis.
  """
  from qttbx.widgets.chat.question_card import QuestionCard
  app = _qapp()
  desc = ("Reuse the existing Claude Code OAuth login so no API key is "
          "needed; this is the default for every profile.")
  card = QuestionCard()
  card.set_request("r1", [{
    "question": "Which backend?",
    "header": "Backend",
    "multiSelect": False,
    "options": [{"label": "claude_code", "description": desc},
                {"label": "anthropic", "description": "Use an API key."}],
  }])
  card.resize(360, 500)
  card.show()
  app.processEvents()
  rendered = _rendered_text(card)
  assert desc in rendered, (
    "option description is not rendered in full anywhere on the card")
  assert "claude_code" in rendered, "option label went missing"


def _rendered_text(widget):
  """Concatenate the text of every child widget that carries any."""
  parts = []
  for child in widget.findChildren(QtWidgets.QWidget):
    getter = getattr(child, "text", None)
    if not callable(getter):
      continue
    try:
      value = getter()
    except TypeError:
      continue
    if isinstance(value, str):
      parts.append(value)
  return "\n".join(parts)


def exercise():
  exercise_unbroken_question_header_does_not_floor_the_card_width()
  exercise_long_option_description_does_not_floor_the_card_width()
  exercise_long_option_description_stays_readable()
  exercise_single_select_emits_chosen_label()
  exercise_single_select_other_overrides_radio()
  exercise_multi_select_collects_all_checked()
  exercise_multi_select_other_is_appended()
  exercise_multiple_questions_in_one_card()
  exercise_duplicate_question_text_keeps_both_answers()
  exercise_header_uses_assistant_label()
  exercise_submit_disables_buttons_and_hides_card()
  exercise_malformed_question_payload_unblocks_worker()
  exercise_non_dict_options_do_not_crash()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
