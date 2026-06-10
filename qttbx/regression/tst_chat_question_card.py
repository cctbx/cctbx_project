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


def exercise():
  exercise_single_select_emits_chosen_label()
  exercise_single_select_other_overrides_radio()
  exercise_multi_select_collects_all_checked()
  exercise_multi_select_other_is_appended()
  exercise_multiple_questions_in_one_card()
  exercise_header_uses_assistant_label()
  exercise_submit_disables_buttons_and_hides_card()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
