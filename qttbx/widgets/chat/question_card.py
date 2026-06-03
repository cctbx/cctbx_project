"""Inline card asking the user one or more multiple-choice questions.

Triggered when the agent calls the ``phenix_ask_user_question`` tool.
The card renders each question with its options as radio buttons (or
checkboxes when ``multiSelect`` is True), an optional ``Other`` textbox
for free-form answers, and a Submit button. On submit it emits
``answered(request_id, answers)`` -- a dict keyed by question text with
the selected label (or list of labels).
"""

from qttbx.qt import QtCore, QtWidgets


class QuestionCard(QtWidgets.QFrame):
  """Single card carrying one or more questions.

  One Submit click emits all answers at once. Pattern mirrors
  ``ToolApprovalCard`` so the chat infrastructure can render it the
  same way.
  """

  answered = QtCore.Signal(str, dict)            # request_id, answers

  def __init__(self, parent=None):
    super().__init__(parent)
    self.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self._request_id = ""
    self._questions = []
    # Per-question state: list of dicts with:
    #   'group' (QButtonGroup or None for multiSelect),
    #   'option_buttons' (list of (QAbstractButton, label_str)),
    #   'other_edit' (QLineEdit), 'multi_select' (bool).
    self._question_state = []
    self._layout = QtWidgets.QVBoxLayout(self)
    self._layout.setContentsMargins(8, 6, 8, 6)
    self._content_widget = None
    self._buttons_widget = None

  # ---- public API ----------------------------------------------------------

  def set_request(self, request_id, questions):
    self._request_id = request_id
    self._questions = list(questions or [])
    self._rebuild()

  def click_submit(self):
    """Collect answers from every question and emit them.

    Programmatic submit (tests use this; the Submit button calls the
    same path). Emits ``answered(request_id, answers)``.
    """
    answers = {}
    for q, state in zip(self._questions, self._question_state):
      text = q.get("question", "")
      if state["multi_select"]:
        picked = [label for btn, label in state["option_buttons"]
                  if btn.isChecked()]
        other = state["other_edit"].text().strip()
        if other:
          picked.append(other)
        answers[text] = picked
      else:
        picked = None
        for btn, label in state["option_buttons"]:
          if btn.isChecked():
            picked = label
            break
        other = state["other_edit"].text().strip()
        if other:
          # Free-form text wins over a radio selection when both are set;
          # the user's intent is whatever they last typed.
          picked = other
        answers[text] = picked
    if self._buttons_widget is not None:
      self._buttons_widget.setEnabled(False)
    self.hide()
    self.answered.emit(self._request_id, answers)

  # ---- UI build ------------------------------------------------------------

  def _rebuild(self):
    if self._content_widget is not None:
      self._content_widget.setParent(None)
    if self._buttons_widget is not None:
      self._buttons_widget.setParent(None)
    self._question_state = []
    self._content_widget = self._build_content()
    self._buttons_widget = self._build_buttons()
    self._layout.addWidget(self._content_widget)
    self._layout.addWidget(self._buttons_widget)

  def _build_content(self):
    box = QtWidgets.QWidget(self)
    layout = QtWidgets.QVBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
    head = QtWidgets.QLabel("Claude needs an answer:", box)
    f = head.font(); f.setBold(True); head.setFont(f)
    layout.addWidget(head)
    for q in self._questions:
      layout.addWidget(self._build_one_question(q, box))
    return box

  def _build_one_question(self, q, parent):
    """Build one question block.

    Records per-question state in ``self._question_state`` so
    ``click_submit`` can read the selections back.

    Parameters
    ----------
    q : dict
        Question spec with ``header``, ``question``, ``multiSelect``,
        and ``options`` keys.
    parent : QtWidgets.QWidget
        Parent widget for the constructed block.

    Returns
    -------
    QtWidgets.QFrame
        The frame holding the rendered question.
    """
    frame = QtWidgets.QFrame(parent)
    frame.setFrameShape(QtWidgets.QFrame.NoFrame)
    v = QtWidgets.QVBoxLayout(frame)
    v.setContentsMargins(0, 4, 0, 4)
    header = q.get("header")
    if header:
      hl = QtWidgets.QLabel("[%s] " % header, frame)
      hl.setStyleSheet("color: palette(mid);")
      v.addWidget(hl)
    text = q.get("question", "")
    if text:
      tl = QtWidgets.QLabel(text, frame)
      tl.setWordWrap(True)
      v.addWidget(tl)
    multi_select = bool(q.get("multiSelect", False))
    options = q.get("options", []) or []
    option_buttons = []
    group = None if multi_select else QtWidgets.QButtonGroup(frame)
    if group is not None:
      group.setExclusive(True)
    for opt in options:
      label = opt.get("label", "")
      desc = opt.get("description", "")
      text_full = "%s -- %s" % (label, desc) if desc else label
      if multi_select:
        btn = QtWidgets.QCheckBox(text_full, frame)
      else:
        btn = QtWidgets.QRadioButton(text_full, frame)
        group.addButton(btn)
      v.addWidget(btn)
      option_buttons.append((btn, label))
    # "Other" free-form input -- always shown so the user can type a
    # different answer if none of the canned options fit.
    other_row = QtWidgets.QHBoxLayout()
    other_row.setContentsMargins(0, 4, 0, 0)
    other_lbl = QtWidgets.QLabel("Other:", frame)
    other_lbl.setStyleSheet("color: palette(mid);")
    other_edit = QtWidgets.QLineEdit(frame)
    other_edit.setPlaceholderText("type a custom answer (optional)")
    other_row.addWidget(other_lbl)
    other_row.addWidget(other_edit, 1)
    v.addLayout(other_row)
    self._question_state.append({
      "group": group,
      "option_buttons": option_buttons,
      "other_edit": other_edit,
      "multi_select": multi_select,
    })
    return frame

  def _build_buttons(self):
    box = QtWidgets.QWidget(self)
    layout = QtWidgets.QHBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addStretch(1)
    submit = QtWidgets.QPushButton("Submit", box)
    submit.clicked.connect(self.click_submit)
    layout.addWidget(submit)
    self._submit_btn = submit
    return box
