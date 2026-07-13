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

  def __init__(self, parent=None, assistant_label="Assistant"):
    super().__init__(parent)
    self.setFrameShape(QtWidgets.QFrame.StyledPanel)
    # Active assistant display name, used in the header ("<name> needs an
    # answer:"); supplied by ConversationView from the current backend.
    self._assistant_label = assistant_label or "Assistant"
    self._request_id = ""
    self._questions = []
    # True once an answer has been emitted (via Submit or the malformed-
    # payload auto-answer), so the two paths never double-emit.
    self._resolved = False
    # Per-question state: list of dicts with:
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
    self._resolved = False
    self._questions = self._normalize_questions(questions)
    self._rebuild()
    if not self._questions:
      # A malformed/empty payload (not a list of question dicts) leaves
      # nothing to render. The agent worker that emitted
      # phenix_ask_user_question is parked on question_queue.get() waiting
      # for an answer a Submit click can no longer produce -- so release it
      # ourselves with a safe error answer instead of hanging the turn.
      # Deferred via singleShot because the caller
      # (ConversationView.add_question_request) connects `answered` only
      # AFTER set_request returns; a synchronous emit would be dropped.
      self.hide()
      QtCore.QTimer.singleShot(0, self._emit_unanswerable)

  @staticmethod
  def _normalize_questions(questions):
    """Coerce a raw questions payload into a list of question dicts.

    The agent is supposed to send a list of ``{question, options, ...}``
    dicts, but a malformed tool call can send a bare dict, a string,
    ``None``, or a list carrying non-dict entries -- any of which would
    otherwise crash ``_build_one_question`` on ``q.get(...)`` and strand
    the parked worker. A bare dict becomes a single-question list; anything
    that is not a list/tuple yields no questions; non-dict entries are
    dropped.

    Parameters
    ----------
    questions : object
        The raw payload from the tool call.

    Returns
    -------
    list of dict
        The questions safe to render (possibly empty).
    """
    if isinstance(questions, dict):
      questions = [questions]
    elif not isinstance(questions, (list, tuple)):
      return []
    return [q for q in questions if isinstance(q, dict)]

  def _emit_unanswerable(self):
    """Emit a safe error answer for a payload that rendered no questions.

    Releases the worker parked in ``AgentSession._await_question_answer``
    (see :meth:`set_request`). Guarded by ``_resolved`` so it never
    double-emits alongside a real Submit.
    """
    if self._resolved:
      return
    self._resolved = True
    self.answered.emit(
      self._request_id, {"error": "no valid questions to display"})

  @staticmethod
  def _unique_answer_key(answers, text):
    """Return a key for ``text`` that does not collide in ``answers``.

    Two questions can carry the same ``question`` string (or both omit it,
    yielding ""), and the model reads the answers keyed by question text (the
    tool contract). Keeping the text as the key but appending an ordinal
    suffix on a collision -- "Same?", "Same? (2)", ... -- keeps that contract
    while ensuring each of the N questions yields its own entry; keying by raw
    text would let the second write clobber the first and silently drop one of
    the user's selections.

    Parameters
    ----------
    answers : dict
        Keys already assigned in this submit.
    text : str
        The question's text (its preferred key).

    Returns
    -------
    str
        ``text`` if free, else ``text`` with the lowest free ordinal suffix.
    """
    if text not in answers:
      return text
    n = 2
    while ("%s (%d)" % (text, n)) in answers:
      n += 1
    return "%s (%d)" % (text, n)

  def click_submit(self):
    """Collect answers from every question and emit them.

    Programmatic submit (tests use this; the Submit button calls the
    same path). Emits ``answered(request_id, answers)``.
    """
    if self._resolved:
      return
    self._resolved = True
    answers = {}
    for q, state in zip(self._questions, self._question_state):
      text = q.get("question", "")
      if state["multi_select"]:
        picked = [label for btn, label in state["option_buttons"]
                  if btn.isChecked()]
        other = state["other_edit"].text().strip()
        if other:
          picked.append(other)
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
      answers[self._unique_answer_key(answers, text)] = picked
    if self._buttons_widget is not None:
      self._buttons_widget.setEnabled(False)
    self.hide()
    self.answered.emit(self._request_id, answers)

  def is_resolved(self):
    """Report whether an answer (or the auto-error) has been emitted.

    Once resolved, the card is hidden and its Submit disabled.
    ``ConversationView`` consults this so its turn-end sweep never
    re-finalizes an already-answered card.

    Returns
    -------
    bool
        ``True`` once an answer has been emitted or the card finalized.
    """
    return self._resolved

  def finalize(self):
    """Disable this card WITHOUT emitting an answer.

    Called when the card's turn ends while still unanswered (the user stopped
    the turn) so a later Submit can't emit a stale answer whose ``request_id``
    the parked worker no longer waits on -- mirroring ``ToolApprovalCard.
    finalize`` and the approval-misroute guard. Marks the card resolved (so
    ``click_submit`` won't fire), disables its buttons, and hides it. A no-op
    once a real answer has been emitted.
    """
    if self._resolved:
      return
    self._resolved = True
    if self._buttons_widget is not None:
      self._buttons_widget.setEnabled(False)
    self.hide()

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
    head = QtWidgets.QLabel("%s needs an answer:" % self._assistant_label, box)
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
    if not isinstance(options, (list, tuple)):
      options = []
    option_buttons = []
    group = None if multi_select else QtWidgets.QButtonGroup(frame)
    if group is not None:
      group.setExclusive(True)
    for opt in options:
      # Tolerate a malformed option: a {label, description} dict is the
      # spec, but a bare string is coerced to its label and anything else
      # is skipped -- never crash on opt.get(...).
      if isinstance(opt, dict):
        label = opt.get("label", "")
        desc = opt.get("description", "")
      elif isinstance(opt, str):
        label, desc = opt, ""
      else:
        continue
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
    return box
