"""Storage layer test: ConversationMeta.agent_session_id round-trips
through _meta_to_dict / _meta_from_dict, and a pre-feature meta.json
(missing the key) loads cleanly as None rather than raising KeyError."""

import json

from libtbx.utils import format_cpu_times
from qttbx.widgets.chat.agent.storage import (
  _json_default, _meta_to_dict, _meta_from_dict)
from qttbx.widgets.chat.agent.conversation import ConversationMeta, now


def _minimal_meta(conv_id):
  """Build a ConversationMeta with only its required fields populated."""
  ts = now()
  return ConversationMeta(id=conv_id, title="t", profile_name="p",
                          model="m", created_at=ts, updated_at=ts)


def _round_trip(doc):
  """Serialize to JSON and back, as ConversationStorage does on disk.

  This converts the datetime values in ``doc`` to ISO strings (the form
  ``_meta_from_dict`` expects), so the test exercises the real persistence
  path rather than tripping over raw datetimes.
  """
  return json.loads(json.dumps(doc, default=_json_default))


def exercise_agent_session_id_round_trips():
  """agent_session_id survives _meta_to_dict -> _meta_from_dict."""
  m = _minimal_meta("c1")
  m.agent_session_id = "sess-123"
  back = _meta_from_dict(_round_trip(_meta_to_dict(m)))
  assert back.agent_session_id == "sess-123", back.agent_session_id


def exercise_missing_agent_session_id_loads_as_none():
  """Old meta.json without the key loads cleanly (None), not KeyError."""
  doc = _round_trip(_meta_to_dict(_minimal_meta("c2")))
  doc.pop("agent_session_id", None)            # simulate pre-feature meta.json
  back = _meta_from_dict(doc)
  assert back.agent_session_id is None, back.agent_session_id


def exercise_public_conv_dir_matches_private():
  """conv_dir() is the public accessor for a conversation's directory and
  resolves under <root>/conversations/."""
  import tempfile
  from qttbx.widgets.chat.agent.storage import ConversationStorage
  from libtbx.utils import null_out
  st = ConversationStorage(project_dir=tempfile.mkdtemp(), log=null_out())
  d = st.conv_dir("abc")
  assert d == st._conv_dir("abc")
  # Not a tautology: the path resolves under <root>/conversations/.
  assert d.parent == st.root / "conversations", d
  assert d.parent.parent == st.root, d


def exercise():
  exercise_agent_session_id_round_trips()
  exercise_missing_agent_session_id_loads_as_none()
  exercise_public_conv_dir_matches_private()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
