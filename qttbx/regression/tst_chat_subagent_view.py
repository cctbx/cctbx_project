"""A finished subagent's transcript is listed and openable as a conversation.

The record was already written under its parent, holding the child's whole
turn-by-turn history -- and nothing surfaced it, so the only way to read a
review was to open the JSON by hand.

It is exposed as a *derived view*, not promoted to a conversation directory:
there stays exactly one copy of the data, deleting the parent takes its reviews
with it, and records written before the feature existed appear with no
migration. The cost of that choice is that the view has no directory, so every
write path must refuse it -- which is what most of this file pins.
"""

import json
import shutil
import tempfile
from pathlib import Path

from libtbx.utils import format_cpu_times, null_out, Sorry
from qttbx.widgets.chat.agent.conversation import (
  Conversation, ContentBlock, Message, SubagentRecord, TokenUsage, now)
from qttbx.widgets.chat.agent.storage import (
  ConversationStorage, split_subagent_view_id, subagent_view_id)
from qttbx.widgets.chat.agent.subagent import subagent_conversation_title


def _storage():
  tmp = tempfile.mkdtemp()
  return tmp, ConversationStorage(Path(tmp), log=null_out())


def _record(sub_id="sa_abc123", final_text="# Adversarial review — m.pdb\n\nbody",
            messages=None, profile_name="phenix_reviewer"):
  return SubagentRecord(
    sub_id=sub_id,
    parent_conversation_id="",
    parent_tool_use_id="tu_1",
    task="You are an adversarial reviewer.",
    profile_name=profile_name,
    model="claude-opus-5",
    started_at=now(),
    finished_at=now(),
    final_text=final_text,
    token_usage=TokenUsage(),
    messages=messages if messages is not None else [
      Message(role="user", timestamp=now(),
              content=[ContentBlock("text", {"text": "brief"})]),
      Message(role="assistant", timestamp=now(),
              content=[ContentBlock("text", {"text": final_text})]),
    ])


def _parent(storage, title="parent"):
  conv = Conversation.new(profile_name="phenix_expert", model="claude-opus-5")
  conv.meta.title = title
  storage.save(conv)
  return conv


def exercise_title_comes_from_the_report_heading():
  """The name is derived, so it needs no model call and no credentials.

  Preamble before the heading is normal -- a child narrates a line or two
  before starting the document -- so the scan cannot simply read line one.
  """
  rec = _record(final_text=(
    "I now have coverage; compiling the report.\n\n"
    "---\n\n"
    "# Adversarial review — model_v2 (2.20 Å)\n\n"
    "## Coverage\n"))
  assert subagent_conversation_title(rec) == \
    "Adversarial review — model_v2 (2.20 Å)", subagent_conversation_title(rec)


def exercise_a_hash_inside_a_fence_is_not_the_title():
  """A shell comment in a fenced block is not a heading.

  Reviews quote commands, and `# phenix.refine ...` at the start of one would
  otherwise name the conversation after a command line.
  """
  rec = _record(final_text=(
    "```\n# phenix.refine model.pdb data.mtz\n```\n\n"
    "# Adversarial review — the real heading\n"))
  assert subagent_conversation_title(rec) == \
    "Adversarial review — the real heading", subagent_conversation_title(rec)


def exercise_a_report_with_no_heading_falls_back_to_naming_the_run():
  """A child that failed, was capped, or answered in prose still needs a name.

  The fallback names the run rather than inventing a description of it: an
  empty or truncated report is exactly the case where a confident-sounding
  title would misrepresent what is inside.

  It must also stay NEUTRAL about what the child was. This layer is shared and
  knows nothing about reviewing; the fallback read "Review (...)" for every
  child, which made a generic one indistinguishable from an adversarial review
  to every consumer of these records -- including phenix's autonomous loop,
  which derives a convergence count from them. The profile name carries that
  distinction, so the literal must not assert it.
  """
  rec = _record(final_text="Denied every tool; nothing measured.")
  title = subagent_conversation_title(rec)
  assert title.startswith("Subagent (phenix_reviewer)"), title
  # And a record with no report at all is still nameable.
  assert subagent_conversation_title(_record(final_text="")).startswith(
    "Subagent (phenix_reviewer)")
  # The neutrality itself, asserted on the record the old literal mislabelled:
  # a child of a NON-reviewer profile must not be called a review here.
  generic = subagent_conversation_title(
    _record(final_text="", profile_name="claude_assistant_subagent"))
  assert "review" not in generic.lower(), generic
  assert "claude_assistant_subagent" in generic, generic


def exercise_a_long_heading_is_elided():
  rec = _record(final_text="# " + "x" * 200)
  title = subagent_conversation_title(rec)
  assert len(title) <= 80, len(title)
  assert title.endswith("…"), title


def exercise_view_ids_round_trip_and_ordinary_ids_are_not_views():
  """`load` decides what an id means from the id alone -- no index lookup.

  A stale or corrupt index is a case the storage layer already handles by
  rebuilding; it must not also make a child transcript unopenable.
  """
  vid = subagent_view_id("parent123", "sa_abc")
  assert split_subagent_view_id(vid) == ("parent123", "sa_abc")
  assert split_subagent_view_id("parent123") is None
  assert split_subagent_view_id("") is None
  # Malformed halves are not views either.
  assert split_subagent_view_id("~sa_abc") is None
  assert split_subagent_view_id("parent~") is None


def exercise_a_stored_subagent_is_listed_as_a_readonly_conversation():
  """The listing is what makes the review reachable at all."""
  tmp, storage = _storage()
  try:
    parent = _parent(storage)
    storage.store_subagent(parent.meta.id, _record())
    metas = storage.list_conversations()
    views = [m for m in metas if m.readonly]
    assert len(views) == 1, [m.id for m in metas]
    view = views[0]
    assert view.title == "Adversarial review — m.pdb", view.title
    assert view.profile_name == "phenix_reviewer"
    # The parent is still listed alongside it, not replaced by it.
    assert any(m.id == parent.meta.id and not m.readonly for m in metas)
  finally:
    shutil.rmtree(tmp)


def exercise_records_written_before_the_feature_need_no_migration():
  """Retroactive by construction: the view is derived at scan time.

  The fixture is the state a real pre-existing project is actually in -- an
  index.json that EXISTS and was written before views were listed. Deleting
  the index instead (the first version of this test) exercised the rebuild
  path and passed while the shipped behavior showed no reviews at all: the
  cached index is preferred over the scan, so a stale one hides every review
  already on disk. The generation stamp is what forces the rebuild.
  """
  tmp, storage = _storage()
  try:
    parent = _parent(storage)
    storage.store_subagent(parent.meta.id, _record())
    # Rewrite index.json exactly as an older build left it: real conversation
    # rows, no views, no generation stamp.
    index_path = Path(tmp) / ".phenix_chat" / "index.json"
    doc = json.loads(index_path.read_text(encoding="utf-8"))
    doc.pop("views_generation", None)
    doc["conversations"] = [d for d in doc["conversations"]
                            if not d.get("readonly")]
    index_path.write_text(json.dumps(doc), encoding="utf-8")
    assert [m.title for m in storage.list_conversations() if m.readonly] == \
      ["Adversarial review — m.pdb"], "a pre-existing index hid the review"
    # And the rebuild stamped the index, so the next listing is served from
    # cache rather than rescanning every launch.
    assert json.loads(index_path.read_text(encoding="utf-8")).get(
      "views_generation"), "the rebuilt index was not stamped"
  finally:
    shutil.rmtree(tmp)


def exercise_opening_a_view_yields_the_childs_messages_and_writes_nothing():
  """Opening a review must not mint a directory under its synthetic id."""
  tmp, storage = _storage()
  try:
    parent = _parent(storage)
    storage.store_subagent(parent.meta.id, _record())
    vid = subagent_view_id(parent.meta.id, "sa_abc123")
    conv = storage.load(vid)
    assert conv.meta.readonly is True
    assert conv.meta.id == vid
    # The parent is read back OUT of the id, not from a stored copy.
    from qttbx.widgets.chat.agent.storage import split_subagent_view_id
    assert split_subagent_view_id(conv.meta.id)[0] == parent.meta.id
    assert [m.role for m in conv.messages] == ["user", "assistant"]
    # The record rides along so the view can show what the child reported.
    assert conv.subagents and conv.subagents[0].sub_id == "sa_abc123"
    conv_root = Path(tmp) / ".phenix_chat" / "conversations"
    assert sorted(p.name for p in conv_root.iterdir()) == [parent.meta.id], \
      "opening a view created a directory"
  finally:
    shutil.rmtree(tmp)


def exercise_a_view_cannot_be_saved():
  """Refused loudly, not skipped silently.

  Saving would fork the child's transcript into a second copy that its
  parent's record no longer matches, and a no-op would let a caller report a
  save that never happened. The window gates on `meta.readonly`; this is the
  layer that makes the gate's absence a failure rather than corruption.
  """
  tmp, storage = _storage()
  try:
    parent = _parent(storage)
    storage.store_subagent(parent.meta.id, _record())
    conv = storage.load(subagent_view_id(parent.meta.id, "sa_abc123"))
    try:
      storage.save(conv)
    except Sorry as e:
      assert "read-only" in str(e), str(e)
    else:
      raise AssertionError("saving a read-only view was allowed")
  finally:
    shutil.rmtree(tmp)


def exercise_a_view_cannot_be_deleted_but_its_parent_takes_it_along():
  """Deleting a review would strip evidence out of the parent's record.

  The parent is the only place that decision belongs -- and deleting it must
  still remove its reviews, or the sidebar would keep listing transcripts
  whose conversation is gone.
  """
  tmp, storage = _storage()
  try:
    parent = _parent(storage)
    storage.store_subagent(parent.meta.id, _record())
    vid = subagent_view_id(parent.meta.id, "sa_abc123")
    try:
      storage.delete(vid)
    except Sorry as e:
      assert "belonging to another conversation" in str(e), str(e)
    else:
      raise AssertionError("deleting a read-only view was allowed")
    storage.delete(parent.meta.id)
    assert storage.list_conversations() == []
  finally:
    shutil.rmtree(tmp)


def exercise_an_unreadable_parent_meta_does_not_hide_its_reviews():
  """The transcript is the artifact that cannot be regenerated.

  A parent's meta.json is cheap to lose and its conversation is recoverable
  from messages.json; a review is neither. Losing the review because the file
  next to it was corrupted would discard the more valuable of the two.
  """
  tmp, storage = _storage()
  try:
    parent = _parent(storage)
    storage.store_subagent(parent.meta.id, _record())
    (Path(tmp) / ".phenix_chat" / "conversations" / parent.meta.id
     / "meta.json").write_text("{ not json", encoding="utf-8")
    (Path(tmp) / ".phenix_chat" / "index.json").unlink()
    metas = storage.list_conversations()
    assert [m.title for m in metas] == ["Adversarial review — m.pdb"], \
      [m.title for m in metas]
  finally:
    shutil.rmtree(tmp)


def exercise_the_readonly_flag_survives_the_index():
  """The index is a cache of metas; a view that loses its flag there would be
  offered to the user as an editable conversation on the next launch."""
  tmp, storage = _storage()
  try:
    parent = _parent(storage)
    storage.store_subagent(parent.meta.id, _record())
    storage.list_conversations()              # writes index.json
    fresh = ConversationStorage(Path(tmp), log=null_out())
    views = [m for m in fresh.list_conversations() if m.readonly]
    assert len(views) == 1, "readonly did not survive the index round-trip"
  finally:
    shutil.rmtree(tmp)


def exercise_subject_digests_are_computed_here_and_round_trip():
  """The identity a later reader checks must be computed, never asserted.

  The agent supplying the paths is the same one whose later claim -- "nothing
  changed since that review" -- the digest exists to test, so a digest it
  wrote into the brief would verify nothing. Hashing in this layer is what
  makes the staleness rule checkable rather than remembered, and it is the
  field the whole re-derivation path depends on.
  """
  import hashlib
  from qttbx.widgets.chat.agent.subagent import (
    digest_subject_paths, normalize_subject_paths)
  tmp, storage = _storage()
  try:
    model = Path(tmp) / "model.pdb"
    model.write_bytes(b"ATOM ...\n")
    expected = hashlib.sha256(b"ATOM ...\n").hexdigest()
    digests = digest_subject_paths([str(model)])
    assert digests[str(model)] == expected, digests
    # Unreadable is recorded as unverifiable -- NOT silently omitted, which a
    # reader would see as "no identity recorded" and NOT as "" which it must
    # not read as unchanged.
    missing = digest_subject_paths([str(Path(tmp) / "gone.pdb")])
    assert list(missing.values()) == [""], missing
    # Survives the record round-trip, or the successor reads nothing.
    parent = _parent(storage)
    rec = _record()
    rec.subject_digests = digests
    storage.store_subagent(parent.meta.id, rec)
    back = storage.load_subagent(parent.meta.id, rec.sub_id)
    assert back.subject_digests == digests, back.subject_digests
    # A malformed argument is dropped, not raised: a bad identity must not
    # block a review that is otherwise fine.
    assert normalize_subject_paths("one/path") == ["one/path"]
    assert normalize_subject_paths({"a": 1}) == []
    assert normalize_subject_paths(None) == []
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_title_comes_from_the_report_heading()
  exercise_a_hash_inside_a_fence_is_not_the_title()
  exercise_a_report_with_no_heading_falls_back_to_naming_the_run()
  exercise_a_long_heading_is_elided()
  exercise_view_ids_round_trip_and_ordinary_ids_are_not_views()
  exercise_a_stored_subagent_is_listed_as_a_readonly_conversation()
  exercise_records_written_before_the_feature_need_no_migration()
  exercise_opening_a_view_yields_the_childs_messages_and_writes_nothing()
  exercise_a_view_cannot_be_saved()
  exercise_a_view_cannot_be_deleted_but_its_parent_takes_it_along()
  exercise_an_unreadable_parent_meta_does_not_hide_its_reviews()
  exercise_the_readonly_flag_survives_the_index()
  exercise_subject_digests_are_computed_here_and_round_trip()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
