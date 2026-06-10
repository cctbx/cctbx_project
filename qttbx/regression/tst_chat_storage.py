"""Storage layer tests: ConversationStorage (atomic JSON, dedup,
subagents, schema version) plus the project-directory / chat-root
path resolvers it depends on."""

import json
import os
import shutil
import tempfile
from pathlib import Path

from libtbx.utils import format_cpu_times, null_out, Sorry
from qttbx.widgets.chat.agent.conversation import (
  Conversation, ContentBlock, Message, TokenUsage, now)
from qttbx.widgets.chat.agent.paths import (
  chat_root_for, resolve_project_dir)
from qttbx.widgets.chat.agent.storage import ConversationStorage


def _new_storage():
  tmp = tempfile.mkdtemp()
  return tmp, ConversationStorage(Path(tmp), log=null_out())


def exercise_lazy_directory_creation():
  """Section 6.3: .phenix_chat/ is NOT created on construction; only on
  first persistence write."""
  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(Path(tmp), log=null_out())
    assert not (Path(tmp) / ".phenix_chat").exists()
    # Listing on an empty/missing root returns empty list, doesn't create.
    assert storage.list_conversations() == []
    assert not (Path(tmp) / ".phenix_chat").exists()
  finally:
    shutil.rmtree(tmp)


def exercise_save_then_load_roundtrip():
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="phenix_expert",
                            model="claude-opus-4-7",
                            title="Test")
    conv.append(Message(role="user",
                        content=[ContentBlock(type="text",
                                              data={"text": "hi"})],
                        timestamp=now()))
    storage.save(conv)
    # Now .phenix_chat/ should exist
    assert (Path(tmp) / ".phenix_chat" / "conversations" / conv.meta.id /
            "messages.json").exists()
    loaded = storage.load(conv.meta.id)
    assert loaded.meta.id == conv.meta.id
    assert loaded.meta.title == "Test"
    assert len(loaded.messages) == 1
    assert loaded.messages[0].content[0].data["text"] == "hi"
  finally:
    shutil.rmtree(tmp)


def exercise_meta_backend_and_per_turn_stamp_roundtrip():
  """The conversation meta carries a backend; each assistant message
  carries the model + backend that produced it. Both survive a
  save/load round-trip (older files without the fields load with empty
  defaults -- covered by the other round-trip tests)."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="phenix_expert",
                            model="claude-opus-4-8",
                            backend="anthropic",
                            title="Backends")
    conv.append(Message(role="user",
                        content=[ContentBlock(type="text",
                                              data={"text": "hi"})],
                        timestamp=now()))
    conv.append(Message(role="assistant",
                        content=[ContentBlock(type="text",
                                              data={"text": "hello"})],
                        timestamp=now(), stop_reason="end_turn",
                        model="claude-opus-4-8", backend="anthropic"))
    storage.save(conv)
    loaded = storage.load(conv.meta.id)
    assert loaded.meta.backend == "anthropic"
    last = loaded.messages[-1]
    assert last.role == "assistant"
    assert last.model == "claude-opus-4-8", last.model
    assert last.backend == "anthropic", last.backend
    # The user message has no stamp; the fields default to None.
    assert loaded.messages[0].model is None
    assert loaded.messages[0].backend is None
  finally:
    shutil.rmtree(tmp)


def exercise_atomic_write_interruption_leaves_prior_intact():
  """If a save is interrupted mid-write, the prior version stays valid.
  We simulate by writing once, then writing a tmp file that we leave behind
  (as if a crash). The .tmp file should not affect load."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="v1")
    storage.save(conv)
    msgs_path = (Path(tmp) / ".phenix_chat" / "conversations" /
                 conv.meta.id / "messages.json")
    # Leave a corrupt .tmp file behind.
    (msgs_path.with_suffix(".json.tmp")).write_text("not valid json")
    loaded = storage.load(conv.meta.id)
    assert loaded.meta.title == "v1"
  finally:
    shutil.rmtree(tmp)


def exercise_attachment_dedup_by_sha256():
  """Section 6.2 / Section 17.7. Same bytes → same sha → one file on disk,
  one Attachment record returned twice."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    data = b"PNG bytes here"
    a1 = storage.store_attachment(conv.meta.id, data, mime="image/png")
    a2 = storage.store_attachment(conv.meta.id, data, mime="image/png")
    assert a1.sha256 == a2.sha256
    # Only one file on disk
    att_dir = (Path(tmp) / ".phenix_chat" / "conversations" /
               conv.meta.id / "attachments")
    files = sorted(p.name for p in att_dir.iterdir())
    assert len(files) == 1
    # load_attachment returns the bytes
    assert storage.load_attachment(conv.meta.id, a1.sha256) == data
  finally:
    shutil.rmtree(tmp)


def exercise_index_listing():
  tmp, storage = _new_storage()
  try:
    c1 = Conversation.new(profile_name="p", model="m", title="First")
    c2 = Conversation.new(profile_name="p", model="m", title="Second")
    storage.save(c1)
    storage.save(c2)
    listing = storage.list_conversations()
    assert len(listing) == 2
    titles = {m.title for m in listing}
    assert titles == {"First", "Second"}
  finally:
    shutil.rmtree(tmp)


def exercise_index_rebuild_when_missing():
  """index.json is a cache — should rebuild from on-disk truth if absent."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="Recovered")
    storage.save(conv)
    index_path = Path(tmp) / ".phenix_chat" / "index.json"
    assert index_path.exists()
    index_path.unlink()
    # Fresh storage instance; should rebuild on first list_conversations
    storage2 = ConversationStorage(Path(tmp), log=null_out())
    listing = storage2.list_conversations()
    assert len(listing) == 1
    assert listing[0].title == "Recovered"
    # And the index file got recreated
    assert index_path.exists()
  finally:
    shutil.rmtree(tmp)


def exercise_subagent_store_and_load():
  from qttbx.widgets.chat.agent.conversation import SubagentRecord
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    rec = SubagentRecord(
      sub_id="sa_test",
      parent_conversation_id=conv.meta.id,
      parent_tool_use_id="toolu_1",
      task="monitor",
      profile_name="phenix_expert_subagent",
      model="claude-opus-4-7",
      started_at=now(),
      finished_at=now(),
      final_text="done",
      token_usage=TokenUsage(input=10, output=5),
      messages=[])
    storage.store_subagent(conv.meta.id, rec)
    loaded = storage.load_subagent(conv.meta.id, "sa_test")
    assert loaded.sub_id == "sa_test"
    assert loaded.final_text == "done"
    assert loaded.token_usage.input == 10
  finally:
    shutil.rmtree(tmp)


def exercise_schema_version_check_rejects_future_version():
  """load() raises Sorry if a JSON file has an unknown schema_version
  (Section 6.6 migration seam)."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="x")
    storage.save(conv)
    # Corrupt the messages.json's schema_version
    msgs_path = (Path(tmp) / ".phenix_chat" / "conversations" /
                 conv.meta.id / "messages.json")
    with open(msgs_path) as fh:
      doc = json.load(fh)
    doc["schema_version"] = "99.0"  # future version
    with open(msgs_path, "w") as fh:
      json.dump(doc, fh)
    try:
      storage.load(conv.meta.id)
    except Sorry as e:
      assert "99.0" in str(e)
      assert "schema_version" in str(e)
    else:
      from libtbx.test_utils import Exception_expected
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_schema_version_default_accepts_current():
  """A doc with the current schema_version loads cleanly."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="ok")
    storage.save(conv)
    # No corruption; load should succeed.
    loaded = storage.load(conv.meta.id)
    assert loaded.meta.title == "ok"
  finally:
    shutil.rmtree(tmp)


# ---- paths: project directory + chat-data root resolution ----------------


def exercise_resolve_project_dir_cli_arg_wins():
  tmp = tempfile.mkdtemp()
  try:
    result = resolve_project_dir(cli_arg=Path(tmp))
    assert result == Path(tmp).resolve()
  finally:
    os.rmdir(tmp)


def exercise_resolve_project_dir_embedded_arg():
  tmp = tempfile.mkdtemp()
  try:
    result = resolve_project_dir(embedded_arg=Path(tmp))
    assert result == Path(tmp).resolve()
  finally:
    os.rmdir(tmp)


def exercise_resolve_project_dir_env_var():
  tmp = tempfile.mkdtemp()
  saved = os.environ.get("PHENIX_PROJECT_DIR")
  os.environ["PHENIX_PROJECT_DIR"] = tmp
  try:
    result = resolve_project_dir()
    assert result == Path(tmp).resolve()
  finally:
    if saved is None:
      del os.environ["PHENIX_PROJECT_DIR"]
    else:
      os.environ["PHENIX_PROJECT_DIR"] = saved
    os.rmdir(tmp)


def exercise_resolve_project_dir_cli_beats_env():
  tmp_cli = tempfile.mkdtemp()
  tmp_env = tempfile.mkdtemp()
  saved = os.environ.get("PHENIX_PROJECT_DIR")
  os.environ["PHENIX_PROJECT_DIR"] = tmp_env
  try:
    result = resolve_project_dir(cli_arg=Path(tmp_cli))
    assert result == Path(tmp_cli).resolve()
  finally:
    if saved is None:
      del os.environ["PHENIX_PROJECT_DIR"]
    else:
      os.environ["PHENIX_PROJECT_DIR"] = saved
    os.rmdir(tmp_cli)
    os.rmdir(tmp_env)


def exercise_resolve_project_dir_falls_back_to_cwd():
  tmp = tempfile.mkdtemp()
  saved_cwd = os.getcwd()
  saved_env = os.environ.get("PHENIX_PROJECT_DIR")
  if saved_env is not None:
    del os.environ["PHENIX_PROJECT_DIR"]
  try:
    os.chdir(tmp)
    result = resolve_project_dir()
    assert result == Path(tmp).resolve()
  finally:
    os.chdir(saved_cwd)
    if saved_env is not None:
      os.environ["PHENIX_PROJECT_DIR"] = saved_env
    os.rmdir(tmp)


def exercise_chat_root_for_default():
  tmp = tempfile.mkdtemp()
  saved_env = os.environ.get("PHENIX_CHAT_HOME")
  if saved_env is not None:
    del os.environ["PHENIX_CHAT_HOME"]
  try:
    root = chat_root_for(Path(tmp))
    assert root == Path(tmp) / ".phenix_chat"
    # Path is not created -- that's lazy.
    assert not root.exists()
  finally:
    if saved_env is not None:
      os.environ["PHENIX_CHAT_HOME"] = saved_env
    os.rmdir(tmp)


def exercise_chat_root_for_env_override():
  tmp_project = tempfile.mkdtemp()
  tmp_override = tempfile.mkdtemp()
  saved = os.environ.get("PHENIX_CHAT_HOME")
  os.environ["PHENIX_CHAT_HOME"] = tmp_override
  try:
    root = chat_root_for(Path(tmp_project))
    assert root == Path(tmp_override)
  finally:
    if saved is None:
      del os.environ["PHENIX_CHAT_HOME"]
    else:
      os.environ["PHENIX_CHAT_HOME"] = saved
    os.rmdir(tmp_project)
    os.rmdir(tmp_override)


def exercise():
  exercise_lazy_directory_creation()
  exercise_save_then_load_roundtrip()
  exercise_meta_backend_and_per_turn_stamp_roundtrip()
  exercise_atomic_write_interruption_leaves_prior_intact()
  exercise_attachment_dedup_by_sha256()
  exercise_index_listing()
  exercise_index_rebuild_when_missing()
  exercise_subagent_store_and_load()
  exercise_schema_version_check_rejects_future_version()
  exercise_schema_version_default_accepts_current()
  exercise_resolve_project_dir_cli_arg_wins()
  exercise_resolve_project_dir_embedded_arg()
  exercise_resolve_project_dir_env_var()
  exercise_resolve_project_dir_cli_beats_env()
  exercise_resolve_project_dir_falls_back_to_cwd()
  exercise_chat_root_for_default()
  exercise_chat_root_for_env_override()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
