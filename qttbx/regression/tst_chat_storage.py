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


def exercise_delete_removes_conversation():
  """delete(conv_id) removes the conversation's on-disk directory so it no
  longer appears in list_conversations, while other conversations survive."""
  tmp, storage = _new_storage()
  try:
    c1 = Conversation.new(profile_name="p", model="m", title="First")
    c2 = Conversation.new(profile_name="p", model="m", title="Second")
    storage.save(c1)
    storage.save(c2)
    assert len(storage.list_conversations()) == 2
    c1_dir = storage.conv_dir(c1.meta.id)
    assert c1_dir.exists()
    storage.delete(c1.meta.id)
    # Its directory is gone...
    assert not c1_dir.exists()
    # ...it is no longer listed, and the other conversation survives.
    listing = storage.list_conversations()
    assert len(listing) == 1, listing
    assert {m.title for m in listing} == {"Second"}
    assert storage.conv_dir(c2.meta.id).exists()
  finally:
    shutil.rmtree(tmp)


def exercise_delete_unknown_raises_sorry():
  """delete() of an unknown conv_id raises Sorry, matching load()'s
  not-found contract rather than leaking a raw error -- and leaves the
  existing conversations untouched."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="x")
    storage.save(conv)
    try:
      storage.delete("no-such-conversation")
    except Sorry:
      pass
    else:
      from libtbx.test_utils import Exception_expected
      raise Exception_expected
    assert len(storage.list_conversations()) == 1
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


def exercise_load_corrupt_meta_json_raises_sorry():
  """A conversation whose meta.json is corrupt (invalid JSON) must surface
  as a Sorry, not a raw json.JSONDecodeError -- load()'s documented
  Sorry-only contract."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="x")
    storage.save(conv)
    meta_path = (Path(tmp) / ".phenix_chat" / "conversations" /
                 conv.meta.id / "meta.json")
    meta_path.write_text("{ this is not valid json")
    try:
      storage.load(conv.meta.id)
    except Sorry:
      pass
    else:
      from libtbx.test_utils import Exception_expected
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_load_missing_created_at_raises_sorry():
  """A meta.json missing the required 'created_at' key must surface as a
  Sorry, not a cryptic KeyError ('created_at')."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="x")
    storage.save(conv)
    meta_path = (Path(tmp) / ".phenix_chat" / "conversations" /
                 conv.meta.id / "meta.json")
    with open(meta_path) as fh:
      doc = json.load(fh)
    del doc["created_at"]
    with open(meta_path, "w") as fh:
      json.dump(doc, fh)
    try:
      storage.load(conv.meta.id)
    except Sorry:
      pass
    else:
      from libtbx.test_utils import Exception_expected
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_load_content_block_not_dict_raises_sorry():
  """A messages.json whose content entry is a bare string (so
  _content_block_from_dict hits AttributeError on .get) must surface as a
  Sorry, not a raw AttributeError -- load()'s documented Sorry-only
  contract."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="x")
    storage.save(conv)
    msgs_path = (Path(tmp) / ".phenix_chat" / "conversations" /
                 conv.meta.id / "messages.json")
    with open(msgs_path, "w") as fh:
      json.dump({"schema_version": "1.0",
                 "messages": [{"role": "user", "timestamp": None,
                               "content": ["text"]}]}, fh)
    try:
      storage.load(conv.meta.id)
    except Sorry:
      pass
    else:
      from libtbx.test_utils import Exception_expected
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_load_usage_wrong_type_raises_sorry():
  """A messages.json whose 'usage' is a bare int (so _token_usage_from_dict
  hits TypeError on the membership test) must surface as a Sorry, not a raw
  TypeError -- load()'s documented Sorry-only contract."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="x")
    storage.save(conv)
    msgs_path = (Path(tmp) / ".phenix_chat" / "conversations" /
                 conv.meta.id / "messages.json")
    with open(msgs_path, "w") as fh:
      json.dump({"schema_version": "1.0",
                 "messages": [{"role": "assistant", "timestamp": None,
                               "content": [], "usage": 5}]}, fh)
    try:
      storage.load(conv.meta.id)
    except Sorry:
      pass
    else:
      from libtbx.test_utils import Exception_expected
      raise Exception_expected
  finally:
    shutil.rmtree(tmp)


def exercise_save_failure_cleans_up_tmp():
  """A save that fails mid-write (an un-serializable value) must not leave
  an orphaned .tmp file behind: the atomic write cleans up its temp file."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="x")
    # An un-serializable value in a message makes messages.json's write
    # raise after meta.json has already been committed.
    conv.append(Message(role="user",
                        content=[ContentBlock(type="text",
                                              data={"text": object()})],
                        timestamp=now()))
    try:
      storage.save(conv)
    except TypeError:
      pass
    else:
      from libtbx.test_utils import Exception_expected
      raise Exception_expected
    conv_dir = Path(tmp) / ".phenix_chat" / "conversations" / conv.meta.id
    orphans = sorted(str(p) for p in conv_dir.rglob("*.tmp"))
    assert orphans == [], orphans
  finally:
    shutil.rmtree(tmp)


def exercise_store_attachment_cleans_up_tmp_on_replace_failure():
  """A store_attachment whose os.replace fails mid-write must not leave an
  orphaned <file>.tmp behind: the content-addressed write uses the same
  try/finally cleanup _atomic_write_json has. The originating OSError still
  propagates; the contract under test is solely the temp cleanup.

  Revert-proof: an inline tmp+os.replace without the cleanup leaves the .tmp
  on disk and the orphan assert fails."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    att_dir = (Path(tmp) / ".phenix_chat" / "conversations" /
               conv.meta.id / "attachments")
    orig_replace = os.replace
    def boom(src, dst):
      raise OSError("simulated os.replace failure")
    os.replace = boom
    raised = False
    try:
      storage.store_attachment(conv.meta.id, b"PNG bytes here", mime="image/png")
    except OSError:
      raised = True
    finally:
      os.replace = orig_replace
    assert raised, "the OSError from os.replace must still propagate out"
    orphans = (sorted(str(p) for p in att_dir.glob("*.tmp"))
               if att_dir.exists() else [])
    assert orphans == [], orphans
  finally:
    shutil.rmtree(tmp)


def exercise_all_token_usage_fields_round_trip():
  """Every TokenUsage field on an assistant message survives save()+load().
  Dataclass equality compares all fields, so a usage field added later is
  covered automatically -- it must not silently drop on persist."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="usage")
    usage = TokenUsage(input=11, output=22, cache_read=33, cache_creation=44)
    conv.append(Message(role="assistant",
                        content=[ContentBlock(type="text",
                                              data={"text": "hi"})],
                        timestamp=now(), stop_reason="end_turn", usage=usage))
    storage.save(conv)
    loaded = storage.load(conv.meta.id)
    assert loaded.messages[-1].usage == usage, loaded.messages[-1].usage
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


def exercise_delete_rejects_unsafe_conv_id():
  """delete() must reject a conv_id that is not a single safe path segment
  (traversal / nested / absolute), raising Sorry rather than touching
  anything outside the per-conversation directory."""
  tmp, storage = _new_storage()
  try:
    for bad in ("..", "../x", "a/b", "a\\b", "/etc/passwd", ".", ""):
      try:
        storage.delete(bad)
      except Sorry:
        pass
      else:
        raise RuntimeError("delete(%r) should have raised Sorry" % (bad,))
  finally:
    shutil.rmtree(tmp)


def exercise_delete_refuses_symlink_escape():
  """A conversations/<id> symlink pointing outside the chat root is refused
  by the resolve()/parents guard, and the symlink target is left intact."""
  tmp, storage = _new_storage()
  outside = tempfile.mkdtemp()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)                       # realize the conversations/ tree
    conv_root = storage.conv_dir(conv.meta.id).parent
    victim = os.path.join(outside, "precious.txt")
    with open(victim, "w") as fh:
      fh.write("keep")
    link = os.path.join(str(conv_root), "evil")
    try:
      os.symlink(outside, link)
    except (OSError, NotImplementedError):
      print("symlink not supported; skipping escape check")
      return
    try:
      storage.delete("evil")
    except Sorry:
      pass
    else:
      raise RuntimeError("delete('evil') should have refused the escape")
    assert os.path.exists(victim), "symlink target must be untouched"
  finally:
    shutil.rmtree(tmp)
    shutil.rmtree(outside, ignore_errors=True)


def exercise():
  exercise_lazy_directory_creation()
  exercise_save_then_load_roundtrip()
  exercise_meta_backend_and_per_turn_stamp_roundtrip()
  exercise_atomic_write_interruption_leaves_prior_intact()
  exercise_attachment_dedup_by_sha256()
  exercise_index_listing()
  exercise_index_rebuild_when_missing()
  exercise_delete_removes_conversation()
  exercise_delete_unknown_raises_sorry()
  exercise_delete_rejects_unsafe_conv_id()
  exercise_delete_refuses_symlink_escape()
  exercise_subagent_store_and_load()
  exercise_schema_version_check_rejects_future_version()
  exercise_schema_version_default_accepts_current()
  exercise_load_corrupt_meta_json_raises_sorry()
  exercise_load_missing_created_at_raises_sorry()
  exercise_load_content_block_not_dict_raises_sorry()
  exercise_load_usage_wrong_type_raises_sorry()
  exercise_save_failure_cleans_up_tmp()
  exercise_store_attachment_cleans_up_tmp_on_replace_failure()
  exercise_all_token_usage_fields_round_trip()
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
