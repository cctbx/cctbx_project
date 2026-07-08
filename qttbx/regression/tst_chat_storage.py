"""Storage layer tests: ConversationStorage (atomic JSON, dedup,
subagents, schema version) plus the project-directory / chat-root
path resolvers it depends on."""

import json
import os
import shutil
import tempfile
import threading
import time
from pathlib import Path

from libtbx.utils import format_cpu_times, null_out, Sorry
from qttbx.widgets.chat.agent.conversation import (
  Conversation, ContentBlock, Message, TokenUsage, now)
from qttbx.widgets.chat.agent.paths import (
  chat_root_for, resolve_project_dir)
from qttbx.widgets.chat.agent.storage import (
  ConversationStorage, LOCK_ACQUIRED, LOCK_HELD, LOCK_UNAVAILABLE)


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


def exercise_atomic_replace_fails_fast_on_posix():
  """On POSIX os.replace never blocks on an open target, so a PermissionError
  there is a PERMANENT condition (a read-only tree / EACCES) that will not clear:
  _atomic_replace must re-raise it IMMEDIATELY, not spin the Windows
  sharing-violation retry for _REPLACE_RETRY_SECONDS (~2s). Guards the
  docstring's 'no-op fast path on POSIX' claim, which an unconditional retry loop
  silently broke -- a 2s GUI-thread freeze on every EPERM turn-end save.

  Revert-proof: an unconditional retry loop calls os.replace many times over ~2s
  before re-raising, so both the call-count and the elapsed assertion fail."""
  if os.name != "posix":
    return                                   # the retry is a Windows-only concession
  from qttbx.widgets.chat.agent.storage import _atomic_replace
  calls = []
  orig_replace = os.replace
  def boom(src, dst):
    calls.append((src, dst))
    raise PermissionError("simulated EPERM")
  os.replace = boom
  raised = False
  start = time.monotonic()
  try:
    _atomic_replace("src.tmp", "dst")
  except PermissionError:
    raised = True
  finally:
    os.replace = orig_replace
  elapsed = time.monotonic() - start
  assert raised, "a POSIX PermissionError must propagate"
  assert calls == [("src.tmp", "dst")], \
    "POSIX must try os.replace exactly once (no Windows retry spin): %r" % calls
  assert elapsed < 0.5, \
    "POSIX must fail fast, not spin the ~2s retry (elapsed=%.3fs)" % elapsed


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


def exercise_read_json_retries_windows_sharing_violation():
  """On Windows a reader's ``open`` can fail with ``PermissionError``
  (ERROR_ACCESS_DENIED) during the brief window a concurrent saver holds the file
  for ``os.replace`` (the atomic-write final step) -- a transient sharing
  violation, not corruption. ``_read_json`` retries it, mirroring
  ``_atomic_replace``'s write-side retry, so ``load`` doesn't raise a false
  'corrupt' Sorry under concurrent saves (the Windows CI failure in
  ``exercise_concurrent_saves_do_not_corrupt``). The Windows-only retry branch is
  forced here so it runs cross-platform.

  Revert-proof: without the retry the first PermissionError propagates and the
  read never returns the document."""
  import builtins
  from qttbx.widgets.chat.agent import storage as _store
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m", title="win")
    storage.save(conv)
    meta_path = str(storage.conv_dir(conv.meta.id) / "meta.json")
    opens = [0]
    real_open = builtins.open
    def flaky_open(file, *a, **k):
      if str(file) == meta_path:
        opens[0] += 1
        if opens[0] <= 2:              # first two opens hit the sharing violation
          raise PermissionError(13, "simulated sharing violation")
      return real_open(file, *a, **k)
    orig_name = os.name
    try:
      os.name = "nt"                   # force the Windows retry branch on any host
      builtins.open = flaky_open
      doc = _store._read_json(meta_path)
    finally:
      builtins.open = real_open
      os.name = orig_name
    assert opens[0] == 3, "expected 2 retries then a success, got %d opens" % opens[0]
    assert doc["id"] == conv.meta.id, doc
  finally:
    shutil.rmtree(tmp)


def exercise_concurrent_saves_do_not_corrupt():
  """Two threads saving the same conversation, while two others read it, must
  never observe or produce a half-written meta.json / messages.json /
  index.json. Guards the unique-temp fix; fails with a shared fixed .tmp path
  (colliding os.replace -> FileNotFoundError, or a truncated file -> load()
  raises Sorry)."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="fake", model="fake", title="Race")
    conv.append(Message(role="user",
                        content=[ContentBlock(type="text",
                                              data={"text": "hi"})],
                        timestamp=now()))
    storage.save(conv)
    # A few other conversations so the index has real content to rewrite.
    for i in range(3):
      storage.save(Conversation.new(profile_name="fake", model="fake",
                                    title="c%d" % i))

    errors = []
    done = threading.Event()

    def _saver():
      try:
        for _ in range(300):
          storage.save(conv)
      except Exception as e:
        errors.append(("save", repr(e)))

    def _reader():
      try:
        while not done.is_set():
          storage.load(conv.meta.id)
          storage.list_conversations()
      except Exception as e:
        errors.append(("read", repr(e)))

    savers = [threading.Thread(target=_saver) for _ in range(2)]
    readers = [threading.Thread(target=_reader) for _ in range(2)]
    for t in savers + readers:
      t.start()
    for t in savers:
      t.join(timeout=30)
    done.set()
    for t in readers:
      t.join(timeout=10)

    assert not errors, "concurrent corruption: %r" % (errors[:5],)
    loaded = storage.load(conv.meta.id)
    assert loaded.meta.id == conv.meta.id
    assert len(loaded.messages) == 1
    assert conv.meta.id in {m.id for m in storage.list_conversations()}
  finally:
    shutil.rmtree(tmp)


def exercise_save_reindex_false_skips_index():
  """save(conv, reindex=False) writes the conversation's own files but does not
  touch index.json, so list_conversations (which reads index.json) doesn't see
  it until a normal reindexing save."""
  tmp, storage = _new_storage()
  try:
    c1 = Conversation.new(profile_name="fake", model="fake")
    storage.save(c1)                              # reindex=True -> index has c1
    c2 = Conversation.new(profile_name="fake", model="fake")
    storage.save(c2, reindex=False)               # c2 dir written, index not
    ids = {m.id for m in storage.list_conversations()}
    assert c1.meta.id in ids
    assert c2.meta.id not in ids                   # not indexed yet
    assert (storage.conv_dir(c2.meta.id) / "meta.json").exists()
    storage.save(c2)                              # reindex=True picks it up
    assert c2.meta.id in {m.id for m in storage.list_conversations()}
  finally:
    shutil.rmtree(tmp)


def exercise_dir_without_meta_is_not_listed():
  """A conversation dir with other files but no meta.json (the mid-turn state:
  only job_history.json, before any save) is skipped by list_conversations --
  the reported empty-sidebar case. A save creates meta.json and it lists."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="fake", model="fake", title="Mid")
    conv_dir = storage.conv_dir(conv.meta.id)
    conv_dir.mkdir(parents=True, exist_ok=True)
    (conv_dir / "job_history.json").write_text("{}", encoding="utf-8")
    assert conv.meta.id not in {m.id for m in storage.list_conversations()}
    storage.save(conv)
    assert conv.meta.id in {m.id for m in storage.list_conversations()}
  finally:
    shutil.rmtree(tmp)


def exercise_stale_tmp_files_are_swept():
  """Orphaned atomic-write temp files -- left by a hard kill: a unique
  ``<name>.<pid>.<uuid>.tmp`` that nothing reuses -- are reaped by the explicit,
  once-per-instance ``sweep_stale_tmp``, tree-wide (a conversation dir OR a
  deeper attachments subdir), so unique names don't accumulate without bound. A
  just-created temp (an in-flight concurrent write) is left untouched. A plain
  ``save`` does NOT sweep: the full-tree walk is deferred off the hot save path
  (the GUI schedules sweep_stale_tmp after the window is up)."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="fake", model="fake")
    storage.save(conv)                             # creates the tree
    conv_dir = storage.conv_dir(conv.meta.id)
    old = time.time() - 10 * 24 * 3600             # 10 days old: past any cutoff
    stale = conv_dir / ("messages.json.%d.deadbeef.tmp" % os.getpid())
    stale.write_text("partial", encoding="utf-8")
    os.utime(stale, (old, old))
    # An orphan in a deeper attachments subdir: only the full-tree walk reaches
    # it, so this pins that sweep_stale_tmp keeps whole-tree coverage.
    att_dir = conv_dir / "attachments"
    att_dir.mkdir(parents=True, exist_ok=True)
    stale_att = att_dir / ("blob.%d.feed.tmp" % os.getpid())
    stale_att.write_text("partial", encoding="utf-8")
    os.utime(stale_att, (old, old))
    fresh = conv_dir / ("messages.json.%d.cafebabe.tmp" % os.getpid())
    fresh.write_text("in-flight", encoding="utf-8")  # current mtime
    # A plain save does NOT sweep now -- the walk is deferred to sweep_stale_tmp.
    storage2 = ConversationStorage(Path(tmp), log=null_out())
    storage2.save(Conversation.new(profile_name="fake", model="fake"))
    assert stale.exists(), "a plain save must not sweep (deferred to sweep_stale_tmp)"
    # The explicit deferred sweep reaps stale orphans tree-wide, keeps in-flight.
    storage2.sweep_stale_tmp()
    assert not stale.exists(), "stale conv-dir orphan should be swept"
    assert not stale_att.exists(), "stale attachment orphan should be swept (full tree)"
    assert fresh.exists(), "in-flight tmp must be left for its writer"
    # The sweep runs at most once per instance: a stale temp planted AFTER the
    # first sweep survives a second sweep (the _swept_stale_tmp guard keeps the
    # full-tree walk to a single run per storage).
    stale2 = conv_dir / ("meta.json.%d.feedface.tmp" % os.getpid())
    stale2.write_text("partial", encoding="utf-8")
    os.utime(stale2, (old, old))
    storage2.sweep_stale_tmp()
    assert stale2.exists(), "second sweep must be a no-op (once per instance)"
  finally:
    shutil.rmtree(tmp)


def exercise_save_trims_trailing_unanswered_tool_use():
  """save() never persists a transcript ending in an assistant tool_use with no
  following tool_result -- an un-resumable provider 400 on the API backends. The
  trim lives in save() so EVERY writer (mid-turn autosave and the GUI
  turn-end/close/rename saves) is orphan-safe. A well-formed transcript is
  written verbatim."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="fake", model="fake")
    conv.append(Message(role="user",
                        content=[ContentBlock(type="text", data={"text": "go"})],
                        timestamp=now()))
    conv.append(Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="text", data={"text": "calling"}),
      ContentBlock(type="tool_use",
                   data={"id": "t1", "name": "echo", "input": {}})]))
    storage.save(conv)                             # orphan tail must be trimmed
    loaded = storage.load(conv.meta.id)
    assert loaded.messages, "committed prefix must survive"
    last = loaded.messages[-1]
    assert not (last.role == "assistant"
                and any(b.type == "tool_use" for b in last.content)), \
      [(m.role, [b.type for b in m.content]) for m in loaded.messages]
    # A well-formed transcript (ending in a tool_result) is written verbatim.
    conv.append(Message(role="user", timestamp=now(), content=[
      ContentBlock(type="tool_result",
                   data={"tool_use_id": "t1", "content": [], "is_error": False})]))
    storage.save(conv)
    loaded2 = storage.load(conv.meta.id)
    assert len(loaded2.messages) == 3, len(loaded2.messages)
    assert loaded2.messages[-1].content[0].type == "tool_result"
  finally:
    shutil.rmtree(tmp)


def exercise_save_trims_trailing_empty_assistant():
  """save() never persists a transcript ending in an EMPTY assistant message
  (content=[]) -- the placeholder run_turn appends before it yields an
  error-before-content (e.g. an auth give-up). Persisting it makes the next send
  post an empty-content assistant to an API backend, a non-recoverable 400. The
  trim lives in save()'s persistable_prefix choke point so EVERY writer
  (turn-end / error / close / autosave) is safe. A non-empty assistant tail is
  written verbatim.

  Revert-proof: with the empty-assistant trim removed, the reloaded transcript
  ends in the assistant placeholder and the role assertion fails."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="fake", model="fake")
    conv.append(Message(role="user",
                        content=[ContentBlock(type="text", data={"text": "go"})],
                        timestamp=now()))
    conv.append(Message(role="assistant", timestamp=now(), content=[]))  # placeholder
    storage.save(conv)                             # empty tail must be trimmed
    loaded = storage.load(conv.meta.id)
    assert [m.role for m in loaded.messages] == ["user"], \
      [(m.role, [b.type for b in m.content]) for m in loaded.messages]
    # A non-empty assistant tail survives verbatim.
    conv.messages.pop()                            # drop the placeholder
    conv.append(Message(role="assistant", timestamp=now(), content=[
      ContentBlock(type="text", data={"text": "hi"})]))
    storage.save(conv)
    loaded2 = storage.load(conv.meta.id)
    assert [m.role for m in loaded2.messages] == ["user", "assistant"], \
      [m.role for m in loaded2.messages]
    assert loaded2.messages[-1].content[0].data["text"] == "hi"
  finally:
    shutil.rmtree(tmp)


def exercise_recency_key_normalizes_null_and_naive():
  """recency_key sorts metas by updated_at, robust to a hand-edited/legacy meta
  whose updated_at is null (-> min, sorts last) or naive (-> UTC-coerced), so
  neither raises TypeError when compared against tz-aware timestamps."""
  from datetime import datetime, timezone
  from qttbx.widgets.chat.agent.conversation import (
    ConversationMeta, now, recency_key)

  def meta(cid, updated_at):
    return ConversationMeta(id=cid, title=cid, profile_name="t", model="m",
                            created_at=now(), updated_at=updated_at)

  aware_new = now()
  aware_old = datetime(2000, 1, 1, tzinfo=timezone.utc)
  naive = datetime(2010, 1, 1)                        # hand-edited, no tzinfo
  metas = [meta("null", None), meta("new", aware_new),
           meta("old", aware_old), meta("naive", naive)]
  ordered = [m.id for m in sorted(metas, key=recency_key, reverse=True)]
  assert ordered[0] == "new", ordered                # newest first
  assert ordered[-1] == "null", ordered              # null sorts last
  # naive (2010) sorts above old (2000) without raising:
  assert ordered.index("naive") < ordered.index("old"), ordered


def _spawn_live_process():
  """A real, live child process whose PID stands in for 'held by another live
  local process' -- portable, unlike PID 1 (which doesn't exist on Windows,
  where the fixed liveness probe would correctly read it as dead). Caller uses
  ``.pid``, then ``kill()`` + ``wait()``."""
  import subprocess
  import sys
  return subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])


def exercise_conversation_open_lock():
  """acquire/release/is_locked_by_other model a per-conversation open marker so
  two phenix.chat processes don't both drive (and clobber) one conversation. A
  marker naming our own PID or a dead PID is NOT a lock; a live foreign PID is,
  and acquire is atomic (can't steal a live lock) but steals a stale one. A
  marker from another host is treated as locked -- its PID can't be probed."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    cid = conv.meta.id
    lock = storage.conv_dir(cid) / ".open"
    assert storage.is_conversation_locked_by_other(cid) is False   # no marker
    assert storage.acquire_conversation_lock(cid) == LOCK_ACQUIRED  # we take it
    assert lock.read_text().strip().startswith("%d@" % os.getpid())
    assert storage.is_conversation_locked_by_other(cid) is False   # ours != other
    assert storage.acquire_conversation_lock(cid) == LOCK_ACQUIRED  # idempotent
    # A live FOREIGN pid reads as locked and cannot be stolen. Use a real child
    # process -- PID 1 doesn't exist on Windows, where the fixed probe correctly
    # reads it as dead.
    proc = _spawn_live_process()
    try:
      foreign = "%d\n" % proc.pid                    # live, same host, not ours
      lock.write_text(foreign, encoding="utf-8")
      assert storage.is_conversation_locked_by_other(cid) is True
      assert storage.acquire_conversation_lock(cid) == LOCK_HELD
      assert lock.read_text() == foreign             # untouched
      # A dead pid is stale -> not a lock, and acquire steals it.
      lock.write_text("2147483646\n", encoding="utf-8")
      assert storage.is_conversation_locked_by_other(cid) is False
      assert storage.acquire_conversation_lock(cid) == LOCK_ACQUIRED
      assert lock.read_text().strip().startswith("%d@" % os.getpid())
      storage.release_conversation_lock(cid)
      assert not lock.exists()                        # released
      # release only removes OUR marker, not a live foreign one.
      lock.write_text(foreign, encoding="utf-8")
      storage.release_conversation_lock(cid)
      assert lock.exists() and lock.read_text() == foreign
    finally:
      proc.kill()
      proc.wait()
    # A marker from ANOTHER host is treated conservatively as locked: we can't
    # probe a remote PID, so assume it is alive rather than risk a clobber, and
    # neither steal nor remove it.
    lock.write_text("1234@some-other-host\n", encoding="utf-8")
    assert storage.is_conversation_locked_by_other(cid) is True
    assert storage.acquire_conversation_lock(cid) == LOCK_HELD
    storage.release_conversation_lock(cid)
    assert lock.read_text().strip() == "1234@some-other-host"        # untouched
  finally:
    shutil.rmtree(tmp)


def exercise_conversation_lock_survives_corrupt_and_readonly():
  """The lock never bricks a launch: a non-UTF-8 or implausible-PID .open marker
  reads as 'not locked' instead of raising, and a read-only project dir (where
  the marker can't be written) still lets acquire succeed so launch / File > New
  proceed -- matching the read-only-must-still-launch behaviour of the save
  path."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    cid = conv.meta.id
    lock = storage.conv_dir(cid) / ".open"
    # Corrupt markers must not crash is_locked and (with atomic publish) are
    # classified stale -- not a lock -- regardless of age.
    lock.write_bytes(b"\xff\xfe\x00nonsense")         # non-UTF-8: no crash
    assert storage.is_conversation_locked_by_other(cid) is False
    lock.write_text("999999999999999999999@host\n", encoding="utf-8")  # huge PID
    assert storage.is_conversation_locked_by_other(cid) is False
    lock.unlink()
    if os.name == "posix" and os.getuid() != 0:
      d = storage.conv_dir(cid)
      mode = d.stat().st_mode
      os.chmod(d, 0o500)                              # r-x: can't create in dir
      try:
        assert storage.acquire_conversation_lock(cid) == LOCK_ACQUIRED  # launches
        assert not lock.exists()                       # nothing was written
      finally:
        os.chmod(d, mode)                              # restore for rmtree
  finally:
    shutil.rmtree(tmp)


def exercise_marker_state_classifier():
  """_marker_state is the single policy behind acquire / is-locked / release:
  absent, ours, live (a live local PID or ANY other-host marker), stale (a dead
  local PID or a corrupt body), or unreadable (a TRANSIENT read fault). Atomic
  publish means there is no mid-write state, so a corrupt marker is stale
  (stealable). A STRUCTURAL read error (.open is a directory) is stale too --
  recoverable -- not a permanent unreadable phantom lock."""
  import errno as _errno
  from qttbx.widgets.chat.agent.storage import (
    _marker_state, _HOSTNAME, _MARKER_ABSENT, _MARKER_OURS, _MARKER_LIVE,
    _MARKER_STALE, _MARKER_UNREADABLE)
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    lock = storage.conv_dir(conv.meta.id) / ".open"
    assert _marker_state(lock) == _MARKER_ABSENT
    lock.write_text("%d@%s\n" % (os.getpid(), _HOSTNAME), encoding="utf-8")
    assert _marker_state(lock) == _MARKER_OURS
    lock.write_text("1234@some-other-host\n", encoding="utf-8")
    assert _marker_state(lock) == _MARKER_LIVE          # another host
    lock.write_text("2147483646\n", encoding="utf-8")   # dead local PID
    assert _marker_state(lock) == _MARKER_STALE
    lock.write_bytes(b"\xff\xfe garbage")               # corrupt -> stale not live
    assert _marker_state(lock) == _MARKER_STALE
    # STRUCTURAL read error (.open is a directory) -> STALE (recoverable), NOT a
    # permanent UNREADABLE phantom lock.
    lock.unlink()
    lock.mkdir()
    assert _marker_state(lock) == _MARKER_STALE
    lock.rmdir()
    # A genuinely TRANSIENT read fault -> UNREADABLE (inconclusive).
    class _EIO:
      def read_text(self, encoding=None):
        raise OSError(_errno.EIO, "input/output error")
    assert _marker_state(_EIO()) == _MARKER_UNREADABLE
  finally:
    shutil.rmtree(tmp)


def exercise_publish_marker_is_atomic_and_cleans_temp():
  """_publish_marker writes the FULL marker (never an empty/partial file) and
  leaves no .open.tmp staging file behind; a second publish onto an existing
  marker reports 'exists' rather than overwriting."""
  from qttbx.widgets.chat.agent.storage import _HOSTNAME
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    d = storage.conv_dir(conv.meta.id)
    lock = d / ".open"
    assert storage._publish_marker(lock) is True
    assert lock.read_text().strip() == "%d@%s" % (os.getpid(), _HOSTNAME)
    assert storage._publish_marker(lock) is False       # already exists
    assert not list(d.glob(".open.*.tmp"))              # staging temp cleaned up
  finally:
    shutil.rmtree(tmp)


def exercise_conversation_lock_three_way_result():
  """acquire returns LOCK_UNAVAILABLE for a TRANSIENT write failure (ENOSPC) so
  the caller warns instead of skipping/claiming, and LOCK_ACQUIRED for a
  PERMANENT one (read-only/EROFS: nothing to coordinate, launch proceeds)."""
  import errno as _errno
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    cid = conv.meta.id
    storage._publish_marker = lambda p: OSError(_errno.ENOSPC, "no space")
    assert storage.acquire_conversation_lock(cid) == LOCK_UNAVAILABLE
    storage._publish_marker = lambda p: OSError(_errno.EROFS, "read-only")
    assert storage.acquire_conversation_lock(cid) == LOCK_ACQUIRED
  finally:
    shutil.rmtree(tmp)


def exercise_unreadable_marker_is_inconclusive_and_held():
  """A transient read fault yields _MARKER_UNREADABLE, treated as possibly-held:
  is_locked -> True, and a write error while the recheck read faults reports HELD
  (not UNAVAILABLE), so restore-on-unavailable can't restore a held conversation
  just because we couldn't read the marker."""
  import errno as _errno
  from qttbx.widgets.chat.agent import storage as _sm
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    cid = conv.meta.id
    class _EIOPath:                                   # a path whose read faults
      def read_text(self, encoding=None):
        raise OSError(_errno.EIO, "input/output error")
    assert _sm._marker_state(_EIOPath()) == _sm._MARKER_UNREADABLE
    orig = _sm._marker_state
    _sm._marker_state = lambda path: _sm._MARKER_UNREADABLE
    try:
      assert storage.is_conversation_locked_by_other(cid) is True
      storage._publish_marker = lambda p: OSError(_errno.EIO, "io")
      assert storage.acquire_conversation_lock(cid) == LOCK_HELD
      # Also the MAIN classify path (marker exists -> UNREADABLE -> HELD), not
      # only the write-error guard: a real planted marker makes publish return
      # False, and the classify still sees the stubbed UNREADABLE.
      del storage._publish_marker            # restore the real method
      (storage.conv_dir(cid) / ".open").write_text(
        "1@some-other-host\n", encoding="utf-8")
      assert storage.acquire_conversation_lock(cid) == LOCK_HELD
    finally:
      _sm._marker_state = orig
  finally:
    shutil.rmtree(tmp)


def exercise_publish_error_with_live_holder_reports_held():
  """A marker write error while another live process holds the conversation must
  report LOCK_HELD, not LOCK_UNAVAILABLE -- otherwise _startup's
  restore-on-unavailable would restore a HELD conversation (two writers)."""
  import errno as _errno
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    cid = conv.meta.id
    (storage.conv_dir(cid) / ".open").write_text(     # a live foreign holder
      "1@some-other-host\n", encoding="utf-8")
    storage._publish_marker = lambda p: OSError(_errno.ENOSPC, "no space")
    assert storage.acquire_conversation_lock(cid) == LOCK_HELD
  finally:
    shutil.rmtree(tmp)


def exercise_clear_conversation_lock_removes_any_marker():
  """clear_conversation_lock backs the user's explicit 'Unlock': it removes the
  marker even when a live foreign holder on another host has it, leaving the
  conversation lock-free (not claimed for us)."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    cid = conv.meta.id
    lock = storage.conv_dir(cid) / ".open"
    lock.write_text("1@some-other-host\n", encoding="utf-8")  # held elsewhere
    assert storage.is_conversation_locked_by_other(cid) is True
    storage.clear_conversation_lock(cid)
    assert not lock.exists()                                  # marker gone
    assert storage.is_conversation_locked_by_other(cid) is False   # lock-free
  finally:
    shutil.rmtree(tmp)


def exercise_own_lock_intact_detects_a_cleared_marker():
  """own_lock_intact is True only while OUR marker is present -- it goes False
  when the marker is cleared (an Unlock elsewhere) or replaced by a foreign one.
  This is the signal a holder polls to know it must republish."""
  tmp, storage = _new_storage()
  try:
    conv = Conversation.new(profile_name="p", model="m")
    storage.save(conv)
    cid = conv.meta.id
    assert storage.acquire_conversation_lock(cid) == LOCK_ACQUIRED
    assert storage.own_lock_intact(cid) is True                # ours
    storage.clear_conversation_lock(cid)                       # unlocked elsewhere
    assert storage.own_lock_intact(cid) is False               # marker gone
    (storage.conv_dir(cid) / ".open").write_text(
      "1@some-other-host\n", encoding="utf-8")                 # foreign marker
    assert storage.own_lock_intact(cid) is False
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
  exercise_atomic_replace_fails_fast_on_posix()
  exercise_all_token_usage_fields_round_trip()
  exercise_read_json_retries_windows_sharing_violation()
  exercise_concurrent_saves_do_not_corrupt()
  exercise_save_reindex_false_skips_index()
  exercise_dir_without_meta_is_not_listed()
  exercise_stale_tmp_files_are_swept()
  exercise_save_trims_trailing_unanswered_tool_use()
  exercise_save_trims_trailing_empty_assistant()
  exercise_conversation_open_lock()
  exercise_conversation_lock_survives_corrupt_and_readonly()
  exercise_marker_state_classifier()
  exercise_publish_marker_is_atomic_and_cleans_temp()
  exercise_conversation_lock_three_way_result()
  exercise_publish_error_with_live_holder_reports_held()
  exercise_unreadable_marker_is_inconclusive_and_held()
  exercise_clear_conversation_lock_removes_any_marker()
  exercise_own_lock_intact_detects_a_cleared_marker()
  exercise_recency_key_normalizes_null_and_naive()
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
