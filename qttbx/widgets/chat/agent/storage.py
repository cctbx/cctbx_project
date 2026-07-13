"""Per-project conversation storage.

Layout::

    <project_dir>/.phenix_chat/
      index.json
      conversations/<conv_id>/
        meta.json
        messages.json
        attachments/sha256-<hex>.<ext>
        subagents/<sub_id>.json

Atomic writes via tmp + rename. Lazy directory creation. Content-addressed
attachments (sha256-keyed; automatic dedup within a conversation).
"""

import errno
import hashlib
import json
import mimetypes
import os
import socket
import sys
import time
import uuid
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

from qttbx.widgets.chat.agent.conversation import (
  Attachment, ContentBlock, Conversation, ConversationMeta, Message,
  SubagentRecord, TokenUsage)
from qttbx.widgets.chat.agent.paths import chat_root_for


_SCHEMA_VERSION = "1.0"

# An atomic-write temp file older than this is treated as an orphan left by a
# hard kill and swept (see ConversationStorage.sweep_stale_tmp). Well above any
# real write, so a concurrent writer's in-flight temp is never hit.
_TMP_STALE_SECONDS = 3600


class ConversationStorage:
  """Read and write conversations under a project's ``.phenix_chat/`` dir.

  Construction does not touch the filesystem. Directories appear lazily on
  first write.

  Parameters
  ----------
  project_dir : str or pathlib.Path
      Project directory whose chat root holds the conversations.
  log : file-like, optional
      Stream for diagnostic messages. Defaults to ``sys.stdout``.
  """

  def __init__(self, project_dir, log=None):
    self.project_dir = Path(project_dir)
    self.root = chat_root_for(self.project_dir)
    self.log = log if log is not None else sys.stdout
    # Orphaned atomic-write temps are swept once per instance by an explicit,
    # deferred sweep_stale_tmp() (the GUI schedules it off the constructor).
    # __init__ stays filesystem-free.
    self._swept_stale_tmp = False

  # ---- conversations -------------------------------------------------------

  def list_conversations(self):
    """Return conversation metadata, rebuilding the index if needed.

    Reads the cached ``index.json``, or rebuilds it from the on-disk
    conversations when the index is missing or corrupt.

    Returns
    -------
    list of ConversationMeta
        Metadata for the stored conversations.
    """
    index_path = self.root / "index.json"
    if index_path.exists():
      try:
        with open(index_path, encoding="utf-8") as fh:
          data = json.load(fh)
        return [_meta_from_dict(d) for d in data.get("conversations", [])]
      except Exception as e:
        print("storage: index.json unreadable (%s); rebuilding"
              % e, file=self.log)
    return self._rebuild_index()

  def load(self, conv_id):
    """Load a conversation's meta and messages.

    Attachments and subagents are loaded on demand, not here.

    Parameters
    ----------
    conv_id : str
        Identifier of the conversation to load.

    Returns
    -------
    Conversation
        The conversation with ``meta`` and ``messages`` populated and
        empty ``attachments`` / ``subagents``.

    Raises
    ------
    libtbx.utils.Sorry
        If the conversation directory does not exist, or a document's
        ``schema_version`` is unsupported.
    """
    from libtbx.utils import Sorry
    conv_dir = self._conv_dir(conv_id)
    if not conv_dir.exists():
      raise Sorry("Conversation not found: %s" % conv_id)
    try:
      meta_doc = _read_json(conv_dir / "meta.json")
      _check_schema_version(meta_doc, str(conv_dir / "meta.json"))
      meta = _meta_from_dict(meta_doc)
      messages_doc = _read_json(conv_dir / "messages.json")
      _check_schema_version(messages_doc, str(conv_dir / "messages.json"))
      messages = [_message_from_dict(m)
                  for m in messages_doc.get("messages", [])]
    except Sorry:
      # Already a clear user-facing error (e.g. unsupported schema_version);
      # surface it unchanged.
      raise
    except (ValueError, KeyError, TypeError, AttributeError, OSError) as e:
      # Corrupt/partial JSON (json.JSONDecodeError is a ValueError), a
      # missing required key (KeyError), a bad timestamp (ValueError), a
      # wrong-typed field (TypeError, e.g. a non-iterable 'usage';
      # AttributeError, e.g. a content block that is a bare string), or a
      # missing/unreadable file (OSError) -- map all to the documented
      # Sorry-only contract instead of leaking a cryptic raw exception.
      raise Sorry("Conversation '%s' could not be loaded; its files are "
                  "missing or corrupt (%s)." % (conv_id, e))
    # Attachments and subagents loaded on demand.
    return Conversation(meta=meta, messages=messages,
                        attachments={}, subagents=[])

  def save(self, conv, reindex=True):
    """Write ``meta.json`` and ``messages.json`` atomically, optionally reindex.

    A transcript that would end in an un-persistable tail -- an assistant
    ``tool_use`` with no following ``tool_result``, or an empty-content assistant
    -- is trimmed to the well-formed prefix first (see ``persistable_prefix``),
    so no caller -- the worker's mid-turn autosave or the GUI turn-end / error /
    close / rename saves -- can freeze an un-resumable conversation on disk. The
    live message list is snapshotted before trimming so a GUI-thread save can run
    while the worker is still appending this turn without racing it.

    Parameters
    ----------
    conv : Conversation
        The conversation to persist.
    reindex : bool, optional
        When True (default) refresh ``index.json`` after writing. Mid-turn
        autosaves pass ``reindex=False`` to keep the O(N) index rescan (and the
        extra ``index.json`` write) out of the hot loop. The conversation is
        already listed, so it never appears/disappears mid-turn; the one field
        that does drift is its cached ``updated_at`` in ``index.json``, which
        lags the freshly written ``meta.json`` (each ``append`` bumps
        ``updated_at``) until the next reindexing save at turn end / close /
        rename. Accepted, self-healing trade-off: the sole visible effect is
        that a hard kill mid-turn can leave this conversation ordered by its
        previous turn-end time on the next launch (sidebar order and startup
        auto-restore) until a later reindexing save heals it -- nothing is lost.
    """
    self._ensure_root()
    conv_dir = self._conv_dir(conv.meta.id)
    conv_dir.mkdir(parents=True, exist_ok=True)
    _atomic_write_json(conv_dir / "meta.json", _meta_to_dict(conv.meta))
    # Snapshot the messages up front, THEN trim: a GUI-thread save (turn-end /
    # error / close) can run while the worker thread is still appending this turn
    # -- a YIELDED error leaves run_turn mid-flight with is_busy() still true, and
    # the session drives the SAME Conversation object the window hands to save().
    # persistable_prefix returns the live list when nothing needs trimming, so
    # without a private copy the serialization below would race conv.append and
    # could persist an orphaned/empty tail. list(...) is an atomic snapshot under
    # the GIL; the trim then drops any un-persistable tail (see persistable_prefix)
    # so no writer freezes an un-resumable transcript on disk.
    messages = persistable_prefix(list(conv.messages))
    _atomic_write_json(conv_dir / "messages.json", {
      "schema_version": _SCHEMA_VERSION,
      "messages": [_message_to_dict(m) for m in messages],
    })
    if reindex:
      self._refresh_index()

  def delete(self, conv_id):
    """Delete a conversation's on-disk directory, then reindex.

    Removes the conversation's directory (its ``meta.json`` /
    ``messages.json`` / ``attachments`` / ``subagents``) so it disappears
    from ``list_conversations``, then rewrites the index from the remaining
    on-disk conversations.

    Parameters
    ----------
    conv_id : str
        Identifier of the conversation to delete.

    Raises
    ------
    libtbx.utils.Sorry
        If ``conv_id`` is unsafe as a path segment, or no such conversation
        exists (matching ``load``'s not-found contract).
    """
    import shutil
    from libtbx.utils import Sorry
    # _conv_dir validates conv_id via _safe_segment (rejects path
    # separators / '..' / absolute ids), so the join cannot escape the
    # chat root.
    conv_dir = self._conv_dir(conv_id)
    if not conv_dir.exists():
      raise Sorry("Conversation not found: %s" % conv_id)
    # Defense in depth before a destructive rmtree: confirm the resolved
    # target really sits directly under this project's conversations dir
    # (a conv id that resolves elsewhere -- e.g. through a symlink -- is
    # refused rather than followed).
    conv_root = (self.root / "conversations").resolve()
    resolved = conv_dir.resolve()
    if conv_root not in resolved.parents:
      raise Sorry("Refusing to delete outside the chat root: %s" % conv_id)
    shutil.rmtree(resolved)
    self._refresh_index()

  # ---- attachments ---------------------------------------------------------

  def store_attachment(self, conv_id, data, mime):
    """Store attachment bytes content-addressed by their sha256.

    Idempotent: writing the same bytes twice produces one file.

    Parameters
    ----------
    conv_id : str
        Conversation the attachment belongs to.
    data : bytes
        The attachment bytes.
    mime : str
        MIME type, used to choose the file extension.

    Returns
    -------
    Attachment
        Reference to the stored bytes, keyed by their sha256.
    """
    self._ensure_root()
    sha = hashlib.sha256(data).hexdigest()
    ext = mimetypes.guess_extension(mime) or ".bin"
    fname = "sha256-%s%s" % (sha, ext)
    att_dir = self._conv_dir(conv_id) / "attachments"
    att_dir.mkdir(parents=True, exist_ok=True)
    fpath = att_dir / fname
    if not fpath.exists():
      _atomic_write_bytes(fpath, data)
    return Attachment(sha256=sha, mime=mime, path=fname)

  def load_attachment(self, conv_id, sha256):
    """Load attachment bytes by their sha256.

    Parameters
    ----------
    conv_id : str
        Conversation the attachment belongs to.
    sha256 : str
        Hex sha256 of the attachment bytes.

    Returns
    -------
    bytes
        The stored attachment bytes.

    Raises
    ------
    libtbx.utils.Sorry
        If ``sha256`` is unsafe as a path segment, or no matching
        attachment exists.
    """
    import glob as _glob
    sha256 = _safe_segment(sha256, "attachment")
    att_dir = self._conv_dir(conv_id) / "attachments"
    # _glob.escape neutralizes any glob metacharacters in the (validated)
    # sha so it is matched literally, not as a pattern.
    matches = sorted(att_dir.glob("sha256-%s.*" % _glob.escape(sha256)))
    if not matches:
      from libtbx.utils import Sorry
      raise Sorry("Attachment not found: sha256=%s" % sha256)
    with open(matches[0], "rb") as fh:
      return fh.read()

  # ---- subagents -----------------------------------------------------------

  def store_subagent(self, conv_id, record):
    """Write a subagent record atomically under the conversation."""
    self._ensure_root()
    sub_dir = self._conv_dir(conv_id) / "subagents"
    sub_dir.mkdir(parents=True, exist_ok=True)
    path = sub_dir / ("%s.json" % _safe_segment(record.sub_id, "subagent"))
    _atomic_write_json(path, _subagent_to_dict(record))

  def load_subagent(self, conv_id, sub_id):
    """Load a subagent record by its id.

    Parameters
    ----------
    conv_id : str
        Conversation the subagent belongs to.
    sub_id : str
        Identifier of the subagent record.

    Returns
    -------
    SubagentRecord
        The loaded subagent record.

    Raises
    ------
    libtbx.utils.Sorry
        If ``sub_id`` is unsafe as a path segment, or no record exists.
    """
    path = (self._conv_dir(conv_id) / "subagents"
            / ("%s.json" % _safe_segment(sub_id, "subagent")))
    if not path.exists():
      from libtbx.utils import Sorry
      raise Sorry("Subagent record not found: %s" % sub_id)
    return _subagent_from_dict(_read_json(path))

  def conv_dir(self, conv_id):
    """Public path to a conversation's directory (creates nothing)."""
    return self._conv_dir(conv_id)

  # ---- open marker (concurrent-launch safety) ------------------------------

  def _lock_path(self, conv_id):
    return self._conv_dir(conv_id) / ".open"

  def acquire_conversation_lock(self, conv_id):
    """Try to publish this process's ``.open`` marker for a conversation.

    Two phenix.chat processes on one project must not both drive a single
    conversation: ``save`` is a whole-file rewrite, so the last to finish a turn
    silently clobbers the other's messages. The marker holds ``<pid>@<hostname>``
    and is published ATOMICALLY -- written to a unique temp, then ``os.link``ed
    into place -- so it never appears empty or partially written, even over NFS
    (no mid-write window, no grace heuristic). A marker left by a crashed process
    on THIS host (dead PID) is stolen; one held by another live local process, or
    by any process on a DIFFERENT host (whose PID we cannot check), is not.

    BEST-EFFORT on network filesystems. The steal path is classify-then-rename,
    with no atomic compare-and-swap on NFS/SMB, so a transient EIO/ESTALE at the
    wrong instant plus a concurrent holder can still let two writers through a
    narrow window -- an inherent limit, not a bug to keep chasing. The real
    data-loss vector is the whole-file-clobbering ``save``; the backstop is that
    the lock is advisory -- a holder republishes a marker cleared from under it,
    and the explicit Unlock frees a genuinely dead holder's marker.
    Eliminating it would take a non-clobbering ``save`` (concurrent-modification
    detection); accepted as-is -- this residual network-FS race is out of scope
    for the concurrent-session workflows we support.

    Likewise on WINDOWS a concurrent marker read (another window's 2s sidebar
    lock scan) can briefly block the stale-marker steal (``os.rename`` below) or a
    ``release`` (``os.unlink``) -- CPython opens without ``FILE_SHARE_DELETE`` --
    surfacing as a transient ``LOCK_UNAVAILABLE`` or a momentarily un-cleared
    marker, both healed by the next poll / acquire. (The whole-file ``save`` /
    autosave replace is retried instead; see ``_atomic_replace``.)

    Returns
    -------
    str
        ``LOCK_ACQUIRED`` -- we now hold it (or already did, or the dir is
        permanently unwritable, so there is nothing to coordinate and launch
        must still proceed); ``LOCK_HELD`` -- another live process holds it (here
        or on another host); ``LOCK_UNAVAILABLE`` -- the marker could not be
        written for a TRANSIENT reason (ENOSPC/EIO/ESTALE/...), so the caller
        should proceed with a warning rather than skip or silently claim it.
    """
    p = self._lock_path(conv_id)
    try:
      self._conv_dir(conv_id).mkdir(parents=True, exist_ok=True)
    except OSError as e:
      print("storage: open-marker dir unavailable for %s (%s)"
            % (conv_id, e), file=self.log)
      return LOCK_ACQUIRED if _permanent_errno(e) else LOCK_UNAVAILABLE
    steal_blocked = False                  # a non-ENOENT steal error occurred
    for _ in range(3):                     # retries after stealing a stale marker
      published = self._publish_marker(p)
      if published is True:
        return LOCK_ACQUIRED
      if published is not False:           # an OSError -> couldn't write it
        # A live holder trumps a write error: report HELD, not UNAVAILABLE, so
        # restore-on-unavailable can't restore a held conversation. The recheck
        # read can ITSELF fault (EIO/ESTALE) -- UNREADABLE is inconclusive, so
        # treat it as possibly-held too rather than fire UNAVAILABLE blind.
        if _marker_state(p) in _HELD_STATES:
          return LOCK_HELD
        print("storage: could not publish open marker for %s (%s)"
              % (conv_id, published), file=self.log)
        return (LOCK_ACQUIRED if _permanent_errno(published)
                else LOCK_UNAVAILABLE)
      # published is False: the marker already exists -- inspect the holder.
      state = _marker_state(p)
      if state == _MARKER_OURS:
        return LOCK_ACQUIRED
      if state in _HELD_STATES:
        return LOCK_HELD                   # held, or a read we couldn't resolve
      if state == _MARKER_ABSENT:
        continue                           # vanished between link and read; retry
      # _MARKER_STALE (dead PID or corrupt): steal atomically. Rename aside so of
      # two racers only the one whose rename succeeds proceeds.
      steal = _unique_tmp(p, "steal.")
      try:
        os.rename(str(p), str(steal))
      except FileNotFoundError:
        continue                           # lost the race; someone moved it
      except OSError:
        steal_blocked = True               # can't rename (AV / perm), not a lock
        continue
      # We renamed something aside. If it turned out NOT stale -- a peer
      # republished a fresh (live) marker in the gap and we stole IT -- put it
      # back and re-loop, which then sees the holder. If that restore itself
      # fails transiently we drop the aside copy: the peer's lock is already off
      # `p` the instant the rename succeeded, so no restore path can fully
      # recover it -- inherent to classify-then-rename without atomic CAS (see
      # the best-effort note on acquire); not chased further.
      if _marker_state(steal) != _MARKER_STALE:
        try:
          os.replace(str(steal), str(p))
        except OSError:
          _unlink_quietly(steal)
        continue
      _unlink_quietly(steal)
      continue                             # genuine stale steal; retry publish
    return LOCK_UNAVAILABLE if steal_blocked else LOCK_HELD

  def _publish_marker(self, p):
    """Atomically publish our marker at ``p`` with full content: write it to a
    unique temp, then ``os.link`` the temp onto ``p``. The target never appears
    empty or partial (NFS-safe). Returns ``True`` (published), ``False`` (``p``
    already exists -- someone holds it), or the ``OSError`` if the marker could
    not be written.

    Falls back to ``O_EXCL`` create+write on ANY ``os.link`` error except
    ``FileExistsError`` -- a link-less or hardlink-denying filesystem (SMB/NAS
    EACCES, FAT/exFAT ENOTSUP/EINVAL) where the temp write just succeeded, so
    'nothing to coordinate' would be false. The fallback's own outcome is
    returned."""
    marker = "%d@%s\n" % (os.getpid(), _HOSTNAME)
    tmp = _unique_tmp(p)
    try:
      tmp.write_text(marker, encoding="utf-8")   # open+write+close; errors caught
    except OSError as e:
      _unlink_quietly(tmp)
      return e
    try:
      os.link(str(tmp), str(p))            # atomic create-with-content
      return True
    except FileExistsError:
      return False
    except OSError:                        # link unsupported/denied -> O_EXCL
      return self._publish_marker_oexcl(p, marker)
    finally:
      _unlink_quietly(tmp)                 # the temp is only a staging file

  def _publish_marker_oexcl(self, p, marker):
    """``O_EXCL`` create+write fallback for filesystems that can't (or won't)
    hardlink. Carries the small empty-marker window ``os.link`` avoids, so it is
    used only when ``os.link`` is unsupported. Same return contract as
    ``_publish_marker``.

    The fd is ALWAYS closed -- no leak even when the write raises, or the
    turn-end retry would walk a full disk to EMFILE and break the very saves the
    lock protects. Ownership caveat: on a failed write the just-created marker is
    still empty, so a racing steal can't be told from our own file; the cleanup
    unlink is accepted as-is on this FAT-tier fallback path."""
    try:
      fd = os.open(str(p), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
      return False
    except OSError as e:
      return e
    err = None
    try:
      os.write(fd, marker.encode("utf-8"))
    except OSError as e:
      err = e
    try:
      os.close(fd)                         # always; can also defer ENOSPC/EIO
    except OSError as e:
      err = err or e
    if err is not None:
      _unlink_quietly(p)                   # drop our partial marker (see caveat)
      return err
    return True

  def clear_conversation_lock(self, conv_id):
    """Unconditionally remove this conversation's ``.open`` marker, whoever holds
    it -- the user's explicit 'Unlock', taken knowing the conversation may be
    live in another window (both would then write, risking corruption). The
    conversation is left LOCK-FREE; this window re-locks it only if the user
    switches to it. Best-effort; never raises."""
    _unlink_quietly(self._lock_path(conv_id))

  def release_conversation_lock(self, conv_id):
    """Remove this process's ``.open`` marker (leaves one owned by another live
    process, or one on a different host, untouched). Best-effort; never
    raises."""
    p = self._lock_path(conv_id)
    try:
      # Narrow accepted race: if another window unlocks this conversation and a
      # third re-acquires it between this classify and the unlink, we could
      # delete that fresh marker. Left as-is -- sub-millisecond, and re-locking
      # on the next switch/poll heals it. (_marker_state == OURS implies the
      # file exists; a vanish before unlink is caught below.)
      if _marker_state(p) == _MARKER_OURS:
        p.unlink()
    except OSError:
      pass

  def is_conversation_locked_by_other(self, conv_id):
    """True if another LIVE process holds this conversation's ``.open`` marker --
    a live local PID, any marker from a different host (whose PID we cannot
    probe), or one we couldn't read (transient fault -- treated as possibly
    held). Our own PID, a dead PID (stale after a crash), a corrupt marker, or no
    marker is NOT a lock."""
    return _marker_state(self._lock_path(conv_id)) in _HELD_STATES

  def own_lock_intact(self, conv_id):
    """True if this conversation's ``.open`` marker still exists AND is ours. A
    holder checks this right after acquiring to learn whether it established a
    real on-disk marker to defend (a read-only dir yields ACQUIRED with no
    marker, so this is False -- there is nothing to coordinate or republish)."""
    return _marker_state(self._lock_path(conv_id)) == _MARKER_OURS

  def own_marker_needs_republish(self, conv_id):
    """True if our ``.open`` marker is neither present-and-ours NOR unreadable --
    i.e. ABSENT (cleared by an Unlock elsewhere) or held by someone else, so a
    holder that established a real marker should try to re-establish it. An
    UNREADABLE fault reading our OWN marker returns False: it is inconclusive, so
    a live session is never demoted on a transient fault reading its own marker
    (the next poll re-checks once the fault clears)."""
    return _marker_state(self._lock_path(conv_id)) not in (
      _MARKER_OURS, _MARKER_UNREADABLE)

  # ---- internal ------------------------------------------------------------

  def _ensure_root(self):
    # 0700 on the root + conversations dir denies other users traversal into
    # every conversation beneath, so the transcripts, the claude_session.jsonl,
    # and attachments are not world-readable on a shared project dir regardless
    # of per-file mode (files are also written 0600 via _atomic_write).
    self.root.mkdir(parents=True, exist_ok=True)
    _restrict(self.root, 0o700)
    convs = self.root / "conversations"
    convs.mkdir(exist_ok=True)
    _restrict(convs, 0o700)

  def sweep_stale_tmp(self):
    """Best-effort, once-per-instance removal of orphaned atomic-write temps.

    ``_atomic_write`` gives each writer a unique ``<name>.<pid>.<uuid>.tmp`` and
    removes its own on failure, but a hard kill (crash / power loss / OOM)
    between ``open`` and ``os.replace`` leaves one behind that nothing reuses --
    so unique names can accumulate across runs. A full-tree ``rglob`` reaps them
    wherever they land -- a conversation dir, or a deeper attachments / subagents
    subdir. Only temps older than ``_TMP_STALE_SECONDS`` are removed, so a
    concurrent writer's in-flight temp is never touched. Never raises --
    housekeeping must not break anything.

    DEFERRED off the save hot path: the walk cost scales with the total stored
    file count, so callers must NOT run it during a save. The GUI schedules it
    once, after the window is up (``ChatWindow.__init__`` via
    ``QTimer.singleShot``), rather than paying an unbounded ``rglob`` inside the
    first save during window construction.
    """
    if self._swept_stale_tmp:
      return
    self._swept_stale_tmp = True
    cutoff = time.time() - _TMP_STALE_SECONDS
    try:
      candidates = list(self.root.rglob("*.tmp"))    # whole tree; see docstring
    except OSError:
      return
    for p in candidates:
      try:
        if p.is_file() and p.stat().st_mtime < cutoff:
          p.unlink()
      except OSError:
        pass

  def _conv_dir(self, conv_id):
    return self.root / "conversations" / _safe_segment(conv_id, "conversation")

  def _refresh_index(self):
    """Re-derive the index from on-disk conversations and write it.

    Done as a full rewrite rather than an incremental update; the index
    is small enough that this is simpler.
    """
    self._ensure_root()
    self._write_index(self._scan_conversation_metas())

  def _rebuild_index(self):
    metas = self._scan_conversation_metas()
    if metas:
      self._ensure_root()
      self._write_index(metas)
    return metas

  def _write_index(self, metas):
    """Atomically (re)write index.json from the given conversation metas.

    Single index writer shared by _refresh_index (unconditional) and
    _rebuild_index (only when on-disk conversations exist); callers own the
    _ensure_root() precondition.
    """
    _atomic_write_json(self.root / "index.json", {
      "schema_version": _SCHEMA_VERSION,
      "conversations": [_meta_to_dict(m) for m in metas],
    })

  def _scan_conversation_metas(self):
    conv_root = self.root / "conversations"
    if not conv_root.exists():
      return []
    out = []
    for conv_dir in sorted(conv_root.iterdir()):
      meta_path = conv_dir / "meta.json"
      if not meta_path.exists():
        continue
      try:
        out.append(_meta_from_dict(_read_json(meta_path)))
      except Exception as e:
        print("storage: skipping unreadable %s (%s)"
              % (meta_path, e), file=self.log)
    return out


# ---- serialization helpers -------------------------------------------------

def _safe_segment(value, kind):
  """Validate an externally-supplied identifier before joining it into a path.

  Conversation ids, attachment sha256s and subagent ids are normally
  uuids / hex, but they reach storage from on-disk conversation files that
  may have been shared or hand-edited. Rejecting empty / ``.`` / ``..`` /
  path separators / drive markers / absolute paths keeps a crafted id from
  escaping the chat root (``pathlib`` resolves ``..`` and drops the left
  operand entirely when the right operand is absolute).

  Parameters
  ----------
  value : str
      Identifier to validate (a conversation id, attachment sha256, or
      subagent id).
  kind : str
      Human-readable category used in the error message, e.g.
      ``"conversation"``, ``"attachment"``, or ``"subagent"``.

  Returns
  -------
  str
      ``str(value)`` unchanged, once validated as a single safe path
      segment.

  Raises
  ------
  libtbx.utils.Sorry
      If ``value`` is empty, ``.`` / ``..``, absolute, or contains a path
      separator, drive marker, or embedded NUL byte.
  """
  s = str(value)
  if (not s or s in (".", "..")
      or "/" in s or "\\" in s or ":" in s or "\x00" in s
      or os.path.isabs(s) or s != os.path.basename(s)):
    from libtbx.utils import Sorry
    raise Sorry("Unsafe %s identifier: %r" % (kind, value))
  return s


def _read_json(path):
  """Read and parse a JSON file.

  On Windows the ``open`` is retried through a transient ``PermissionError``
  (ERROR_ACCESS_DENIED): a concurrent writer's ``os.replace`` -- the atomic-write
  final step -- briefly holds the target, so a reader that opens in that window
  fails, the mirror image of the ``_atomic_replace`` write-side sharing violation
  (another phenix.chat process saving, or the worker's autosave racing a
  ``load`` / sidebar scan). A short bounded backoff (up to
  ``_REPLACE_RETRY_SECONDS``) wins the gap; a genuinely unreadable path (a real
  EACCES / read-only tree) still surfaces once the deadline passes, so ``load``'s
  corrupt-file contract is unchanged. POSIX ``os.replace`` swaps atomically with
  no such window, so POSIX reads once (a POSIX EACCES is a real permission
  error)."""
  if os.name != "nt":
    with open(path, encoding="utf-8") as fh:
      return json.load(fh)
  delay = 0.001
  deadline = time.monotonic() + _REPLACE_RETRY_SECONDS
  while True:
    try:
      with open(path, encoding="utf-8") as fh:
        return json.load(fh)
    except PermissionError:
      if time.monotonic() >= deadline:
        raise
      time.sleep(delay)
      delay = min(delay * 2, 0.05)


_HOSTNAME = socket.gethostname() or "localhost"

# PIDs at/above this are implausible; a marker naming one is treated as garbage
# and never handed to os.kill (which would raise OverflowError) / ctypes.
_MAX_PID = 2 ** 31

# acquire_conversation_lock outcomes -- three-way so callers can distinguish
# "held by another" from "couldn't write the marker" (disk full / storage error,
# or a stale lock that couldn't be cleared). A read-only dir maps to
# LOCK_ACQUIRED (nothing to coordinate), NOT LOCK_UNAVAILABLE.
LOCK_ACQUIRED = "acquired"
LOCK_HELD = "held"
LOCK_UNAVAILABLE = "unavailable"

# _marker_state classifications -- the single source of truth for
# acquire / is-locked / release, so the policy can't drift across copies.
_MARKER_ABSENT = "absent"
_MARKER_OURS = "ours"
_MARKER_LIVE = "live"
_MARKER_STALE = "stale"
_MARKER_UNREADABLE = "unreadable"         # present but a transient read fault
# States meaning "don't touch this conversation" -- a live holder, or a marker
# we couldn't read (inconclusive -> possibly held). Centralised so the three
# dispatch sites (acquire's guard + main classify, is-locked) can't drift.
_HELD_STATES = frozenset((_MARKER_LIVE, _MARKER_UNREADABLE))


def _parse_marker(raw):
  """Parse a stripped ``.open`` marker string into ``(pid, host)``, or
  ``(None, None)`` if malformed or an implausible PID. New markers are
  ``<pid>@<hostname>``; a legacy bare ``<pid>`` parses to ``(pid, None)`` -- host
  unknown, treated as the local host by callers."""
  pid_str, sep, host = raw.partition("@")
  try:
    pid = int(pid_str)
  except ValueError:
    return (None, None)
  if not 0 < pid < _MAX_PID:               # implausible -> garbage, not a PID
    return (None, None)
  return (pid, host if sep else None)


def _unlink_quietly(path):
  """Best-effort unlink; ignores a missing / undeletable file."""
  try:
    os.unlink(str(path))
  except OSError:
    pass


def _unique_tmp(p, infix=""):
  """A unique sibling temp path for ``p`` that ENDS in ``.tmp`` so the stale-temp
  sweeper (``rglob('*.tmp')``) reaps crash-orphans. ``infix`` tags the kind
  (e.g. ``"steal."``). Centralised so the .tmp-suffix invariant can't be broken
  by editing one copy."""
  return p.with_name(".open.%s%d.%s.tmp"
                     % (infix, os.getpid(), uuid.uuid4().hex))


def _marker_state(path):
  """Classify a conversation's ``.open`` marker in a SINGLE read -- the one source
  of truth for acquire / is-locked / release, so the three can't drift.

  Returns ``_MARKER_ABSENT`` (no marker), ``_MARKER_OURS`` (this process),
  ``_MARKER_LIVE`` (a live local PID, or any marker from a different host whose
  PID we cannot probe), ``_MARKER_STALE`` (a dead local PID, or a corrupt body),
  or ``_MARKER_UNREADABLE`` (present but a TRANSIENT read fault -- EIO/ESTALE on
  a network mount; see _transient_errno). ``FileNotFoundError`` alone is ABSENT
  -- no second ``exists()`` probe, so a marker published in the gap isn't misread
  as STALE and stolen. UNREADABLE is INCONCLUSIVE: callers stay conservative
  (treat as possibly-held) rather than steal or clobber a marker that might be a
  live holder. A STRUCTURAL read error (EISDIR -- ``.open`` is a directory --
  ENOTDIR, ...) is NOT inconclusive: it maps to STALE so the stray is renamed
  aside and recovered, never a permanent phantom lock.

  Atomic publish (``os.link``) means our markers never appear empty or partial.
  A KNOWN residual: an older-release peer on a shared mount still uses
  create-then-write, so its momentary 0-byte marker classifies STALE and could
  be stolen -- documented rather than guarded, since round-4 showed clock-based
  grace windows are NFS-skew-fragile."""
  try:
    raw = path.read_text(encoding="utf-8").strip()
  except FileNotFoundError:
    return _MARKER_ABSENT
  except ValueError:                       # non-UTF-8 decode -> genuinely corrupt
    return _MARKER_STALE
  except OSError as e:                      # read fault: transient vs structural
    return _MARKER_UNREADABLE if _transient_errno(e) else _MARKER_STALE
  pid, host = _parse_marker(raw)
  if pid is None:                          # present but malformed -> corrupt
    return _MARKER_STALE
  if pid == os.getpid() and host in (None, _HOSTNAME):
    return _MARKER_OURS
  if host is not None and host != _HOSTNAME:
    return _MARKER_LIVE                    # another host: can't probe, assume alive
  return _MARKER_LIVE if _pid_alive(pid) else _MARKER_STALE


def _transient_errno(e):
  """True if an OSError READING the ``.open`` marker is a genuinely transient /
  inconclusive fault (a flaky network mount) -- so the marker is treated as
  possibly-held (UNREADABLE) rather than stolen. A STRUCTURAL error (EISDIR --
  ``.open`` is a directory from a botched sync; ENOTDIR; ...) is NOT transient: it
  maps to STALE so the stray can be renamed aside and recovered
  in-app. An ALLOW-list, mirroring _permanent_errno on the write
  side."""
  transient = {errno.EIO, errno.EAGAIN, errno.EINTR, errno.EBUSY}
  for name in ("ESTALE", "ETIMEDOUT", "EWOULDBLOCK", "ENETUNREACH",
               "EHOSTUNREACH", "ECONNRESET", "EREMOTEIO", "ENOLINK"):
    code = getattr(errno, name, None)
    if code is not None:
      transient.add(code)
  return e.errno in transient


def _permanent_errno(e):
  """True if an OSError writing the ``.open`` marker is a PERMANENT 'nothing to
  coordinate' condition -- a read-only tree, a permission denial, or a stray file
  where the dir should be. An ALLOW-list (not a deny-list), so an unlisted or
  transient error (ENOSPC/EIO/ESTALE/EMFILE/...) is NOT treated as permanent:
  acquire then reports LOCK_UNAVAILABLE (proceed with a warning) instead of
  silently claiming the lock and letting two writers through once the mount
  recovers."""
  permanent = (errno.EROFS, errno.EACCES, errno.EPERM,
               errno.EEXIST, errno.ENOTDIR)
  return e.errno in permanent


def _pid_alive(pid):
  """Best-effort liveness check for a PID **on this host**.

  POSIX uses ``os.kill(pid, 0)``. Windows must NOT: there ``os.kill`` maps to
  ``TerminateProcess`` and would kill the target, so a non-destructive
  ``OpenProcess`` probe is used instead. On any other platform, and whenever the
  check is inconclusive, the PID is treated as alive -- the conservative choice,
  since a false 'dead' risks two sessions clobbering one conversation whereas a
  false 'alive' only keeps a conversation locked a while longer."""
  if pid is None or pid <= 0:
    return False
  if os.name == "posix":
    try:
      os.kill(pid, 0)
    except ProcessLookupError:
      return False
    except (OSError, OverflowError, ValueError):
      return True                          # e.g. EPERM: alive, another user
    return True
  if os.name == "nt":
    return _pid_alive_windows(pid)
  return True                              # unknown platform: assume alive


def _pid_alive_windows(pid):
  """Non-destructive Windows liveness probe (never ``os.kill``, which there
  would ``TerminateProcess``). Opens the process with minimal rights; if that
  fails, distinguishes a genuinely dead PID (``OpenProcess`` sets
  ``ERROR_INVALID_PARAMETER``) from one we simply may not open
  (``ERROR_ACCESS_DENIED`` and anything else -> alive), so a crashed session's
  marker is recoverable instead of a permanent lockout.

  Uses a PRIVATE ``WinDLL(use_last_error=True)`` so ``get_last_error()`` is
  reliable and the shared ``ctypes.windll`` cache is not mutated by the
  ``argtypes`` assignments below."""
  try:
    import ctypes
    from ctypes import wintypes
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
  except Exception:
    return True
  ERROR_INVALID_PARAMETER = 87
  PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
  SYNCHRONIZE = 0x00100000
  STILL_ACTIVE = 259
  WAIT_TIMEOUT = 0x00000102
  k32.OpenProcess.restype = wintypes.HANDLE
  k32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
  k32.GetExitCodeProcess.restype = wintypes.BOOL
  k32.GetExitCodeProcess.argtypes = [wintypes.HANDLE,
                                     ctypes.POINTER(wintypes.DWORD)]
  k32.WaitForSingleObject.restype = wintypes.DWORD
  k32.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
  k32.CloseHandle.restype = wintypes.BOOL
  k32.CloseHandle.argtypes = [wintypes.HANDLE]
  handle = k32.OpenProcess(
    PROCESS_QUERY_LIMITED_INFORMATION | SYNCHRONIZE, False, pid)
  if not handle:
    # No such PID -> ERROR_INVALID_PARAMETER means dead; anything else (e.g.
    # ERROR_ACCESS_DENIED) means it exists but we can't query it -> alive.
    return ctypes.get_last_error() != ERROR_INVALID_PARAMETER
  try:
    code = wintypes.DWORD()
    if k32.GetExitCodeProcess(handle, ctypes.byref(code)):
      if code.value != STILL_ACTIVE:
        return False                       # a real exit code -> dead
      # STILL_ACTIVE (259) is ambiguous (a process may exit with 259); confirm
      # with a zero wait -- WAIT_TIMEOUT means still running.
      return k32.WaitForSingleObject(handle, 0) == WAIT_TIMEOUT
    return True
  finally:
    k32.CloseHandle(handle)


# Longest the final replace retries a Windows sharing violation before giving up.
# On Windows a concurrent reader holding the target open -- another phenix.chat
# process's load() / list_conversations() on the same project -- blocks the
# replace: CPython's open() grants FILE_SHARE_READ|WRITE but NOT
# FILE_SHARE_DELETE, so MoveFileEx can't take the delete access os.replace needs
# and raises PermissionError. The reader's hold is brief (open -> json.load ->
# close), so a short backoff wins the gap. POSIX renames over an open file
# unconditionally, so this never fires there.
_REPLACE_RETRY_SECONDS = 2.0


def _atomic_replace(tmp, path):
  """``os.replace(tmp, path)`` with a bounded retry on a Windows sharing
  violation (see ``_REPLACE_RETRY_SECONDS``). A concurrent reader's handle on
  ``path`` is held only momentarily, so a short exponential backoff almost always
  wins; a genuinely permanent ``PermissionError`` (e.g. a read-only tree) still
  surfaces once the deadline passes, preserving the prior raise-on-failure
  contract.

  POSIX renames over an open target unconditionally, so the sharing violation
  the retry defends against cannot occur there: a ``PermissionError`` on POSIX is
  a PERMANENT condition (a read-only tree / EACCES / EPERM) that will not clear
  in ``_REPLACE_RETRY_SECONDS``. So POSIX takes a fast path that re-raises at
  once -- spinning the retry would only freeze the GUI thread ~2s on every failed
  save before surfacing the same error. The retry is a Windows-only concession."""
  if os.name != "nt":
    os.replace(tmp, path)                    # POSIX: no sharing violation -> fail fast
    return
  delay = 0.001
  deadline = time.monotonic() + _REPLACE_RETRY_SECONDS
  while True:
    try:
      os.replace(tmp, path)
      return
    except PermissionError:
      if time.monotonic() >= deadline:
        raise
      time.sleep(delay)
      delay = min(delay * 2, 0.05)


def _restrict(path, mode):
  """Best-effort tighten permissions (POSIX; a near-no-op on Windows). Chat
  transcripts / session files must not be world-readable on a shared project
  dir. Swallow errors: a perms tweak must never break a save."""
  try:
    os.chmod(str(path), mode)
  except OSError:
    pass


def _atomic_write(path, mode, write_fn):
  # Write to a UNIQUE sibling .tmp via write_fn(handle), then atomically rename
  # it into place. The pid+uuid tag means two concurrent writers -- two threads
  # (a worker-thread mid-turn autosave racing the GUI's turn-end save) or even
  # two phenix.chat processes on the same project dir -- never share a tmp, so
  # each replace stays atomic and it is last-writer-wins. The replace goes through
  # _atomic_replace so a Windows reader holding the target open only delays the
  # write briefly instead of failing it. If the write (or replace) fails partway,
  # remove the temp file so a mid-write failure cannot leave an orphaned .tmp
  # behind.
  tmp = path.with_suffix(
    path.suffix + ".%d.%s.tmp" % (os.getpid(), uuid.uuid4().hex))
  replaced = False
  try:
    with open(tmp, mode) as fh:
      write_fn(fh)
    _restrict(tmp, 0o600)          # the rename carries the tightened mode in
    _atomic_replace(tmp, path)
    replaced = True
  finally:
    if not replaced:
      try:
        os.remove(tmp)
      except OSError:
        pass


def _atomic_write_json(path, obj):
  _atomic_write(path, "w",
                lambda fh: json.dump(obj, fh, indent=2, default=_json_default))


def _atomic_write_bytes(path, data):
  _atomic_write(path, "wb", lambda fh: fh.write(data))


def _json_default(o):
  if isinstance(o, datetime):
    return o.isoformat()
  # Defense in depth: a single value a backend failed to normalize to a
  # JSON-native form (e.g. a provider SDK's pydantic result object that slipped
  # through) must not make the whole conversation permanently unsaveable.
  # Coerce a pydantic model (mode="json" -> fully JSON-native) instead of
  # raising and losing every turn from here on. A genuinely unknown object
  # still raises, so a real serialization bug is surfaced rather than masked.
  dump = getattr(o, "model_dump", None)
  if callable(dump):
    return dump(mode="json")
  raise TypeError("Not JSON serializable: %r" % o)


def _parse_dt(s):
  if s is None:
    return None
  # fromisoformat handles offset-aware ISO timestamps in 3.11+.
  return datetime.fromisoformat(s)


def _check_schema_version(doc, source):
  """Validate a document's ``schema_version`` against the supported one.

  Currently v1 only (no migrations needed); this is the seam future
  migrations plug into.

  Parameters
  ----------
  doc : dict
      The loaded JSON document. Non-dict values are accepted and skipped.
  source : str
      Path or label of the document, used in the error message.

  Raises
  ------
  libtbx.utils.Sorry
      If ``schema_version`` is from a version this client cannot migrate.
  """
  if not isinstance(doc, dict):
    return
  version = doc.get("schema_version", _SCHEMA_VERSION)
  if version != _SCHEMA_VERSION:
    from libtbx.utils import Sorry
    raise Sorry(
      "%s: schema_version '%s' is not supported by this client "
      "(expected '%s'). Update PhenixChat or use a compatible version."
      % (source, version, _SCHEMA_VERSION))


def _meta_to_dict(m):
  return {
    "id": m.id,
    "title": m.title,
    "profile_name": m.profile_name,
    "model": m.model,
    "backend": m.backend,
    "created_at": m.created_at,
    "updated_at": m.updated_at,
    "archived": m.archived,
    "pinned": m.pinned,
    "summary": m.summary,
    "agent_session_id": m.agent_session_id,
    "schema_version": m.schema_version,
  }


def _meta_from_dict(d):
  return ConversationMeta(
    id=d["id"],
    title=d.get("title", ""),
    profile_name=d.get("profile_name", ""),
    model=d.get("model", ""),
    backend=d.get("backend", ""),
    created_at=_parse_dt(d["created_at"]),
    updated_at=_parse_dt(d["updated_at"]),
    archived=d.get("archived", False),
    pinned=d.get("pinned", False),
    summary=d.get("summary", ""),
    agent_session_id=d.get("agent_session_id"),
    schema_version=d.get("schema_version", _SCHEMA_VERSION),
  )


def _content_block_to_dict(b):
  # Recursively serialize nested ContentBlocks inside tool_result.content.
  data = {}
  for k, v in b.data.items():
    if isinstance(v, list) and v and isinstance(v[0], ContentBlock):
      data[k] = [_content_block_to_dict(inner) for inner in v]
    else:
      data[k] = v
  return {"type": b.type, "data": data}


def _content_block_from_dict(d):
  data = {}
  for k, v in d.get("data", {}).items():
    if isinstance(v, list) and v and isinstance(v[0], dict) and "type" in v[0]:
      data[k] = [_content_block_from_dict(inner) for inner in v]
    else:
      data[k] = v
  return ContentBlock(type=d["type"], data=data)


def _token_usage_from_dict(d):
  """Build a stored ``TokenUsage`` from its persisted dict.

  Reads every field declared on ``TokenUsage`` (the canonical field list),
  so a usage field added there is loaded automatically; absent keys fall
  back to the dataclass defaults. The single read-side counterpart of
  ``TokenUsage`` event ``to_stored()``.
  """
  from dataclasses import fields
  d = d or {}
  return TokenUsage(**{f.name: d[f.name]
                       for f in fields(TokenUsage) if f.name in d})


def persistable_prefix(messages):
  """Return ``messages`` truncated so it never ends in an un-persistable assistant
  tail: an assistant ``tool_use`` with no following ``tool_result``, OR an
  empty-content assistant message.

  Both are transient mid-turn states, never committed ones, and both are a
  non-recoverable provider 400 on the API backends if a crash makes them the
  newest on-disk state:

  - an unanswered ``tool_use`` -- the in-progress assistant right after it
    streamed a ``tool_use`` but before the session dispatched it (the
    orphaned-tool_use hazard ``AgentSession``'s cancel/turn-cap paths guard
    against). A completed turn never legitimately ends here (``run_turn`` always
    appends the answering ``tool_result``).
  - an EMPTY-content assistant -- the placeholder ``run_turn`` appends before it
    yields an error-before-content (an auth give-up, a rate-limit/network drop
    before any delta). Backends reject an empty-content assistant message, so
    persisting it bricks the conversation; a real assistant turn always carries
    content.

  Applied here, in the single choke point ``save`` uses, so EVERY writer -- the
  worker's mid-turn autosave and the GUI-thread turn-end / error / close / rename
  saves alike -- is safe. Returns the input list unchanged when nothing needs
  trimming.
  """
  msgs = messages
  while msgs and getattr(msgs[-1], "role", None) == "assistant" \
      and (not msgs[-1].content
           or any(b.type == "tool_use" for b in msgs[-1].content)):
    msgs = msgs[:-1]
  return msgs


def _message_to_dict(m):
  out = {
    "role": m.role,
    "timestamp": m.timestamp,
    "content": [_content_block_to_dict(b) for b in m.content],
  }
  if m.stop_reason is not None:
    out["stop_reason"] = m.stop_reason
  if m.usage is not None:
    out["usage"] = asdict(m.usage)
  if m.model is not None:
    out["model"] = m.model
  if m.backend is not None:
    out["backend"] = m.backend
  return out


def _message_from_dict(d):
  usage = None
  if "usage" in d and d["usage"] is not None:
    usage = _token_usage_from_dict(d["usage"])
  return Message(
    role=d["role"],
    content=[_content_block_from_dict(b) for b in d.get("content", [])],
    timestamp=_parse_dt(d["timestamp"]),
    stop_reason=d.get("stop_reason"),
    usage=usage,
    model=d.get("model"),
    backend=d.get("backend"),
  )


def _subagent_to_dict(r):
  return {
    "schema_version": _SCHEMA_VERSION,
    "sub_id": r.sub_id,
    "parent_conversation_id": r.parent_conversation_id,
    "parent_tool_use_id": r.parent_tool_use_id,
    "task": r.task,
    "profile_name": r.profile_name,
    "model": r.model,
    "started_at": r.started_at,
    "finished_at": r.finished_at,
    "final_text": r.final_text,
    "token_usage": asdict(r.token_usage),
    "messages": [_message_to_dict(m) for m in r.messages],
  }


def _subagent_from_dict(d):
  return SubagentRecord(
    sub_id=d["sub_id"],
    parent_conversation_id=d.get("parent_conversation_id", ""),
    parent_tool_use_id=d.get("parent_tool_use_id", ""),
    task=d.get("task", ""),
    profile_name=d.get("profile_name", ""),
    model=d.get("model", ""),
    started_at=_parse_dt(d["started_at"]),
    finished_at=_parse_dt(d["finished_at"]),
    final_text=d.get("final_text", ""),
    token_usage=_token_usage_from_dict(d.get("token_usage", {})),
    messages=[_message_from_dict(m) for m in d.get("messages", [])],
  )
