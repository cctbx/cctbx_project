"""Chat session logging with secret redaction.

One log file per launcher invocation, plus a redaction helper applied to
any payload-bearing log call.

Rotation policy: keep up to 7 most-recent session logs under
``<chat_root>/logs/``. Older logs are deleted on launch.
"""

import datetime
import re

LOG_KEEP = 7

# Regexes are anchored loosely so the helper catches the common header
# spellings without being clever about quoting.
_REDACT_PATTERNS = (
  (re.compile(r'(sk-ant-[A-Za-z0-9_\-]+)'), '<REDACTED>'),
  (re.compile(r'(Bearer\s+)[A-Za-z0-9_\-\.]+', re.IGNORECASE), r'\1<REDACTED>'),
  (re.compile(r'(x-api-key\s*[:=]\s*)\S+', re.IGNORECASE), r'\1<REDACTED>'),
  (re.compile(r'("api_key"\s*:\s*)"[^"]+"'), r'\1"<REDACTED>"'),
)


def session_log_path(chat_root):
  """Return a timestamped session log path under ``chat_root/logs``."""
  ts = datetime.datetime.now().strftime("%Y%m%dT%H%M%S%f")
  return chat_root / "logs" / ("chat-%s.log" % ts)


def open_session_log(chat_root):
  """Open a redacting session log, creating the dir and pruning old logs.

  The caller owns the returned handle and must close it.

  Parameters
  ----------
  chat_root : pathlib.Path
      Chat root directory; logs live under ``chat_root/logs``.

  Returns
  -------
  tuple of (_RedactingLog, pathlib.Path)
      The open, line-buffered, redacting log handle and its path.
  """
  log_dir = chat_root / "logs"
  log_dir.mkdir(parents=True, exist_ok=True)
  # Prune before creating today's log so LOG_KEEP counts only pre-existing
  # files; the new log we're about to open is excluded from the count.
  _prune_old_logs(log_dir, "chat-*.log")
  path = session_log_path(chat_root)
  return _RedactingLog(open(path, "w", buffering=1)), path


def debug_log_path(chat_root):
  """Return a timestamped debug log path under ``chat_root/logs``."""
  ts = datetime.datetime.now().strftime("%Y%m%dT%H%M%S%f")
  return chat_root / "logs" / ("debug-%s.log" % ts)


def open_debug_log(chat_root):
  """Open a per-launch debug log alongside the session log.

  Distinct from :func:`open_session_log`: only created when the launcher
  is invoked with ``--debug``. Collects what would otherwise be lost
  when ``phenix.chat`` runs without a visible terminal -- uncaught
  exceptions, ``AgentError`` events surfaced by the runner, and (for the
  claude_code backend) the SDK's ``debug_stderr`` stream. Same
  ``LOG_KEEP`` rotation as session logs but tracked separately so debug
  history doesn't push out routine session logs.

  The caller owns the returned handle and must close it.

  Parameters
  ----------
  chat_root : pathlib.Path
      Chat root directory; logs live under ``chat_root/logs``.

  Returns
  -------
  tuple of (_RedactingLog, pathlib.Path)
      The open, line-buffered, redacting log handle and its path.
  """
  log_dir = chat_root / "logs"
  log_dir.mkdir(parents=True, exist_ok=True)
  _prune_old_logs(log_dir, "debug-*.log")
  path = debug_log_path(chat_root)
  return _RedactingLog(open(path, "w", buffering=1)), path


def redact_secrets(text):
  """Apply all redaction patterns to ``text``. Idempotent.

  Parameters
  ----------
  text : str
      Text that may contain secrets (API keys, bearer tokens, etc.).

  Returns
  -------
  str
      ``text`` with any matched secrets replaced by ``<REDACTED>``.
  """
  for pat, repl in _REDACT_PATTERNS:
    text = pat.sub(repl, text)
  return text


class _RedactingLog:
  """File-like wrapper that runs :func:`redact_secrets` over every write.

  This is the seam that makes the redaction helper actually do something:
  the session and debug logs collect tracebacks, surfaced ``AgentError``
  messages, tool output and (for backends that route it here) SDK debug
  streams, any of which could carry an API key / OAuth token. Wrapping the
  handle scrubs every payload-bearing write regardless of which caller
  produced it.

  Redaction is applied per :meth:`write` call, so a secret split across
  two separate ``write()`` calls is not matched -- callers log whole
  messages, so this is not a concern in practice.

  ``fileno`` deliberately raises (see :meth:`fileno`): the wrapper only
  redacts Python-level ``write()`` calls, so it must not be handed to a
  subprocess as a stdout/stderr fd -- the child would write to the raw
  descriptor and bypass redaction. Other non-write attributes (``flush``
  / ``close`` / ...) delegate to the wrapped handle.

  Parameters
  ----------
  fh : file-like
      The underlying writable text handle to forward (redacted) writes to.
      Ownership transfers to the wrapper: closing the wrapper closes ``fh``.
  """

  def __init__(self, fh):
    self._fh = fh

  def write(self, text):
    """Redact secrets in ``text`` and forward it to the wrapped handle."""
    return self._fh.write(redact_secrets(text))

  def writelines(self, lines):
    """Write each line through :meth:`write` so every line is redacted."""
    for line in lines:
      self.write(line)

  def __getattr__(self, name):
    """Delegate any non-overridden attribute to the wrapped handle."""
    return getattr(self._fh, name)

  def __enter__(self):
    """Return self so the wrapper can be used as a context manager."""
    return self

  def __exit__(self, *exc):
    """Close the wrapped handle on context-manager exit."""
    self._fh.close()
    return False

  def fileno(self):
    """Refuse to expose a raw file descriptor.

    The wrapper redacts only Python-level ``write()`` calls, so handing it
    to a subprocess as a stdout/stderr fd would let the child write to the
    raw descriptor and bypass redaction. Raising keeps that misuse from
    happening silently.

    Raises
    ------
    io.UnsupportedOperation
        Always; route output through :meth:`write` instead.
    """
    import io
    raise io.UnsupportedOperation(
      "_RedactingLog exposes no fileno: a raw fd would bypass secret "
      "redaction; write() through the wrapper instead")


def _prune_old_logs(log_dir, pattern):
  logs = sorted(log_dir.glob(pattern),
                key=lambda p: p.stat().st_mtime, reverse=True)
  for p in logs[LOG_KEEP:]:
    try:
      p.unlink()
    except OSError:
      pass
