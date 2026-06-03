"""Path resolution for project directory and chat-data root.

The active project directory is where conversations and Phenix job dirs
live; ``.phenix_chat/`` is created lazily on first persistence write so
cwd-as-project doesn't leave a hidden dir behind for chats that never
sent a message.
"""

import os
from pathlib import Path


def resolve_project_dir(cli_arg=None, embedded_arg=None):
  """Resolve the active project directory.

  Resolution order: the CLI ``--project-dir`` argument, then
  ``embedded_arg``, then the ``PHENIX_PROJECT_DIR`` env var, then the
  current working directory (default — cwd-as-project).

  Parameters
  ----------
  cli_arg : str, optional
      The CLI ``--project-dir`` argument passed by the launcher.
  embedded_arg : str, optional
      A directory passed by an embedding GUI's constructor.

  Returns
  -------
  pathlib.Path
      The resolved (absolute) active project directory.
  """
  if cli_arg is not None:
    return Path(cli_arg).resolve()
  if embedded_arg is not None:
    return Path(embedded_arg).resolve()
  env = os.environ.get("PHENIX_PROJECT_DIR")
  if env:
    return Path(env).resolve()
  return Path.cwd().resolve()


def chat_root_for(project_dir):
  """Resolve the chat-data root inside a project directory.

  Resolution order: the ``PHENIX_CHAT_HOME`` env var (test/debug
  override), then ``<project_dir>/.phenix_chat/``.

  The directory is NOT created here. Storage code creates it lazily on the
  first persistence write so cwd-as-project doesn't leave a hidden dir
  behind for chats that never sent a message.

  Parameters
  ----------
  project_dir : str or pathlib.Path
      The active project directory.

  Returns
  -------
  pathlib.Path
      The chat-data root (which may not yet exist on disk).
  """
  override = os.environ.get("PHENIX_CHAT_HOME")
  if override:
    return Path(override)
  return Path(project_dir) / ".phenix_chat"
