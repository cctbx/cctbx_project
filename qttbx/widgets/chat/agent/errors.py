"""Error and cancellation types for the chat session.

``CancelToken`` is a thread-safe flag the GUI sets when the user clicks
Stop. The ``AgentSession`` checks it between events, API calls, and
tool dispatches. ``TurnCancelled`` is raised when ``_await_approval``
is woken by a ``_Cancelled`` sentinel pushed into the approval queue
by the GUI on Stop.
"""

import threading
from dataclasses import dataclass


@dataclass
class AgentEvent:
  """Marker base for events streamed from ``Agent.stream_turn``."""


@dataclass
class AgentError(AgentEvent):
  """Yielded by ``Agent.stream_turn`` on an API / network / auth error.

  ``recoverable=True`` means the UI can offer a retry button (rate limit,
  5xx, network); ``recoverable=False`` is a hard fail (auth, bad request).

  ``kind`` is an optional short classifier (e.g. ``'auth'``,
  ``'rate_limit'``, ``'network'``, ``'internal'``) that the GUI can use to
  pick an icon or copy.
  """
  message: str = ""
  recoverable: bool = True
  kind: str = None


def map_httpx_sdk_error(sdk_module, exc):
  """Map an httpx-based SDK exception (openai or anthropic) onto an AgentError.

  The openai and anthropic SDKs share an exception hierarchy in which
  ``AuthenticationError`` / ``RateLimitError`` / ``APIConnectionError`` /
  ``APIStatusError`` are all subclasses of ``APIError``. Both chat backends map
  that hierarchy onto ``AgentError`` with one policy, parameterized here by the
  SDK module so a single implementation serves both:

    AuthenticationError -> recoverable=False, kind="auth"
    RateLimitError      -> recoverable=True,  kind="rate_limited"
    APIConnectionError  -> recoverable=True   ("Network error: ...")
    APIStatusError      -> recoverable=(status_code >= 500) ("API <code>: ...")
    APIError (catch-all)-> recoverable=False  ("API error: ...")

  Anything that is not one of these five SDK types is re-raised, so a non-SDK
  bug is never silently turned into an AgentError.

  Parameters
  ----------
  sdk_module : module
      The provider SDK (``openai`` or ``anthropic``) supplying the exception
      classes the isinstance ladder tests against.
  exc : Exception
      The exception raised by the SDK call.

  Returns
  -------
  AgentError
      The mapped event.

  Raises
  ------
  Exception
      ``exc`` itself, when it matches none of the five SDK types.
  """
  if isinstance(exc, sdk_module.AuthenticationError):
    return AgentError(message=str(exc), recoverable=False, kind="auth")
  if isinstance(exc, sdk_module.RateLimitError):
    return AgentError(message=str(exc), recoverable=True, kind="rate_limited")
  if isinstance(exc, sdk_module.APIConnectionError):
    return AgentError(message="Network error: %s" % exc, recoverable=True)
  if isinstance(exc, sdk_module.APIStatusError):
    code = getattr(exc, "status_code", 0)
    return AgentError(message="API %s: %s" % (code, exc),
                      recoverable=(code >= 500))
  if isinstance(exc, sdk_module.APIError):
    return AgentError(message="API error: %s" % exc, recoverable=False)
  raise exc


class CancelToken:
  """Thread-safe boolean flag.

  Set by the GUI's Stop handler; checked by the worker thread between
  events and tool dispatches. Setting it does NOT wake a thread parked on
  a ``queue.Queue`` — that needs the ``_Cancelled`` sentinel push (see
  ``tools.py``).
  """

  def __init__(self):
    self._event = threading.Event()

  def set(self):
    self._event.set()

  def is_set(self):
    return self._event.is_set()


class TurnCancelled(Exception):
  """Raised when the approval queue delivers a ``_Cancelled`` sentinel.

  Raised by ``AgentSession._await_approval`` and caught in
  ``_dispatch_and_build_results`` to end the turn cleanly.
  """
