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
  """Marker base for events streamed from Agent.stream_turn."""


@dataclass
class AgentError(AgentEvent):
  """Yielded by Agent.stream_turn when an API / network / auth error occurs.
  recoverable=True means the UI can offer a retry button (rate limit, 5xx,
  network); recoverable=False is a hard fail (auth, bad request).

  kind is an optional short classifier (e.g. 'auth', 'rate_limit', 'network',
  'internal') that the GUI can use to pick an icon or copy."""
  message: str = ""
  recoverable: bool = True
  kind: str = None


class CancelToken:
  """Thread-safe boolean flag. Set by the GUI's Stop handler; checked
  by the worker thread between events and tool dispatches. Setting it
  does NOT wake a thread parked on a ``queue.Queue`` — that needs the
  ``_Cancelled`` sentinel push (see ``tools.py``)."""

  def __init__(self):
    self._event = threading.Event()

  def set(self):
    self._event.set()

  def is_set(self):
    return self._event.is_set()


class TurnCancelled(Exception):
  """Raised by ``AgentSession._await_approval`` when the approval queue
  delivers a ``_Cancelled`` sentinel. Caught in
  ``_dispatch_and_build_results`` to end the turn cleanly."""
