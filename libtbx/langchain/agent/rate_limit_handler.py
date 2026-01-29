"""
Rate limit handling with exponential backoff and decay.

This module provides utilities for handling API rate limits with:
- Exponential backoff on rate limit errors
- Automatic decay of backoff after successful calls
- Configurable retry limits and delays

Usage:
    from libtbx.langchain.agent.rate_limit_handler import RateLimitHandler, call_with_retry

    # Simple usage with default handler
    result = call_with_retry(lambda: api_call(), log_func)

    # Or create a named handler for specific APIs
    handler = RateLimitHandler.get_handler("google_api")
    result = handler.call_with_retry(lambda: api_call(), log_func)
"""

from __future__ import absolute_import, division, print_function

import time
import threading


class RateLimitHandler:
    """
    Handles rate limiting with exponential backoff and decay.

    The handler tracks rate limit state and adjusts delays accordingly:
    - On rate limit: exponential backoff (2s -> 4s -> 8s -> ...)
    - On success: decay back to base delay after decay_time seconds

    Handlers are cached by name so the same handler is reused across calls.
    """

    _handlers = {}
    _lock = threading.Lock()

    def __init__(self, name="default", max_retries=3, base_delay=2.0,
                 max_delay=60.0, decay_time=300.0):
        """
        Initialize rate limit handler.

        Args:
            name: Handler name for identification
            max_retries: Maximum retry attempts (default: 3)
            base_delay: Initial delay in seconds (default: 2.0)
            max_delay: Maximum delay cap in seconds (default: 60.0)
            decay_time: Seconds after success to reset delay (default: 300 = 5 min)
        """
        self.name = name
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.decay_time = decay_time

        # Current state
        self._current_delay = base_delay
        self._last_success_time = None
        self._consecutive_rate_limits = 0
        self._state_lock = threading.Lock()

    @classmethod
    def get_handler(cls, name="default", **kwargs):
        """
        Get or create a named handler.

        Handlers are cached by name so rate limit state persists across calls.

        Args:
            name: Handler name
            **kwargs: Passed to __init__ if creating new handler

        Returns:
            RateLimitHandler instance
        """
        with cls._lock:
            if name not in cls._handlers:
                cls._handlers[name] = cls(name=name, **kwargs)
            return cls._handlers[name]

    @classmethod
    def reset_handler(cls, name="default"):
        """Reset a handler's state."""
        with cls._lock:
            if name in cls._handlers:
                handler = cls._handlers[name]
                with handler._state_lock:
                    handler._current_delay = handler.base_delay
                    handler._last_success_time = None
                    handler._consecutive_rate_limits = 0

    @classmethod
    def reset_all_handlers(cls):
        """Reset all handlers' state."""
        with cls._lock:
            for name in list(cls._handlers.keys()):
                cls.reset_handler(name)

    def _check_decay(self):
        """Check if enough time has passed to decay the delay."""
        with self._state_lock:
            if self._last_success_time is not None:
                elapsed = time.time() - self._last_success_time
                if elapsed >= self.decay_time:
                    # Reset to base delay after decay period
                    self._current_delay = self.base_delay
                    self._consecutive_rate_limits = 0

    def _record_success(self):
        """Record a successful call."""
        with self._state_lock:
            self._last_success_time = time.time()
            # Gradually reduce delay on success (but not immediately to base)
            if self._consecutive_rate_limits > 0:
                self._consecutive_rate_limits = max(0, self._consecutive_rate_limits - 1)
                # Reduce delay by half, but not below base
                self._current_delay = max(self.base_delay, self._current_delay / 2)

    def _record_rate_limit(self):
        """Record a rate limit hit and return the delay to use."""
        with self._state_lock:
            self._consecutive_rate_limits += 1
            # Exponential backoff
            self._current_delay = min(
                self.max_delay,
                self.base_delay * (2 ** self._consecutive_rate_limits)
            )
            return self._current_delay

    def get_current_delay(self):
        """Get the current delay (for informational purposes)."""
        self._check_decay()
        with self._state_lock:
            return self._current_delay

    def call_with_retry(self, call_func, log_func=None):
        """
        Call a function with retry on rate limit errors.

        Args:
            call_func: Function to call (should raise exception on rate limit)
            log_func: Optional logging function (takes string message)

        Returns:
            Result from call_func

        Raises:
            Exception: If max retries exceeded or non-rate-limit error
        """
        if log_func is None:
            log_func = lambda msg: None

        # Check if we should decay the delay
        self._check_decay()

        last_error = None

        for attempt in range(self.max_retries + 1):
            try:
                result = call_func()
                self._record_success()
                return result

            except Exception as e:
                last_error = e

                if is_rate_limit_error(e):
                    if attempt < self.max_retries:
                        delay = self._record_rate_limit()
                        log_func("Rate limit hit (attempt %d/%d), waiting %.1fs before retry" % (
                            attempt + 1, self.max_retries, delay))
                        time.sleep(delay)
                        continue
                    else:
                        log_func("Rate limit - max retries (%d) exceeded" % self.max_retries)
                        raise
                else:
                    # Non-rate-limit error, don't retry
                    raise

        # Should not reach here, but just in case
        if last_error:
            raise last_error
        return None


def is_rate_limit_error(error):
    """
    Check if an exception is a rate limit error.

    Args:
        error: Exception to check

    Returns:
        bool: True if this appears to be a rate limit error
    """
    error_str = str(error).lower()
    error_type = type(error).__name__.lower()

    # Check for common rate limit indicators
    rate_limit_indicators = [
        "429",
        "503",  # Service unavailable / overloaded
        "rate limit",
        "rate_limit",
        "ratelimit",
        "resource exhausted",
        "resourceexhausted",
        "quota exceeded",
        "quotaexceeded",
        "too many requests",
        "toomanyrequests",
        "throttl",
        "retry after",
        "retry-after",
        "slow down",
        "capacity",
        "overloaded",
        "unavailable",  # Model unavailable / overloaded
    ]

    for indicator in rate_limit_indicators:
        if indicator in error_str or indicator in error_type:
            return True

    # Check for specific exception types
    rate_limit_exceptions = [
        "ratelimiterror",
        "ratelimitexceeded",
        "toomanyrequestserror",
        "resourceexhausted",
        "quotaexceedederror",
    ]

    for exc_type in rate_limit_exceptions:
        if exc_type in error_type:
            return True

    return False


def call_with_retry(call_func, log_func=None, handler_name="default", **handler_kwargs):
    """
    Convenience function to call with retry using a named handler.

    Args:
        call_func: Function to call
        log_func: Optional logging function
        handler_name: Name of handler to use (handlers are cached)
        **handler_kwargs: Passed to handler creation if new

    Returns:
        Result from call_func
    """
    handler = RateLimitHandler.get_handler(handler_name, **handler_kwargs)
    return handler.call_with_retry(call_func, log_func)


# Pre-configured handlers for common APIs
def get_google_handler():
    """Get handler configured for Google APIs."""
    return RateLimitHandler.get_handler(
        "google_api",
        max_retries=5,
        base_delay=2.0,
        max_delay=120.0,  # Up to 2 minutes
        decay_time=300.0  # 5 minutes
    )


def get_openai_handler():
    """Get handler configured for OpenAI APIs."""
    return RateLimitHandler.get_handler(
        "openai_api",
        max_retries=5,
        base_delay=1.0,
        max_delay=120.0,  # Up to 2 minutes
        decay_time=300.0
    )


def get_anthropic_handler():
    """Get handler configured for Anthropic APIs."""
    return RateLimitHandler.get_handler(
        "anthropic_api",
        max_retries=5,
        base_delay=1.0,
        max_delay=120.0,  # Up to 2 minutes
        decay_time=300.0
    )
