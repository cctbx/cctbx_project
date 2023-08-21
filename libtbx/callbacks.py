"""
Using GUI callbacks in CCTBX
----------------------------

The libtbx.callbacks module provides (or intends to) a transparent way of
sending information (in the form of pickled Python objects) from a running
non-graphical process to a graphical interface.  The intent is that any
program in CCTBX and derived works can import this module and call methods in
the global manager instance (libtbx.call_back), without worrying whether or
how the data will be sent to the parent process (if any - when the program is
being run by itself independently of a GUI or another process, the callbacks
will be ignored).  An example:

import libtbx.callbacks
libtbx.call_back(message="r_free",
  data=fmodel.r_free())
if (fmodel.r_free() > r_free_start):
  libtbx.warn("R-free increased during refinement")

For the callbacks, the message argument should be a short string (possibly
an actual message, or a statistic name, or a GUI method to call, or some kind
ofID), and the data argument can be any pickle-able object from a number to a
custom class.  There is no limit on the size of this object, but in practice it
will be faster not to pickle the contents of entire files if the filesystem is
being used for communication.  Since some callback data may be important as
part of a sequence, or needs to be shown to the user regardless of
communication latency, the accumulate and cached arguments toggle how
individual callbacks are handled.  The utility warn method sends a specific
type of message, a generic warning to the user, which in Phenix will be
prominently displayed in various windows.

The actual communication is handled either via connection objects from the
multiprocessing module (implemented in libtbx.thread_utils) or intermediate
files (Python pickles, as implemented in libtbx.runtime_utils).  These
modules handle the instantiation of processes (by whatever method) and the
interception of all output and status information, and its transmission to the
parent process.  (Note that standard output and exceptions are propagated
automatically, and do not require any additional callbacks.)  The GUI will
define its own methods to handle the various types of information,
which are called as necessary during runtime.  Ultimately the callbacks are
sent to the GUI as a specific type of event, which can be handled by custom
methods or ignored.  The module wxtbx.process_control provides an example of
how this could be done in wxPython (and is used in the Phenix GUI).

In theory it would be possible to implement a network-transparent callback
system , but we have not attempted this due to the difficulty of making it
work across a wide range of installations, varying levels of network security,
etc.

In practice, the callbacks used in Phenix include refinement statistics,
entire models or maps, plot data, progress bar control, and calls to window
methods.  They may be relatively frequent (as implied by the use of process
bars), but should not be so frequent or so large as to overwhelm the windowing
system (or filesystem, if intermediate files are used).  Standard output is
again handled automatically and buffered to limit the frequency.  Parallel
code should never use callbacks, as the behavior is undefined (but likely to
break).

A final note: while this is of course intended to be called from Python code,
it is possible to send callbacks from C++ code as well, and this is
implemented in the Phaser software.
"""
from __future__ import absolute_import, division, print_function

import six

from libtbx import str_utils
from libtbx import group_args, adopt_init_args
from libtbx.utils import to_unicode
import os
import sys

class manager(object):
  def __init__(self):
    self.handlers = []
    self._pid = os.getpid()
    self._log = None

  def set_warning_log(self, log):
    self._log = log

  def add_piped_callback(self, connection):
    cb = piped_callback(connection)
    self.register_handler(cb)

  def register_handler(self, handler):
    assert hasattr(handler, "__call__")
    self.handlers.append(handler)

  def __call__(self, message, data, accumulate=True, cached=True):
    for handler in self.handlers :
      handler(message=message,
              data=data,
              accumulate=accumulate,
              cached=cached)

  def add_citation(self, citation_info):
    self.__call__("citation", citation_info, accumulate=True)

  def showwarning(self, message, category, filename, lineno,
      *args, **kwds):
    # XXX ADDED *args because some calls come in here incorrectly with extra
    #   arguments, this way we just throw those away as this is just a warning
    self.warn(message)

  def warn(self, message):
    if (not isinstance(message, (six.text_type, six.binary_type))):
      message = to_unicode(message)
    self.__call__(message="warn", data=message, accumulate=True, cached=True)
    if (self._log is not None):
      log = self._log
    else :
      log = sys.stdout
    msg = "WARNING: %s\n" % message
    print("", file=log)
    for line in str_utils.line_breaker(msg, 72):
      print("  " + line, file=log)
    #print >> log, str_utils.wordwrap(msg, 79)

class piped_callback(object):
  def __init__(self, connection):
    adopt_init_args(self, locals())

  def __call__(self, message, data, accumulate=False, cached=True):
    self.connection.send(
      group_args(message=message,
                 data=data,
                 accumulate=accumulate,
                 cached=cached))

import libtbx
if not hasattr(libtbx, "call_back"):
  libtbx.call_back = manager()
  libtbx.warn = libtbx.call_back.warn
