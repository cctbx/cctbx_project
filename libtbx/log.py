"""A simple logging interface."""
from __future__ import absolute_import, division, print_function
import logging
import sys

class manager(object):
  def __init__(self, log=None):
    self.log      = log
    self.messages = []
    self.prefix   = ""

  def set_prefix(self, prefix):
    self.prefix=prefix

  def add(self, msg):
    msg = "%s%s"%(self.prefix, msg)
    self.messages.append(msg)

  def add_and_show(self, msg):
    self.add(msg)
    self.show_last()

  def show(self):
    print(msg, file=self.log)

  def show_last(self):
    if len(self.messages)>0:
      print(self.messages[-1], file=self.log)

  def show_all(self):
    for msg in self.messages:
      print(msg, file=self.log)


class logger(object):
  """A basic wrapper over Python standard library logging module.

  This class is designed to be used as a logging singleton, with a few
  basic convenience methods:
    set_quiet       Suppress output except errors.
    set_logfile     Set an output log file.
    set_stdout      Set log output to stdout.

  It may also be written to as a file-like object.

  Examples:

  # Basic usage
  >>> import libtbx.log
  >>> libtbx.log.info("Helpful message")

  # Write as a file handle
  >>> print >> libtbx.log.logger, "As file handle"

  # Set debug, or quiet states
  >>> libtbx.log.logger.set_debug(True)
  >>> libtbx.log.debug("Debugging output.")
  >>> libtbx.log.logger.set_quiet(True)
  >>> libtbx.log.info("This should be muted.")

  # Output file
  >>> libtbx.log.logger.set_logfile("test.log")
  >>> libtbx.log.info("Redirected to test.log")

  # Experimental sys.stdout redirection
  >>> libtbx.log.logger.set_logfile("test.log")
  >>> sys.stdout = libtbx.log.logger
  >>> print "stdout redirection to test.log"

  """

  def __init__(self):
    # create logger
    self.log = logging.getLogger('phenix')
    self.log.setLevel(logging.INFO)

    # create formatter
    self.formatter = logging.Formatter()

    # Default is to print to stdout
    self.ch = None
    self.set_stdout()

  def write(self, data):
    self.log.info(data)

  def flush(self):
    try:
      self.ch.flush()
    except Exception:
      pass

  def close(self):
    try:
      self.ch.close()
    except Exception:
      pass

  def set_quiet(self, state=True):
    level = logging.INFO
    if state:
      level = logging.ERROR
    self.log.setLevel(level)

  def set_debug(self, state=True):
    level = logging.INFO
    if state:
      level = logging.DEBUG
    self.log.setLevel(level)

  def set_logfile(self, filename, mode='w'):
    ch = logging.FileHandler(filename, mode=mode)
    ch.setFormatter(self.formatter)
    if self.ch:
      self.log.removeHandler(self.ch)
    self.log.addHandler(ch)
    self.ch = ch

  def set_fileobj(self, fileobj):
    pass

  def set_stdout(self):
    ch = logging.StreamHandler()
    ch.setFormatter(self.formatter)
    if self.ch:
      self.log.removeHandler(self.ch)
    self.log.addHandler(ch)
    self.ch = ch

# Singleton
logger = logger()

# Set the module functions.
log = logger.log.info
debug = logger.log.debug
info = logger.log.info
warn = logger.log.warn
error = logger.log.error
critical = logger.log.critical

##### Testing #####

import unittest
class TestLog(unittest.TestCase):
  def test_f(self):
    logger.set_debug(True)
    tests = [debug, info, warn, error, critical]
    for test in tests:
      logger.set_logfile("test.log")
      value = "test %s"%test
      test(value)
      logger.flush()
      logger.close()
      with open("test.log") as f:
        data = f.read()
        assert value in data

  def test_debug(self):
    print("Check set_debug")
    logger.set_logfile("test.log")
    logger.set_debug(False)
    debug("debug: muted")
    logger.set_debug(True)
    debug("debug: ok")
    logger.set_debug(False)
    with open("test.log") as f:
      data = f.read()
      assert "muted" not in data
      assert "ok" in data

  def test_quiet(self):
    print("Check set_quiet")
    logger.set_logfile("test.log")
    logger.set_quiet(False)
    info("quiet: ok")
    logger.set_quiet(True)
    info("quiet: muted")
    logger.set_quiet(False)
    with open("test.log") as f:
      data = f.read()
      assert "muted" not in data
      assert "ok" in data

  def test_logfile(self):
    print("Check output file")
    logger.set_logfile("test.log")
    info("logfile: ok")
    with open("test.log") as f:
      data = f.read()
      assert "ok" in data

  def test_stdout(self):
    oldsys = sys.stdout
    print("Checking sys.stdout redirect")
    sys.stdout = logger
    print("redirect: ok")
    sys.stdout = oldsys
    print("... and back again")

if __name__ == "__main__":
  unittest.main(verbosity=0)
