from __future__ import absolute_import, division, print_function

import os
import subprocess
import sys
import threading
import signal
from six.moves import range


def _show_lines(lines, out, prefix):
  if (out is None): out = sys.stdout
  for line in lines:
    print(prefix+line, file=out)

def macos_dyld():
  '''
  Convenience function for returning either DYLD_LIBRARY_PATH or
  DYLD_FALLBACK_LIBRARY_PATH (for conda environments)
  '''
  dyld_options = ['DYLD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH']
  for dyld in dyld_options:
    dyld_path = os.environ.get(dyld)
    if (dyld_path is not None):
      return '%s="%s"' % (dyld, dyld_path)
  return 'DYLD_LIBRARY_PATH= '

class fully_buffered_base(object):

  def format_errors_if_any(self):
    assert not self.join_stdout_stderr
    if (len(self.stderr_lines) != 0):
      msg = ["child process stderr output:"]
      msg.append("  command: " + repr(self.command))
      for line in self.stderr_lines:
        msg.append("  " + line)
      return "\n".join(msg)
    if (self.return_code != 0):
      return "non-zero return code: %s"%(self.return_code)
    return None

  def raise_if_errors(self, Error=RuntimeError):
    assert not self.join_stdout_stderr
    msg = self.format_errors_if_any()
    if (msg is not None):
      raise Error(msg)
    return self

  def raise_if_output(self, show_output_threshold=10, Error=RuntimeError):
    def start_msg():
      result = ["unexpected child process output:"]
      result.append("  command: " + repr(self.command))
      return result
    if (self.stdout_buffer is not None):
      if (len(self.stdout_buffer) != 0):
        msg = start_msg()
        msg.append("  length of output: %d bytes" % len(self.stdout_buffer))
        raise Error("\n".join(msg))
    elif (len(self.stdout_lines) != 0):
      msg = start_msg()
      for line in self.stdout_lines[:show_output_threshold]:
        msg.append("  " + line)
      n = len(self.stdout_lines)
      if (n > show_output_threshold):
        if (n <= show_output_threshold+2):
          for line in self.stdout_lines[show_output_threshold:n]:
            msg.append("  " + line)
        else:
          msg.append("  ...")
          msg.append("  remaining %d lines omitted."
            % (n-show_output_threshold))
      raise Error("\n".join(msg))
    return self

  def raise_if_errors_or_output(self, Error=RuntimeError):
    self.raise_if_errors(Error=Error)
    self.raise_if_output(Error=Error)
    return self

  def show_stderr(self, out=None, prefix=""):
    _show_lines(lines=self.stderr_lines, out=out, prefix=prefix)

  def show_stdout(self, out=None, prefix=""):
    assert self.stdout_lines is not None
    _show_lines(lines=self.stdout_lines, out=out, prefix=prefix)

class fully_buffered_simple(fully_buffered_base):
  """\
Executes command, sends stdin_lines (str or sequence), then reads
stdout_lines first, stderr_lines second (if join_stdout_stderr
is False).

The constructor may deadlock if the I/O buffers are too small to allow
the blocking write and reads in the given sequence. Specifically,
stdin_lines may be too big, or there may be too many stderr_lines,
but there can be any number of stdout_lines. The tests below are
known to work under Mac OS X, Windows XP, IRIX, and Tru64 Unix with
stdin_lines up to 1000000, stderr_lines up to 500. I.e. this simple
implementation should cover most practical situations.
  """

  def __init__(self,
        command,
        stdin_lines=None,
        join_stdout_stderr=False,
        stdout_splitlines=True,
        bufsize=-1):
    self.command = command
    self.join_stdout_stderr = join_stdout_stderr
    if (join_stdout_stderr):
      child_stdin, child_stdout = os.popen4(command, "t", bufsize)
      child_stderr = None
    else:
      child_stdin, child_stdout, child_stderr = os.popen3(command,"t",bufsize)
    if (stdin_lines is not None):
      if (not isinstance(stdin_lines, str)):
        stdin_lines = '\n'.join(stdin_lines)
        if (len(stdin_lines) != 0):
          stdin_lines += '\n'
      child_stdin.write(stdin_lines)
    child_stdin.close()
    if (stdout_splitlines):
      self.stdout_buffer = None
      self.stdout_lines = child_stdout.read().splitlines()
    else:
      self.stdout_buffer = child_stdout.read()
      self.stdout_lines = None
    if (child_stderr is not None):
      self.stderr_lines = child_stderr.read().splitlines()
    else:
      self.stderr_lines = []
    child_stdout.close()
    if (child_stderr is not None):
      child_stderr.close()
    self.return_code = None

class fully_buffered_subprocess(fully_buffered_base):
  "This implementation is supposed to never block."

  def __init__(self,
        command,
        timeout=None,
        stdin_lines=None,
        join_stdout_stderr=False,
        stdout_splitlines=True,
        bufsize=-1):
    def target(process, lines, result):
      o, e = process.communicate(input=lines)
      result[0] = o
      result[1] = e

    self.command = command
    self.join_stdout_stderr = join_stdout_stderr
    if (not isinstance(command, str)):
      command = subprocess.list2cmdline(command)
    # Timeout functionality based on:
    # https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    # https://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true
    if (sys.platform == 'darwin'):   # bypass SIP on OS X 10.11
      command = ('%s exec ' % macos_dyld()) + command
    if (stdin_lines is not None):
      if (not isinstance(stdin_lines, str)):
        stdin_lines = '\n'.join(stdin_lines)
        if (len(stdin_lines) != 0):
          stdin_lines += '\n'
    if (join_stdout_stderr):
      stderr = subprocess.STDOUT
    else:
      stderr = subprocess.PIPE
    p = subprocess.Popen(
      args=command,
      shell=True,
      bufsize=bufsize,
      stdin=subprocess.PIPE,
      stdout=subprocess.PIPE,
      stderr=stderr,
      universal_newlines=True,
      close_fds=(sys.platform != 'win32'),
      preexec_fn=os.setsid if sys.platform != 'win32' else None)
    if timeout is not None:
      if sys.platform != 'win32':
        r = [None, None]
        thread = threading.Thread(target=target, args=(p, stdin_lines, r))
        thread.start()
        thread.join(timeout)
        if thread.is_alive():
          os.killpg(os.getpgid(p.pid), signal.SIGTERM)
          thread.join()
        o, e = r[0], r[1]
      else: # sys.platform == 'win32'
        # don't respect timeout for now
        o, e = p.communicate(input=stdin_lines)
    else:
      o, e = p.communicate(input=stdin_lines)
    if (stdout_splitlines):
      self.stdout_buffer = None
      self.stdout_lines = o.splitlines()
    else:
      self.stdout_buffer = o
      self.stdout_lines = None
    if (join_stdout_stderr):
      self.stderr_lines = []
    else:
      self.stderr_lines = e.splitlines()
    self.return_code = p.returncode

fully_buffered = fully_buffered_subprocess

def go(command, stdin_lines=None,join_stdout_stderr=True):
  return fully_buffered(
    command=command,
    stdin_lines=stdin_lines,
    join_stdout_stderr=join_stdout_stderr)

def call(command):
  """
  Wraps subprocess.call to run a command.

  Parameters
  ----------
  command : str

  Returns
  -------
  int
      Exit code of subprocess.

  Examples
  --------
  >>> from libtbx.easy_run import call
  >>> ret = call("echo 1")
  1
  >>> print ret
  0
  """
  for s in [sys.stdout, sys.stderr]:
    flush = getattr(s, "flush", None)
    if (flush is not None): flush()
  if (sys.platform == 'darwin'):   # bypass SIP on OS X 10.11
    command = ('%s exec ' % macos_dyld()) + command
  return subprocess.call(args=command, shell=True)

def exercise(args=None):
  from six.moves import cStringIO as StringIO
  if (args is None): args = sys.argv[1:]
  verbose = "--verbose" in args
  #
  if ("--simple" in args):
    fb = fully_buffered_simple
  else:
    fb = fully_buffered
  #
  for command in ["echo hello world", ("echo", "hello", "world")]:
    for result in [fb(command=command).raise_if_errors(),
                   fb(command=command, join_stdout_stderr=True),
                   go(command=command)]:
      if verbose: print(result.stdout_lines)
      assert result.stdout_lines == ["hello world"]
  #
  if (os.path.isfile("/bin/ls")):
    for command in ["/bin/ls /bin", ("/bin/ls", "/bin")]:
      result = fb(command=command).raise_if_errors()
      if verbose: print(result.stdout_lines)
      assert "ls" in result.stdout_lines
  if (os.path.isfile("/usr/bin/wc")):
    for command in ["/usr/bin/wc -l", ("/usr/bin/wc", "-l")]:
      result = fb(command=command).raise_if_errors()
      if verbose: print(result.stdout_lines)
      assert [s.strip() for s in result.stdout_lines] == ["0"]
      result = fb(command=command, stdin_lines=["hello"]) \
        .raise_if_errors()
      if verbose: print(result.stdout_lines)
      assert [s.strip() for s in result.stdout_lines] == ["1"]
      result = fb(command=command, stdin_lines=["hello", "world"]) \
        .raise_if_errors()
      if verbose: print(result.stdout_lines)
      assert [s.strip() for s in result.stdout_lines] == ["2"]
      result = fb(command=command, stdin_lines="hello\nworld\nbye\n") \
        .raise_if_errors()
      if verbose: print(result.stdout_lines)
      assert [s.strip() for s in result.stdout_lines] == ["3"]
  #
  if (os.name == "nt"):
    result = fb(command="dir").raise_if_errors()
    if verbose: print(result.stdout_lines)
    assert len(result.stdout_lines) > 0
    windir = os.environ.get("windir", None)
    if (windir is not None and windir.find(" ") < 0):
      result = fb(command="dir "+windir).raise_if_errors()
      if verbose: print(result.stdout_lines)
      assert len(result.stdout_lines) > 0
  #
  pyexe = sys.executable
  assert pyexe.count('"') == 0
  pyexe = '"' + pyexe + '"'
  if (os.name == "nt"):
    pyexe = "call " + pyexe
  #
  if ("PYTHONPATH" in os.environ):
    if (not hasattr(os, "unsetenv")):
      os.environ["PYTHONPATH"] = ""
    else:
      del os.environ["PYTHONPATH"]
  if (os.name == "nt"):
    result = fb(command="set").raise_if_errors()
  elif (os.path.isfile("/usr/bin/printenv")):
    result = fb(command="/usr/bin/printenv").raise_if_errors()
  else:
    result = None
  if (result is not None):
    if verbose: print(result.stdout_lines)
    for line in result.stdout_lines:
      assert not line.startswith("PYTHONPATH") or line == "PYTHONPATH="
  #
  for stdout_splitlines in [True, False]:
    result = fb(
      command="%s -V" % pyexe,
      stdout_splitlines=stdout_splitlines)
    # python -V outputs to stdout or stderr depending on version
    # https://bugs.python.org/issue18338
    if (len(result.stderr_lines) > 0):
      if verbose: print(result.stderr_lines)
      assert result.stderr_lines[0].startswith(
        "Python " + sys.version.split()[0])
      if (stdout_splitlines):
        assert result.stdout_buffer is None
        assert result.stdout_lines == []
      else:
        assert result.stdout_buffer == ""
        assert result.stdout_lines is None
    else:
      if verbose: print(result.stdout_lines)
      if (stdout_splitlines):
        assert result.stdout_buffer is None
        assert result.stdout_lines[0].startswith(
          "Python " + sys.version.split()[0])
      else:
        assert result.stdout_buffer.startswith(
          "Python " + sys.version.split()[0])
        assert result.stdout_lines is None
  result = go(command="%s -V" % pyexe)
  if verbose: print(result.stdout_lines)
  assert result.stdout_lines[0].startswith("Python " + sys.version.split()[0])
  result = fb(
    command='%s -c "print(3+4)"' % pyexe).raise_if_errors()
  if verbose: print(result.stdout_lines)
  assert result.stdout_lines == ["7"]
  command = command = pyexe \
    + ' -c "import sys; print(len(list(filter(bool, sys.stdin.read().splitlines()))))"'
  result = fb(command=command).raise_if_errors()
  if verbose: print(result.stdout_lines)
  assert result.stdout_lines == ["0"]
  result = fb(command=command, stdin_lines=["hello"]) \
    .raise_if_errors()
  if verbose: print(result.stdout_lines)
  assert result.stdout_lines == ["1"]
  result = fb(command=command, stdin_lines=["hello", "world"]) \
    .raise_if_errors()
  if verbose: print(result.stdout_lines)
  assert result.stdout_lines == ["2"]
  result = fb(command=command, stdin_lines="hello\nworld\nbye\n") \
    .raise_if_errors()
  if verbose: print(result.stdout_lines)
  assert result.stdout_lines == ["3"]
  if ("--quick" in args):
    n_lines_o = 10000
  else:
    n_lines_o = 1000000
  if (fb is fully_buffered_simple):
    n_lines_e = 500 # Windows blocks if this value is greater than 701
  else:
    n_lines_e = 10000
  result = fb(
    command=command, stdin_lines=[str(i) for i in range(n_lines_o)]) \
    .raise_if_errors()
  if verbose: print(result.stdout_lines)
  assert result.stdout_lines == [str(n_lines_o)]
  command = pyexe \
    + ' -c "import sys; sys.stderr.write(sys.stdin.read())"'
  result = fb(command=command, stdin_lines="Hello\nWorld\nBye\n") \
    .raise_if_output()
  s = StringIO()
  result.show_stderr(out=s, prefix="%(")
  if verbose: sys.stdout.write(s.getvalue())
  assert s.getvalue() == """\
%(Hello
%(World
%(Bye
"""
  cat_command = command = pyexe \
    + ' -c "import sys; sys.stdout.write(sys.stdin.read())"'
  result = fb(command=command, stdin_lines="hello\nworld\nbye\n") \
    .raise_if_errors()
  s = StringIO()
  result.show_stdout(out=s, prefix=">:")
  if verbose: sys.stdout.write(s.getvalue())
  assert s.getvalue() == """\
>:hello
>:world
>:bye
"""
  result = fb(
    command=command, stdin_lines=[str(i) for i in range(n_lines_o)]) \
    .raise_if_errors()
  result.stdout_lines = list(filter(bool, result.stdout_lines))
  if verbose: print(result.stdout_lines[:5], result.stdout_lines[-5:])
  assert len(result.stdout_lines) == n_lines_o
  assert result.stdout_lines[:5] == ["0","1","2","3","4"]
  assert result.stdout_lines[-5:] == [str(s)
    for s in range(n_lines_o-5, n_lines_o)]
  command = pyexe \
    + ' -c "import sys; sys.stderr.write(sys.stdin.read())"'
  result = fb(
    command=command, stdin_lines=[str(i) for i in range(n_lines_e,0,-1)])
  assert len(result.stdout_lines) == 0
  result.stderr_lines = list(filter(bool, result.stderr_lines))
  if verbose: print(result.stderr_lines[:5], result.stderr_lines[-5:])
  assert len(result.stderr_lines) == n_lines_e
  assert result.stderr_lines[:5] == [str(s)
    for s in range(n_lines_e, n_lines_e-5, -1)]
  assert result.stderr_lines[-5:] == ["5","4","3","2","1"]
  command = pyexe + "; ".join((''' -c "\
import sys, os
lines = sys.stdin.read()
sys.stdout.write(lines)
sys.stdout.flush()
lines = list(filter(bool, lines.splitlines()))[:%d]
lines.reverse()
nl = chr(%d)
sys.stderr.write(nl.join(lines)+nl)
sys.stderr.flush()"''' % (n_lines_e, ord("\n"))).splitlines())
  result = fb(
    command=command, stdin_lines=[str(i) for i in range(n_lines_o)])
  result.stdout_lines = list(filter(bool, result.stdout_lines))
  result.stderr_lines = list(filter(bool, result.stderr_lines))
  if verbose: print(result.stdout_lines[:5], result.stdout_lines[-5:])
  if verbose: print(result.stderr_lines[:5], result.stderr_lines[-5:])
  assert len(result.stdout_lines) == n_lines_o
  assert result.stdout_lines[:5] == ["0","1","2","3","4"]
  assert result.stdout_lines[-5:] == [str(s)
    for s in range(n_lines_o-5, n_lines_o)]
  assert len(result.stderr_lines) == n_lines_e
  assert result.stderr_lines[:5] == [str(s)
    for s in range(n_lines_e-1, n_lines_e-6, -1)]
  assert result.stderr_lines[-5:] == ["4","3","2","1","0"]
  result = go(
    command=command, stdin_lines=[str(i) for i in range(n_lines_o)])
  result.stdout_lines = list(filter(bool, result.stdout_lines))
  if verbose: print(result.stdout_lines[:5], result.stdout_lines[-5:])
  assert len(result.stdout_lines) == n_lines_o + n_lines_e
  assert result.stdout_lines[:5] == ["0","1","2","3","4"]
  assert result.stdout_lines[-5:] == ["4","3","2","1","0"]
  #
  try: fb(command="C68649356116218352").raise_if_errors()
  except RuntimeError as e:
    if verbose: print(e)
    # Just check for RuntimeError; there are now additional
    # specific error messages.
    pass
    # assert str(e).startswith("child process stderr output:\n")
  else: raise Exception_expected
  #
  for stdout_splitlines in [True, False]:
    for n,b in [(10,20),(11,23),(12,26),(13,29)]:
      try:
        fb(
          command=cat_command,
          stdin_lines=[str(i) for i in range(n)],
          stdout_splitlines=stdout_splitlines).raise_if_output()
      except RuntimeError as e:
        if verbose: print(e)
        assert str(e).startswith("unexpected child process output:\n")
        if (stdout_splitlines):
          if (n != 13):
            assert str(e).endswith(str(n-1))
          else:
            assert str(e).endswith("  remaining 3 lines omitted.")
        else:
          assert str(e).endswith("  length of output: %d bytes" % b)
      else: raise Exception_expected
  #
  fb(command=cat_command).raise_if_errors_or_output()
  #
  result = fb(command=["nslookup", "localhost"])
  if verbose:
    print(result.stdout_lines)
    print(result.stderr_lines)
  #
  while ("--forever" in args): pass
  #
  print("OK")

if (__name__ == "__main__"):
  exercise()
