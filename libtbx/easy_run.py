import sys
# XXX how early a version can we get away with using the built-in module?
if (sys.version_info[1] >= 7) :
  import subprocess
else :
  try: import subprocess_with_fixes as subprocess
  except ImportError: import subprocess
import sys, os

def _show_lines(lines, out, prefix):
  if (out is None): out = sys.stdout
  for line in lines:
    print >> out, prefix+line

class fully_buffered_base(object):

  def format_errors_if_any(self):
    assert not self.join_stdout_stderr
    if (len(self.stderr_lines) != 0):
      msg = ["child process stderr output:"]
      msg.append("  command: " + repr(self.command))
      for line in self.stderr_lines:
        msg.append("  " + line)
      return "\n".join(msg)
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
        stdin_lines = os.linesep.join(stdin_lines)
        if (len(stdin_lines) != 0):
          stdin_lines += os.linesep
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
        stdin_lines=None,
        join_stdout_stderr=False,
        stdout_splitlines=True,
        bufsize=-1):
    self.command = command
    self.join_stdout_stderr = join_stdout_stderr
    if (not isinstance(command, str)):
      command = subprocess.list2cmdline(command)
    if (stdin_lines is not None):
      if (not isinstance(stdin_lines, str)):
        stdin_lines = os.linesep.join(stdin_lines)
        if (len(stdin_lines) != 0):
          stdin_lines += os.linesep
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
      close_fds=not subprocess.mswindows)
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

def go(command, stdin_lines=None):
  return fully_buffered(
    command=command,
    stdin_lines=stdin_lines,
    join_stdout_stderr=True)

def call(command):
  for s in [sys.stdout, sys.stderr]:
    flush = getattr(s, "flush", None)
    if (flush is not None): flush()
  return subprocess.call(args=command, shell=True)

def exercise(args=None):
  from cStringIO import StringIO
  import sys
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
      if (verbose): print result.stdout_lines
      assert result.stdout_lines == ["hello world"]
  #
  if (os.path.isfile("/bin/ls")):
    for command in ["/bin/ls /bin", ("/bin/ls", "/bin")]:
      result = fb(command=command).raise_if_errors()
      if (verbose): print result.stdout_lines
      assert "ls" in result.stdout_lines
  if (os.path.isfile("/usr/bin/wc")):
    for command in ["/usr/bin/wc -l", ("/usr/bin/wc", "-l")]:
      result = fb(command=command).raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["0"]
      result = fb(command=command, stdin_lines=["hello"]) \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["1"]
      result = fb(command=command, stdin_lines=["hello", "world"]) \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["2"]
      result = fb(command=command, stdin_lines="hello\nworld\nbye\n") \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["3"]
  #
  if (os.name == "nt"):
    result = fb(command="dir /?").raise_if_errors()
    if (verbose): print result.stdout_lines
    assert len(result.stdout_lines) > 0
    windir = os.environ.get("windir", None)
    if (windir is not None and windir.find(" ") < 0):
      result = fb(command="dir "+windir).raise_if_errors()
      if (verbose): print result.stdout_lines
      assert len(result.stdout_lines) > 0
  #
  pyexe = sys.executable
  assert pyexe.count('"') == 0
  pyexe = '"' + pyexe + '"'
  if (os.name == "nt"):
    pyexe = "call " + pyexe
  #
  if (os.environ.has_key("PYTHONPATH")):
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
    if (verbose): print result.stdout_lines
    for line in result.stdout_lines:
      assert not line.startswith("PYTHONPATH") or line == "PYTHONPATH="
  #
  for stdout_splitlines in [True, False]:
    result = fb(
      command="%s -V" % pyexe,
      stdout_splitlines=stdout_splitlines).raise_if_output()
    if (verbose): print result.stderr_lines
    assert result.stderr_lines == ["Python " + sys.version.split()[0]]
    if (stdout_splitlines):
      assert result.stdout_buffer is None
      assert result.stdout_lines == []
    else:
      assert result.stdout_buffer == ""
      assert result.stdout_lines is None
  result = go(command="%s -V" % pyexe)
  if (verbose): print result.stdout_lines
  assert result.stdout_lines == ["Python " + sys.version.split()[0]]
  result = fb(
    command='%s -c "print 3+4"' % pyexe).raise_if_errors()
  if (verbose): print result.stdout_lines
  assert result.stdout_lines == ["7"]
  command = command = pyexe \
    + ' -c "import sys; print len(sys.stdin.read().splitlines())"'
  result = fb(command=command).raise_if_errors()
  if (verbose): print result.stdout_lines
  assert result.stdout_lines == ["0"]
  result = fb(command=command, stdin_lines=["hello"]) \
    .raise_if_errors()
  if (verbose): print result.stdout_lines
  assert result.stdout_lines == ["1"]
  result = fb(command=command, stdin_lines=["hello", "world"]) \
    .raise_if_errors()
  if (verbose): print result.stdout_lines
  assert result.stdout_lines == ["2"]
  result = fb(command=command, stdin_lines="hello\nworld\nbye\n") \
    .raise_if_errors()
  if (verbose): print result.stdout_lines
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
    command=command, stdin_lines=[str(i) for i in xrange(n_lines_o)]) \
    .raise_if_errors()
  if (verbose): print result.stdout_lines
  assert result.stdout_lines == [str(n_lines_o)]
  command = pyexe \
    + ' -c "import sys; sys.stderr.write(sys.stdin.read())"'
  result = fb(command=command, stdin_lines="Hello\nWorld\nBye\n") \
    .raise_if_output()
  s = StringIO()
  result.show_stderr(out=s, prefix="%(")
  if (verbose): sys.stdout.write(s.getvalue())
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
  if (verbose): sys.stdout.write(s.getvalue())
  assert s.getvalue() == """\
>:hello
>:world
>:bye
"""
  result = fb(
    command=command, stdin_lines=[str(i) for i in xrange(n_lines_o)]) \
    .raise_if_errors()
  if (verbose): print result.stdout_lines[:5], result.stdout_lines[-5:]
  assert len(result.stdout_lines) == n_lines_o
  assert result.stdout_lines[:5] == ["0","1","2","3","4"]
  assert result.stdout_lines[-5:] == [str(s)
    for s in xrange(n_lines_o-5, n_lines_o)]
  command = pyexe \
    + ' -c "import sys; sys.stderr.write(sys.stdin.read())"'
  result = fb(
    command=command, stdin_lines=[str(i) for i in xrange(n_lines_e,0,-1)])
  assert len(result.stdout_lines) == 0
  if (verbose): print result.stderr_lines[:5], result.stderr_lines[-5:]
  assert len(result.stderr_lines) == n_lines_e
  assert result.stderr_lines[:5] == [str(s)
    for s in xrange(n_lines_e, n_lines_e-5, -1)]
  assert result.stderr_lines[-5:] == ["5","4","3","2","1"]
  command = pyexe + "; ".join((''' -c "\
import sys, os
lines = sys.stdin.read()
sys.stdout.write(lines)
sys.stdout.flush()
lines = lines.splitlines()[:%d]
lines.reverse()
nl = chr(%d)
sys.stderr.write(nl.join(lines)+nl)
sys.stderr.flush()"''' % (n_lines_e, ord("\n"))).splitlines())
  result = fb(
    command=command, stdin_lines=[str(i) for i in xrange(n_lines_o)])
  if (verbose): print result.stdout_lines[:5], result.stdout_lines[-5:]
  if (verbose): print result.stderr_lines[:5], result.stderr_lines[-5:]
  assert len(result.stdout_lines) == n_lines_o
  assert result.stdout_lines[:5] == ["0","1","2","3","4"]
  assert result.stdout_lines[-5:] == [str(s)
    for s in xrange(n_lines_o-5, n_lines_o)]
  assert len(result.stderr_lines) == n_lines_e
  assert result.stderr_lines[:5] == [str(s)
    for s in xrange(n_lines_e-1, n_lines_e-6, -1)]
  assert result.stderr_lines[-5:] == ["4","3","2","1","0"]
  result = go(
    command=command, stdin_lines=[str(i) for i in xrange(n_lines_o)])
  if (verbose): print result.stdout_lines[:5], result.stdout_lines[-5:]
  assert len(result.stdout_lines) == n_lines_o + n_lines_e
  assert result.stdout_lines[:5] == ["0","1","2","3","4"]
  assert result.stdout_lines[-5:] == ["4","3","2","1","0"]
  #
  try: fb(command="C68649356116218352").raise_if_errors()
  except RuntimeError, e:
    if (verbose): print e
    assert str(e).startswith("child process stderr output:\n")
  else: raise Exception_expected
  #
  for stdout_splitlines in [True, False]:
    for n,b in [(10,20),(11,23),(12,26),(13,29)]:
      try:
        fb(
          command=cat_command,
          stdin_lines=[str(i) for i in xrange(n)],
          stdout_splitlines=stdout_splitlines).raise_if_output()
      except RuntimeError, e:
        if (verbose): print e
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
  if (verbose):
    print result.stdout_lines
    print result.stderr_lines
  #
  while ("--forever" in args): pass
  #
  print "OK"

if (__name__ == "__main__"):
  exercise()
