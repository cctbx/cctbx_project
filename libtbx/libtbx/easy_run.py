import subprocess_with_fixes as subprocess
import os

class fully_buffered_base(object):

  def raise_if_errors(self):
    assert not self.join_stdout_stderr
    if (len(self.stderr_lines) != 0):
      msg = ["child process stderr output:"]
      msg.append("  command: " + repr(self.command))
      for line in self.stderr_lines:
        msg.append("  " + line)
      raise RuntimeError("\n".join(msg))
    return self

  def raise_if_output(self, show_output_threshold=10):
    if (len(self.stdout_lines) != 0):
      msg = ["unexpected child process output:"]
      msg.append("  command: " + repr(self.command))
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
      raise RuntimeError("\n".join(msg))
    return self

  def raise_if_errors_or_output(self):
    self.raise_if_errors()
    self.raise_if_output()
    return self

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
    self.stdout_lines = child_stdout.read().splitlines()
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
    self.stdout_lines = o.splitlines()
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

def exercise(args=None):
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
  result = fb(command="%s -V" % pyexe).raise_if_output()
  if (verbose): print result.stderr_lines
  assert result.stderr_lines == ["Python " + sys.version.split()[0]]
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
  if ("--quick" in args):
    n_lines_o = 10000
  else:
    n_lines_o = 1000000
  if (fb is fully_buffered_simple):
    n_lines_e = 500 # Windows blocks if this value is greater than 701
  else:
    n_lines_e = 10000
  assert result.stdout_lines == ["3"]
  result = fb(
    command=command, stdin_lines=[str(i) for i in xrange(n_lines_o)]) \
    .raise_if_errors()
  if (verbose): print result.stdout_lines
  assert result.stdout_lines == [str(n_lines_o)]
  cat_command = command = pyexe \
    + ' -c "import sys; sys.stdout.write(sys.stdin.read())"'
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
  else: raise RuntimeError("Exception expected.")
  #
  for n in [10,11,12,13]:
    try:
      fb(
        command=cat_command,
        stdin_lines=[str(i) for i in xrange(n)]).raise_if_output()
    except RuntimeError, e:
      if (verbose): print e
      assert str(e).startswith("unexpected child process output:\n")
      if (n != 13):
        assert str(e).endswith(str(n-1))
      else:
        assert str(e).endswith("remaining 3 lines omitted.")
    else: raise RuntimeError("Exception expected.")
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
