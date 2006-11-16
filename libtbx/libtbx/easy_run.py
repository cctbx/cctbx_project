import os

class easy_run(object):
  """\
Executes command, sends stdin_lines (str or sequence), then reads
stdout_lines first, stderr_lines second.

The constructor may deadlock if the I/O buffers are too small to allow
the blocking write and reads in the given sequence. Specifically,
stdin_lines may be too big, or there may be too many stderr_lines,
but there can be any number of stdout_lines. Tested under Linux,
Mac OS X, Windows XP, IRIX, Tru64 Unix. The tests below are known
to work portably with stdin_lines up to 1000000, stderr_lines up
to 701. I.e. this simple implementation should cover most practical
situations.

If you'd like to contribute a completely robust and portable
alternative implementation, please email rwgk@yahoo.com.
  """

  def __init__(self, command, stdin_lines=None, mode="text", bufsize=-1):
    assert mode in ["text", "binary"]
    self.command = command
    child_stdin, child_stdout, child_stderr = os.popen3(
      command, mode[0], bufsize)
    if (stdin_lines is not None):
      if (not isinstance(stdin_lines, str)):
        stdin_lines = "\n".join(stdin_lines)
        if (len(stdin_lines) != 0):
          stdin_lines += "\n"
      child_stdin.write(stdin_lines)
    child_stdin.close()
    self.stdout_lines = child_stdout.read().splitlines()
    self.stderr_lines = child_stderr.read().splitlines()
    child_stdout.close()
    child_stderr.close()

  def raise_if_errors(self):
    if (len(self.stderr_lines) != 0):
      msg = ["easy_run errors:"]
      msg.append("  command: " + repr(self.command))
      for line in self.stderr_lines:
        msg.append("  " + line)
      raise RuntimeError("\n".join(msg))
    return self

def exercise():
  import sys
  verbose = "--verbose" in sys.argv[1:]
  #
  if (os.path.isfile("/bin/ls")):
    for command in ["/bin/ls /bin", ("/bin/ls", "/bin")]:
      result = easy_run(command=command).raise_if_errors()
      assert "ls" in result.stdout_lines
      if (verbose): print result.stdout_lines
  if (os.path.isfile("/usr/bin/wc")):
    for command in ["/usr/bin/wc -l", ("/usr/bin/wc", "-l")]:
      result = easy_run(command=command).raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["0"]
      result = easy_run(command=command, stdin_lines=["hello"]) \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["1"]
      result = easy_run(command=command, stdin_lines=["hello", "world"]) \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["2"]
      result = easy_run(
        command=command, stdin_lines="hello\nworld\nbye\n") \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert [s.strip() for s in result.stdout_lines] == ["3"]
  #
  if (os.name == "nt"):
    result = easy_run(command="dir /?").raise_if_errors()
    if (verbose): print result.stdout_lines
    assert len(result.stdout_lines) > 0
    windir = os.environ.get("windir", None)
    if (windir is not None and windir.find(" ") < 0):
      result = easy_run(command="dir "+windir).raise_if_errors()
      if (verbose): print result.stdout_lines
      assert len(result.stdout_lines) > 0
  pyexe = sys.executable
  if (pyexe.find(" ") < 0):
    if (os.environ.has_key("PYTHONPATH")):
      del os.environ["PYTHONPATH"]
    if (os.name == "nt"):
      result = easy_run(command="set").raise_if_errors()
    elif (os.path.isfile("/usr/bin/printenv")):
      result = easy_run(command="/usr/bin/printenv").raise_if_errors()
    else:
      result = None
    if (result is not None):
      if (verbose): print result.stdout_lines
      for line in result.stdout_lines:
        assert not line.startswith("PYTHONPATH")
    result = easy_run(command="%s -V" % pyexe)
    if (verbose): print result.stderr_lines
    assert result.stderr_lines == ["Python " + sys.version.split()[0]]
    result = easy_run(command='%s -c "print 3+4"' % pyexe).raise_if_errors()
    if (verbose): print result.stdout_lines
    assert result.stdout_lines == ["7"]
    command = pyexe \
      + ' -c "import sys; print len(sys.stdin.read().splitlines())"'
    result = easy_run(command=command).raise_if_errors()
    if (verbose): print result.stdout_lines
    assert result.stdout_lines == ["0"]
    result = easy_run(command=command, stdin_lines=["hello"]) \
      .raise_if_errors()
    if (verbose): print result.stdout_lines
    assert result.stdout_lines == ["1"]
    result = easy_run(command=command, stdin_lines=["hello", "world"]) \
      .raise_if_errors()
    if (verbose): print result.stdout_lines
    assert result.stdout_lines == ["2"]
    result = easy_run(
      command=command, stdin_lines="hello\nworld\nbye\n") \
      .raise_if_errors()
    if (verbose): print result.stdout_lines
    if ("--quick" in sys.argv[1:]):
      n_lines_o = 10000
    else:
      n_lines_o = 1000000
    n_lines_e = 500 # Windows blocks if this value is greater than 701
    assert result.stdout_lines == ["3"]
    result = easy_run(
      command=command, stdin_lines=[str(i) for i in xrange(n_lines_o)]) \
      .raise_if_errors()
    if (verbose): print result.stdout_lines
    assert result.stdout_lines == [str(n_lines_o)]
    command = pyexe \
      + ' -c "import sys; sys.stdout.write(sys.stdin.read())"'
    result = easy_run(
      command=command, stdin_lines=[str(i) for i in xrange(n_lines_o)]) \
      .raise_if_errors()
    if (verbose): print result.stdout_lines[:5], result.stdout_lines[-5:]
    assert len(result.stdout_lines) == n_lines_o
    assert result.stdout_lines[:5] == ["0","1","2","3","4"]
    assert result.stdout_lines[-5:] == [str(s)
      for s in xrange(n_lines_o-5, n_lines_o)]
    command = pyexe \
      + ' -c "import sys; sys.stderr.write(sys.stdin.read())"'
    result = easy_run(
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
sys.stderr.write(os.linesep.join(lines)+os.linesep)
sys.stderr.flush()"''' % n_lines_e).splitlines())
    result = easy_run(
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
  #
  try: easy_run(command="C68649356116218352").raise_if_errors()
  except RuntimeError, e:
    assert str(e).startswith("easy_run errors:\n")
  else: raise RuntimeError("Exception expected.")
  #
  if (os.path.isfile("/usr/bin/nslookup")):
    result = easy_run(command=["/usr/bin/nslookup", "localhost"])
    if (verbose):
      print result.stdout_lines
      print result.stderr_lines
  #
  while ("--forever" in sys.argv[1:]): pass
  #
  print "OK"

if (__name__ == "__main__"):
  exercise()
