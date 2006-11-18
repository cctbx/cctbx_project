try:
  import subprocess
except ImportError, e:
  if (str(e) != "No module named subprocess"): raise
  subprocess = None
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
    p = Popen(
      args=command,
      shell=True,
      bufsize=bufsize,
      stdin=subprocess.PIPE,
      stdout=subprocess.PIPE,
      stderr=stderr,
      universal_newlines=True, # splitlines takes care of this
      close_fds=not subprocess.mswindows)
    o, e = p.communicate(input=stdin_lines)
    self.stdout_lines = o.splitlines()
    if (join_stdout_stderr):
      self.stderr_lines = []
    else:
      self.stderr_lines = e.splitlines()
    self.return_code = p.returncode

if (subprocess is None):
  fully_buffered = fully_buffered_simple
else:
  fully_buffered = fully_buffered_subprocess

  class Popen(subprocess.Popen):
    if (not subprocess.mswindows):
      def _communicate(self, input):
            """Copy of Python 2.5 subprocess.py with patch to fix O(N**2)
               problem:
--- subprocess.py.orig  2006-09-20 11:26:11.000000000 -0700
+++ subprocess.py       2006-11-16 11:57:30.000000000 -0800
@@ -1105,6 +1105,7 @@
                 read_set.append(self.stderr)
                 stderr = []

+            input_done = 0
             while read_set or write_set:
                 rlist, wlist, xlist = select.select(read_set, write_set, [])

@@ -1112,9 +1113,10 @@
                     # When select has indicated that the file is writable,
                     # we can write up to PIPE_BUF bytes without risk
                     # blocking.  POSIX defines PIPE_BUF >= 512
-                    bytes_written = os.write(self.stdin.fileno(), input[:512])
-                    input = input[bytes_written:]
-                    if not input:
+                    bytes_written = os.write(self.stdin.fileno(),
+                                             input[input_done:input_done+512])
+                    input_done += bytes_written
+                    if input_done >= len(input):
                         self.stdin.close()
                         write_set.remove(self.stdin)

            """
            from subprocess import select
            read_set = []
            write_set = []
            stdout = None # Return
            stderr = None # Return

            if self.stdin:
                # Flush stdio buffer.  This might block, if the user has
                # been writing to .stdin in an uncontrolled fashion.
                self.stdin.flush()
                if input:
                    write_set.append(self.stdin)
                else:
                    self.stdin.close()
            if self.stdout:
                read_set.append(self.stdout)
                stdout = []
            if self.stderr:
                read_set.append(self.stderr)
                stderr = []

            input_done = 0
            while read_set or write_set:
                rlist, wlist, xlist = select.select(read_set, write_set, [])

                if self.stdin in wlist:
                    # When select has indicated that the file is writable,
                    # we can write up to PIPE_BUF bytes without risk
                    # blocking.  POSIX defines PIPE_BUF >= 512
                    bytes_written = os.write(self.stdin.fileno(),
                                             input[input_done:input_done+512])
                    input_done += bytes_written
                    if input_done >= len(input):
                        self.stdin.close()
                        write_set.remove(self.stdin)

                if self.stdout in rlist:
                    data = os.read(self.stdout.fileno(), 1024)
                    if data == "":
                        self.stdout.close()
                        read_set.remove(self.stdout)
                    stdout.append(data)

                if self.stderr in rlist:
                    data = os.read(self.stderr.fileno(), 1024)
                    if data == "":
                        self.stderr.close()
                        read_set.remove(self.stderr)
                    stderr.append(data)

            # All data exchanged.  Translate lists into strings.
            if stdout is not None:
                stdout = ''.join(stdout)
            if stderr is not None:
                stderr = ''.join(stderr)

            # Translate newlines, if requested.  We cannot let the file
            # object do the translation: It is based on stdio, which is
            # impossible to combine with select (unless forcing no
            # buffering).
            if self.universal_newlines and hasattr(file, 'newlines'):
                if stdout:
                    stdout = self._translate_newlines(stdout)
                if stderr:
                    stderr = self._translate_newlines(stderr)

            self.wait()
            return (stdout, stderr)

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
