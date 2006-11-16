import os

class easy_run(object):

  def __init__(self, command, stdin_lines=None, mode="text"):
    assert mode in ["text", "binary"]
    self.command = command
    child_stdin, child_stdout, child_stderr = os.popen3(command, mode[0])
    if (stdin_lines is not None):
      if (not isinstance(stdin_lines, str)):
        stdin_lines = "\n".join(stdin_lines)
        if (len(stdin_lines) != 0):
          stdin_lines += "\n"
      child_stdin.write(stdin_lines)
    child_stdin.close()
    self.stderr_lines = child_stderr.read().splitlines()
    self.stdout_lines = child_stdout.read().splitlines()
    child_stderr.close()
    child_stdout.close()

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
  if (os.path.isfile("/bin/ls")):
    for command in ["/bin/ls /bin", ("/bin/ls", "/bin")]:
      result = easy_run(command=command).raise_if_errors()
      assert "ls" in result.stdout_lines
      if (verbose): print result.stdout_lines
  if (os.path.isfile("/usr/bin/wc")):
    for command in ["/usr/bin/wc -l", ("/usr/bin/wc", "-l")]:
      result = easy_run(command=command).raise_if_errors()
      if (verbose): print result.stdout_lines
      assert result.stdout_lines[0].strip() == "0"
      result = easy_run(command=command, stdin_lines=["hello"]) \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert result.stdout_lines[0].strip() == "1"
      result = easy_run(command=command, stdin_lines=["hello", "world"]) \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert result.stdout_lines[0].strip() == "2"
      result = easy_run(
        command=command, stdin_lines="hello\nworld\nbye\n") \
        .raise_if_errors()
      if (verbose): print result.stdout_lines
      assert result.stdout_lines[0].strip() == "3"
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
