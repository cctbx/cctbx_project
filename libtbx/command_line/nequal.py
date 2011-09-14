from libtbx.utils import Usage
import fileinput
import sys, os

def run(command_name=os.environ.get(
          "LIBTBX_DISPATCHER_NAME", "libtbx.nequal")):
  if (sys.argv[1:] in [["-h"], ["--help"]]):
    raise Usage("%s [file ...]" % command_name + """
  Similar to the Unix uniq command, but each output line is
  prefixed with the number of identical consecutive lines.
  Example command:
    grep Warning log_file | sort | %s
  Example output:
    12: Warning: missing file'
     9: Warning: missing directory'
     1: Warning: unknown file
    Number of lines shown: 3
    Sum of counts: 22""" % command_name)
  buffer = []
  prev = None
  n = 0
  for line in fileinput.input():
    if (prev is None):
      prev = line
      n = 1
    elif (line != prev):
      buffer.append((n, prev))
      prev = line
      n = 1
    else:
      n += 1
  if (n != 0):
    buffer.append((n, prev))
  if (len(buffer) != 0):
    def cmp_buffer_entries(a, b):
      result = cmp(b[0], a[0])
      if (result == 0):
        result = cmp(a[1], b[1])
      return result
    buffer.sort(cmp_buffer_entries)
    sum_n = 0
    n_fmt = "%%%dd: " % len("%d" % buffer[0][0])
    for n,line in buffer:
      sys.stdout.write(n_fmt % n + line)
      sum_n += n
    print "Number of lines shown:", len(buffer)
    print "Sum of counts:", sum_n

if (__name__ == "__main__"):
  run()
