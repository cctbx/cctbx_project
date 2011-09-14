def run(args):
  assert len(args) == 0
  from libtbx import easy_run
  qstat_buffer = easy_run.fully_buffered(command="qstat -g c")
  el = qstat_buffer.stderr_lines
  ol = qstat_buffer.stdout_lines
  if (len(el) != 0):
    print -1
  elif (len(ol) < 3):
    print -2
  elif (   " ".join(ol[0].split())
        != "CLUSTER QUEUE CQLOAD USED AVAIL TOTAL aoACDS cdsuE"):
    print -3
  elif (not ol[1].startswith("----------")):
    print -4
  else:
    sum_available = 0
    for line in ol[2:]:
      flds = line.split()
      assert len(flds) == 7
      sum_available += int(flds[3])
    print sum_available

if (__name__ == "__main__"):
  import sys
  run(sys.argv[1:])
