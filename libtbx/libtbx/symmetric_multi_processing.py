""" Python threading is useless in the context of symmetric multiprocessor
(SMP) machines because the interpreter data structures are not thread-safe.
As a result, a global interpreter lock (GIL) ensures that one and only one
thread can run the interpreter at any time. As a result, it is impossible
to concurrently run Python functions on several processors. The BDL has made
it quite clear the situation will not change in the foreseeable future.

However multiprocessor machines are now mundane. Hence this module. For the
time being, only UNIX systems are supported, through the use of 'fork'.
"""

from __future__ import generators, division
from libtbx.utils import format_cpu_times
from libtbx.utils import tupleize
import sys, os, select
import cStringIO


class parallelized_function(object):

  def __init__(self, func, max_children, timeout=0.001, output=sys.stdout):
    if 'fork' not in os.__dict__:
      raise NotImplementedError("Only works on platforms featuring 'fork',"
                                "i.e. UNIX")
    self.func = func
    self.max_children = max_children
    self.timeout = timeout
    self.output = output
    self.child_pid_for_out_fd = {}
    self.debug = False

  def child_out_fd(self):
    return self.child_pid_for_out_fd.keys()
  child_out_fd = property(child_out_fd)

  def poll(self, block):
    if block and self.debug: print "** waiting for a child to finish **"
    if block:
      inputs, outputs, errors = select.select([], self.child_out_fd, [])
    else:
      inputs, outputs, errors = select.select([], self.child_out_fd, [],
                                              self.timeout)
    if not outputs:
      if self.debug: print "\tno file descriptor is ready"
    else:
      if self.debug: print "\tfile descriptor%s %s %s ready" % (
        ('','s')[len(outputs)>1],
        ', '.join(["%i" % pid for pid in outputs]),
        ('is', 'are')[len(outputs)>1])
    for fd in outputs:
      f = os.fdopen(fd)
      msg = f.read() # block only shortly because of the child buffering
      print >> self.output, msg
      pid = self.child_pid_for_out_fd.pop(fd)
      os.waitpid(pid, 0) # to make sure the child is fully cleaned up

  def run_on_one_processor(self, iterable):
    stdout, sys.stdout = sys.stdout, self.output
    for x in iterable:
      args = tupleize(x)
      yield self.func(*args)
    sys.stdout = stdout

  def run_on_multiple_processors(self, iterable):
    for x in iterable:
      args = tupleize(x)
      r,w = os.pipe()
      pid = os.fork()
      if pid > 0:
        # parent
        if self.debug: print "parent: from %i, child %i: to %i" % (r, pid, w)
        os.close(w)
        self.child_pid_for_out_fd[r] = pid
        self.poll(block=len(self.child_out_fd) == self.max_children)
      else:
        # child
        os.close(r)
        # fully buffer the print-out of the wrapped function
        sys.stdout = cStringIO.StringIO()
        self.func(*args)
        if self.debug: print "\tchild writing to %i" % w
        f = os.fdopen(w, 'w')
        print >> f, sys.stdout.getvalue()
        try:sys.stdout.close()
        except KeyboardInterrupt: raise
        except:pass
        os._exit(0)
    while self.child_out_fd: self.poll(block=False)

def exercise():
  import math, time
  from libtbx.test_utils import approx_equal
  def func(i):
    s = 0
    n = 8e5
    h = math.pi/2/n
    for j in xrange(int(n)):
      s += math.cos(i*j*h)
    s = abs(s*h)
    print "iteration #%i:" % i
    time.sleep(0.05)
    print "\tS=%.2f" % s
    return s

  ref_results = [ abs(math.sin(i*math.pi/2)/i) for i in xrange(1,11) ]
  ref_printout = '\n'.join([ "iteration #%i:\n\tS=%.2f" % (i,s)
                             for i,s in zip(xrange(1,11), ref_results) ])
  output = cStringIO.StringIO()
  try:
    f = parallelized_function(func, max_children=2, output=output)
  except NotImplementedError:
    print "Skipped!"
    return
  results = list(f.run_on_one_processor(xrange(1,11)))
  assert approx_equal(results, ref_results, eps=0.001)
  assert output.getvalue().strip() == ref_printout
  f.output = cStringIO.StringIO()
  f.run_on_multiple_processors(xrange(1,11))
  assert output.getvalue().strip() == ref_printout

  print format_cpu_times()

if __name__ == '__main__':
  exercise()
