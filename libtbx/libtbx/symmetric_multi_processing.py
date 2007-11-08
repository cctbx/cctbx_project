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
from libtbx import forward_compatibility
from libtbx.utils import tupleize
from libtbx import easy_pickle
import sys, os, select
import cStringIO, cPickle


class parallelized_function(object):

  def __init__(self, func, max_children, timeout=0.001, printout=None):
    if (printout is None): printout = sys.stdout
    if 'fork' not in os.__dict__:
      raise NotImplementedError("Only works on platforms featuring 'fork',"
                                "i.e. UNIX")
    self.func = func
    self.max_children = max_children
    self.timeout = timeout
    self.printout = printout
    self.child_pid_for_out_fd = {}
    self.index_of_pid = {}
    self.next_i = 0
    self.result_buffer = {}
    self.debug = False

  def poll(self):
    child_out_fd = self.child_pid_for_out_fd.keys()
    if len(child_out_fd) == self.max_children:
      if self.debug: print "** waiting for a child to finish **"
      inputs, outputs, errors = select.select(child_out_fd, [], [])
    else:
      inputs, outputs, errors = select.select(child_out_fd, [], [],
                                              self.timeout)
    for fd in inputs:
      f = os.fdopen(fd)
      msg = f.read() # block only shortly because of the child buffering
      l = int(msg[:128], 16)
      result = cPickle.loads(msg[128:128+l])
      print >> self.printout, msg[128+l:],
      pid = self.child_pid_for_out_fd.pop(fd)
      j = self.index_of_pid.pop(pid)
      os.waitpid(pid, 0) # to make sure the child is fully cleaned up
      if self.debug: print "#%i -> %s" % (j, result)
      self.result_buffer[j] = result
    for i in sorted(self.result_buffer.keys()):
      if i == self.next_i:
        self.next_i += 1
        x = self.result_buffer.pop(i)
        yield x
      else:
        break

  def __call__(self, iterable):
    for i,x in enumerate(iterable):
      args = tupleize(x)
      r,w = os.pipe()
      pid = os.fork()
      if pid > 0:
        # parent
        if self.debug: print "spawn child %i" % pid
        os.close(w)
        self.child_pid_for_out_fd[r] = pid
        self.index_of_pid[pid] = i
        for result in self.poll():
          yield result
      else:
        # child
        os.close(r)
        # fully buffer the print-out of the wrapped function
        output = cStringIO.StringIO()
        sys.stdout, stdout = output, sys.stdout
        result = self.func(*args)
        sys.stdout = stdout
        f = os.fdopen(w, 'w')
        pickled_result = easy_pickle.dumps(result)
        msg = ''.join((
          "%128x" % len(pickled_result),
          pickled_result,
          output.getvalue()))
        print >> f, msg,
        try:f.close()
        except KeyboardInterrupt: raise
        except:pass
        os._exit(0)
    while self.child_pid_for_out_fd:
      for result in self.poll():
        yield result


def exercise(max_children):
  import math, re
  from libtbx.test_utils import approx_equal
  from libtbx.utils import show_times_at_exit, show_times
  show_times_at_exit()
  def func(i):
    s = 0
    n = 1e6/(i*i % 5 + 1)
    h = math.pi/2/n
    for j in xrange(int(n)):
      s += math.cos(i*j*h)
    s = abs(s*h)
    print "%i: S=%.2f" % (i,s)
    return s

  try:
    low, high = 1, 35
    result_pat = re.compile("(\d+): S=(\d\.\d\d)")
    ref_results = [ abs(math.sin(i*math.pi/2)/i) for i in xrange(low,high) ]
    ref_printout = zip(xrange(low,high), [ "%.2f" % s for s in ref_results ])
    f = parallelized_function(func, max_children=max_children,
                              printout=cStringIO.StringIO())
    if 0:
      f.debug = True
      f.timeout = 0.1
    results = list(f(xrange(low,high)))
    if 0:
      print results
      print f.printout.getvalue()
    assert approx_equal(results, ref_results, eps=1e-3)
    matches = [ result_pat.search(li)
                for li in f.printout.getvalue().strip().split('\n') ]
    printout = [ (int(m.group(1)), m.group(2)) for m in matches ]
    printout.sort()
    assert printout == ref_printout
    sys.stdout.flush()
  except NotImplementedError:
    print "Skipped!"

  print "OK"

def run():
  if len(sys.argv[1:]): max_children = int(sys.argv[1])
  else: max_children = 2
  exercise(max_children)

if __name__ == '__main__':
  run()
