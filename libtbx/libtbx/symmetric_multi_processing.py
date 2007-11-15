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
import sys, os, select, signal
import cStringIO, cPickle

os_fork = getattr(os, "fork", None)
is_available = os_fork is not None

dict_pop = getattr(dict, "pop", None)
if (dict_pop is None):
  def dict_pop(d, key): # Python 2.2 compatibility
    result = d[key]
    del d[key]
    return result

class logging(object):

  def log(self, msg, newline=True):
    if not self.debug: return
    if newline: print msg
    else: print msg,
    sys.stdout.flush()


class parent(logging):

  def __init__(self, polling_timeout, stdout=None, debug=False):
    self.polling_timeout = polling_timeout
    if stdout is None: stdout = sys.stdout
    self.stdout = stdout
    self.children = []
    self.debug = debug

  def spawn(self, func):
    parent_read, child_write = os.pipe()
    child_read, parent_write = os.pipe()
    if (os_fork is None):
      raise NotImplementedError("os.fork not available.")
    pid = os_fork()
    if pid == 0:
      os.close(parent_read)
      os.close(parent_write)
      c = child(func, child_read, child_write, self.debug)
      c.loop()
    os.close(child_read)
    os.close(child_write)
    self.children.append(child_proxy(pid, parent_read, parent_write,
                                     debug=self.debug))
    self.log("spawned %i" % pid)

  def each(self, iterable):
    self.idle_children = list(self.children)
    self.iterable_index_handled_by_pid = {}
    self.result_buffer = {}
    self.next_i = 0
    for i,x in enumerate(iterable):
      block = len(self.idle_children) == 1
      proxy = self.idle_children.pop()
      proxy.send(x)
      self.iterable_index_handled_by_pid[proxy] = i
      for result in self.poll(block):
        yield result
    while len(self.idle_children) != len(self.children):
      for result in self.poll(block=True):
        yield result
    del self.idle_children
    del self.iterable_index_handled_by_pid
    del self.result_buffer
    del self.next_i

  def poll(self, block):
    potential_reads = [ proxy.fin for proxy in self.children ]
    if block:
      self.log("parent waiting for result from a child")
      reads, writes, errors = select.select(potential_reads, [], [])
    else:
      reads, writes, errors = select.select(potential_reads, [], [],
                                            self.polling_timeout)
    for fin in reads:
      proxy = self.children[potential_reads.index(fin)]
      proxy.receive()
      self.idle_children.append(proxy)
      iterable_idx = self.iterable_index_handled_by_pid[proxy]
      self.result_buffer[iterable_idx] = proxy.result
      self.stdout.write(proxy.printout)
      self.log("result #%i: %s from child %i"
               % (iterable_idx, proxy.result, proxy.pid))
    for i in sorted(self.result_buffer.keys()):
      if i == self.next_i:
        self.next_i += 1
        x = self.result_buffer[i]; del self.result_buffer[i]
        yield x
      else:
        break


class child(logging):

  class input_closed: pass

  def __init__(self, func, read_fd, write_fd, debug):
    self.func = func
    self.fin = os.fdopen(read_fd)
    self.fout = os.fdopen(write_fd, 'w')
    self.debug = debug
    signal.signal(signal.SIGPIPE, self.handle_pipe_other_end_closing)

  def handle_pipe_other_end_closing(self, signum, frame):
    self.log("Broken pipe!!")
    code = frame.f_code
    self.log("file %s at line %i in %s"
             % (code.co_filename, frame.f_lineno, code.co_name))
    self.exit()

  def log(self, msg, newline=True):
    super(child, self).log(("child > %i: " % self.fout.fileno())
                           + msg, newline)

  def read(self, length):
    msg = self.fin.read(length)
    if not msg: raise child.input_closed
    assert len(msg) == length
    return msg

  def loop(self):
    while 1:
      try:
        l = int(self.read(128), 16)
        args = tupleize(cPickle.loads(self.read(l)))
        output = cStringIO.StringIO()
        sys.stdout, stdout = output, sys.stdout
        result = self.func(*args)
        sys.stdout = stdout
        pickled_result = easy_pickle.dumps(result)
        printout = output.getvalue()
        msg = ("%128x%s"*2) % (len(pickled_result), pickled_result,
                               len(printout), printout)
        self.fout.write(msg),
        self.fout.flush()
      except child.input_closed:
        self.log("Input closed!!")
        self.exit()

  def exit(self):
    self.fin.close()
    self.fout.close()
    os._exit(0)


class child_proxy(logging):

  def __init__(self, pid, read_fd, write_fd, debug):
    self.pid = pid
    self.fin = os.fdopen(read_fd)
    self.fout = os.fdopen(write_fd, 'w')
    self.debug = debug

  def __hash__(self):
    return hash(self.pid)

  def __eq__(self, other):
    return self.pid == other.pid

  def send(self, x):
    self.log("feeding %s to child %i" % (x, self.pid))
    args = easy_pickle.dumps(x)
    self.fout.write("%128x%s" % (len(args), args))
    self.fout.flush()

  def receive(self):
    l = int(self.fin.read(128), 16)
    self.result = cPickle.loads(self.fin.read(l))
    l = int(self.fin.read(128), 16)
    self.printout = self.fin.read(l)



class parallelized_function(parent):

  def __init__(self, func, n_workers, timeout=0.001, stdout=None,
               debug = False):
    super(parallelized_function, self).__init__(timeout, stdout, debug)
    for i in xrange(n_workers):
      self.spawn(func)

  def __call__(self, iterable):
    for result in super(parallelized_function, self).each(iterable):
      yield result
