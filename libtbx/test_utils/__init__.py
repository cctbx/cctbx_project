from libtbx.option_parser import option_parser
from libtbx.utils import Sorry
from libtbx.str_utils import show_string
from libtbx import easy_run
from libtbx import introspection
import libtbx.load_env
import difflib
from stdlib import math
import sys, os
import time

try:
  import threading
except ImportError:
  threading = None
else:
  import Queue

diff_function = getattr(difflib, "unified_diff", difflib.ndiff)

def compose_tmp_suffix(suffix, frames_back=0):
  from libtbx.introspection import caller_location
  caller = caller_location(frames_back=frames_back+1)
  s = os.path.basename(caller.file_name)
  if (s.lower().endswith(".py")): s = s[:-3]
  s += "_" + str(caller.line_number)
  return "_" + s + suffix

def open_tmp_file(suffix="", mode="w"):
  import tempfile
  (fd, name) = tempfile.mkstemp(
    suffix=compose_tmp_suffix(suffix, frames_back=1),
    dir=".")
  os.close(fd)
  return open(name, mode)

def open_tmp_directory(suffix=""):
  import tempfile
  return tempfile.mkdtemp(
    suffix=compose_tmp_suffix(suffix, frames_back=1),
    dir=".")

class pickle_detector(object):
  def __init__(O):
    O.unpickled_counter = None
    O.counter = 0
  def __getstate__(O):
    O.counter += 1
    return {"counter": O.counter}
  def __setstate__(O, state):
    assert len(state) == 1
    O.unpickled_counter = state["counter"]
    O.counter = 0

Exception_expected = RuntimeError("Exception expected.")
Exception_not_expected = RuntimeError("Exception not expected.")

class Default: pass

def run_tests(build_dir, dist_dir, tst_list, display_times=False):
  if display_times:
    t0=time.time()
    start = time.asctime()
  libtbx.env.full_testing = True
  args = [arg.lower() for arg in sys.argv[1:]]
  command_line = (option_parser(
    usage="run_tests [-j n]",
    description="Run several threads in parallel, each picking and then"
                " running tests one at a time.")
    .option("-j", "--threads",
      action="store",
      type="int",
      default=1,
      help="number of threads",
      metavar="INT")
    .option("-v", "--verbose",
      action="store_true",
      default=False)
    .option("-q", "--quick",
      action="store_true",
      default=False)
    .option("-g", "--valgrind",
      action="store_true",
      default=False)
  ).process(args=args, max_nargs=0)
  co = command_line.options
  if (threading is None or co.threads == 1):
    for cmd in iter_tests_cmd(co, build_dir, dist_dir, tst_list):
      print cmd
      sys.stdout.flush()
      easy_run.call(command=cmd)
      print
      sys.stderr.flush()
      sys.stdout.flush()
  else:
    cmd_queue = Queue.Queue()
    for cmd in iter_tests_cmd(co, build_dir, dist_dir, tst_list):
      cmd_queue.put(cmd)
    threads_pool = []
    log_queue = Queue.Queue()
    interrupted = threading.Event()
    for i in xrange(co.threads):
      working_dir = os.tempnam()
      os.mkdir(working_dir)
      t = threading.Thread(
        target=make_pick_and_run_tests(working_dir, interrupted,
                                       cmd_queue, log_queue))
      threads_pool.append(t)
    for t in threads_pool:
      t.start()
    finished_threads = 0
    while(1):
      try:
        log = log_queue.get()
        if isinstance(log, tuple):
          msg = log[0]
          print >> sys.stderr, "\n +++ thread %s +++\n" % msg
          finished_threads += 1
          if finished_threads == co.threads: break
        else:
          print log
      except KeyboardInterrupt:
        print
        print "********************************************"
        print "** Received Keyboard Interrupt            **"
        print "** Waiting for running tests to terminate **"
        print "********************************************"
        interrupted.set()
        break
  if display_times:
    print "TIME (%s) (%s) %7.2f %s" % (start,
                                       time.asctime(),
                                       time.time()-t0,
                                       tst_list,
                                       )

def make_pick_and_run_tests(working_dir, interrupted,
                            cmd_queue, log_queue):
  def func():
    while(not interrupted.isSet()):
      try:
        cmd = cmd_queue.get(block=True, timeout=0.2)
        log = "\n%s\n" % cmd
        proc = easy_run.subprocess.Popen(cmd,
                                  shell=True,
                                  stdout=easy_run.subprocess.PIPE,
                                  stderr=easy_run.subprocess.STDOUT,
                                  cwd=working_dir)
        proc.wait()
        log += "\n%s" % proc.stdout.read()
        log_queue.put(log)
      except Queue.Empty:
        log_queue.put( ("done",) )
        break
      except KeyboardInterrupt:
        break
  return func

def iter_tests_cmd(co, build_dir, dist_dir, tst_list):
  for tst in tst_list:
    cmd_args = ""
    if (type(tst) == type([])):
      if (co.verbose):
        cmd_args = " " + " ".join(["--Verbose"] + tst[1:])
      elif (co.quick):
        cmd_args = " " + " ".join(tst[1:])
      tst = tst[0]
    elif (co.verbose):
      continue
    if (tst.startswith("$B")):
      tst_path = tst.replace("$B", build_dir)
    else:
      tst_path = tst.replace("$D", dist_dir)
    assert tst_path.find("$") < 0
    assert tst_path.find('"') < 0
    tst_path = os.path.normpath(tst_path)
    cmd = ""
    if (tst_path.endswith(".py")):
      if (co.valgrind):
        cmd = "libtbx.valgrind "
      cmd += "libtbx.python "
    else:
      if (co.valgrind):
        cmd = os.environ.get(
          "LIBTBX_VALGRIND", "valgrind --tool=memcheck") + " "
    cmd += '"' + tst_path + '"'
    cmd += cmd_args
    yield cmd

def approx_equal_core(a1, a2, eps, multiplier, out, prefix):
  if isinstance(a1, str) or isinstance(a1, unicode):
    return a1 == a2
  if hasattr(a1, "__len__"): # traverse list
    if (len(a1) != len(a2)):
      raise AssertionError(
        "approx_equal ERROR: len(a1) != len(a2): %d != %d" % (
          len(a1), len(a2)))
    for i in xrange(len(a1)):
      if not approx_equal_core(
                a1[i], a2[i], eps, multiplier, out, prefix+"  "):
        return False
    return True
  is_complex_1 = isinstance(a1, complex)
  is_complex_2 = isinstance(a2, complex)
  if (is_complex_1 and is_complex_2): # complex & complex
    if not approx_equal_core(
              a1.real, a2.real, eps, multiplier, out, prefix+"real "):
      return False
    if not approx_equal_core(
              a1.imag, a2.imag, eps, multiplier, out, prefix+"imag "):
      return False
    return True
  elif (is_complex_1): # complex & number
    if not approx_equal_core(
              a1.real, a2, eps, multiplier, out, prefix+"real "):
      return False
    if not approx_equal_core(
              a1.imag, 0, eps, multiplier, out, prefix+"imag "):
      return False
    return True
  elif (is_complex_2): # number & complex
    if not approx_equal_core(
              a1, a2.real, eps, multiplier, out, prefix+"real "):
      return False
    if not approx_equal_core(
              0, a2.imag, eps, multiplier, out, prefix+"imag "):
      return False
    return True
  ok = True
  d = a1 - a2
  if (abs(d) > eps):
    if (multiplier is None):
      ok = False
    else:
      am = max(a1,a2) * multiplier
      d = (am - d) - am
      if (d != 0):
        ok = False
  if (out is not None):
    annotation = ""
    if (not ok):
      annotation = " approx_equal ERROR"
    print >> out, prefix + str(a1) + annotation
    print >> out, prefix + str(a2) + annotation
    print >> out, prefix.rstrip()
    return True
  return ok

def approx_equal(a1, a2, eps=1.e-6, multiplier=1.e10, out=Default, prefix=""):
  ok = approx_equal_core(a1, a2, eps, multiplier, None, prefix)
  if (not ok and out is not None):
    if (out is Default): out = sys.stdout
    print >> out, prefix + "approx_equal eps:", eps
    print >> out, prefix + "approx_equal multiplier:", multiplier
    assert approx_equal_core(a1, a2, eps, multiplier, out, prefix)
  return ok

def not_approx_equal(a1, a2, eps=1.e-6, multiplier=1.e10):
  return not approx_equal(a1, a2, eps, multiplier, out=None)

def eps_eq_core(a1, a2, eps, out, prefix):
  if (hasattr(a1, "__len__")): # traverse list
    assert len(a1) == len(a2)
    for i in xrange(len(a1)):
      if (not eps_eq_core(a1[i], a2[i], eps, out, prefix+"  ")):
        return False
    return True
  if (isinstance(a1, complex)): # complex numbers
    if (not eps_eq_core(a1.real, a2.real, eps, out, prefix+"real ")):
      return False
    if (not eps_eq_core(a1.imag, a2.imag, eps, out, prefix+"imag ")):
      return False
    return True
  ok = True
  if (a1 == 0 or a2 == 0):
    if (abs(a1-a2) > eps):
      ok = False
  else:
    l1 = round(math.log(abs(a1)))
    l2 = round(math.log(abs(a2)))
    m = math.exp(-max(l1, l2))
    if (abs(a1*m-a2*m) > eps):
      ok = False
  if (out is not None):
    annotation = ""
    if (not ok):
      annotation = " eps_eq ERROR"
    print >> out, prefix + str(a1) + annotation
    print >> out, prefix + str(a2) + annotation
    print >> out, prefix.rstrip()
    return True
  return ok

def eps_eq(a1, a2, eps=1.e-6, out=Default, prefix=""):
  ok = eps_eq_core(a1, a2, eps, None, prefix)
  if (not ok and out is not None):
    if (out is Default): out = sys.stdout
    print >> out, prefix + "eps_eq eps:", eps
    assert eps_eq_core(a1, a2, eps, out, prefix)
  return ok

def not_eps_eq(a1, a2, eps=1.e-6):
  return not eps_eq(a1, a2, eps, None)

def is_below_limit(
      value,
      limit,
      eps=1.e-6,
      info_low_eps=None,
      out=Default,
      info_prefix="INFO LOW VALUE: "):
  if (isinstance(value, (int, float)) and value < limit + eps):
    if (info_low_eps is not None and value < limit - info_low_eps):
      if (out is not None):
        if (out is Default): out = sys.stdout
        introspection.show_stack(
          frames_back=1, reverse=True, out=out, prefix=info_prefix)
        print >> out, \
          "%sis_below_limit(value=%s, limit=%s, info_low_eps=%s)" % (
            info_prefix, str(value), str(limit), str(info_low_eps))
    return True
  if (out is not None):
    if (out is Default): out = sys.stdout
    print >> out, "ERROR:", \
      "is_below_limit(value=%s, limit=%s, eps=%s)" % (
        str(value), str(limit), str(eps))
  return False

def is_above_limit(
      value,
      limit,
      eps=1.e-6,
      info_high_eps=None,
      out=Default,
      info_prefix="INFO HIGH VALUE: "):
  if (isinstance(value, (int, float)) and value > limit - eps):
    if (info_high_eps is not None and value > limit + info_high_eps):
      if (out is not None):
        if (out is Default): out = sys.stdout
        introspection.show_stack(
          frames_back=1, reverse=True, out=out, prefix=info_prefix)
        print >> out, \
          "%sis_above_limit(value=%s, limit=%s, info_high_eps=%s)" % (
            info_prefix, str(value), str(limit), str(info_high_eps))
    return True
  if (out is not None):
    if (out is Default): out = sys.stdout
    print >> out, "ERROR:", \
      "is_above_limit(value=%s, limit=%s, eps=%s)" % (
        str(value), str(limit), str(eps))
  return False

def precision_approx_equal(self,other,precision=24):
  # Use concepts from IEEE-754 to determine if the difference between
  # two numbers is within floating point error. Precision is expressed in
  # bits (single precision=24; double precision=53). Not within scope to
  # do this for double precision literals; only interested in the case
  # where the data are from a ~single precision digital-analog converter.
  if precision > 52: raise ValueError()
  if self==other:
    return True
  if (self > 0.) != (other > 0.) :
    return False
  #compute the exponent
  import math
  T = abs(self)
  Np = math.floor(precision-math.log(T,2))
  significand = int(T * 2**Np)
  val1 = significand/(2**Np) # nearest floating point representation of self
  val2 = (significand+1)/(2**Np) # next-nearest
  return abs(T-abs(other)) <= abs(val1-val2)

def show_diff(a, b, selections=None, expected_number_of_lines=None):
  if (not isinstance(a, str)):
    a = "\n".join(a)+"\n"
  if (selections is None):
    assert expected_number_of_lines is None
  else:
    a_lines = a.splitlines(1)
    a = []
    for selection in selections:
      for i in selection:
        if (i < 0): i += len(a_lines)
        a.append(a_lines[i])
      a.append("...\n")
    a = "".join(a[:-1])
  if (a != b):
    if (not a.endswith("\n") or not b.endswith("\n")):
      a += "\n"
      b += "\n"
    print "".join(diff_function(b.splitlines(1), a.splitlines(1)))
    return True
  if (    expected_number_of_lines is not None
      and len(a_lines) != expected_number_of_lines):
    print \
      "show_diff: expected_number_of_lines != len(a.splitlines()): %d != %d" \
        % (expected_number_of_lines, len(a_lines))
    return True
  return False

def block_show_diff(lines, expected, last_startswith=False):
  if (isinstance(lines, str)):
    lines = lines.splitlines()
  if (isinstance(expected, str)):
    expected = expected.splitlines()
  assert len(expected) > 1
  def raise_not_found():
    print "block_show_diff() lines:"
    print "-"*80
    print "\n".join(lines)
    print "-"*80
    raise AssertionError('Expected line not found: "%s"' % eline)
  eline = expected[0]
  for i,line in enumerate(lines):
    if (line == eline):
      lines = lines[i:]
      break
  else:
    raise_not_found()
  eline = expected[-1]
  for i,line in enumerate(lines):
    if (last_startswith):
      if (line.startswith(eline)):
        lines = lines[:i]
        expected = expected[:-1]
        break
    else:
      if (line == eline):
        lines = lines[:i+1]
        break
  else:
    raise_not_found()
  lines = "\n".join(lines)+"\n"
  expected = "\n".join(expected)+"\n"
  return show_diff(lines, expected)

def blocks_show_diff(lines, expected, block_sperator="\n...\n"):
  assert isinstance(expected, str)
  result = False
  expected_blocks = expected.split(block_sperator)
  for expected in expected_blocks:
    if (block_show_diff(lines=lines, expected=expected)):
      result = True
  return result

def anchored_block_show_diff(lines, anchor_index, expected):
  if (isinstance(lines, str)):
    lines = lines.splitlines()
  if (isinstance(expected, str)):
    expected = expected.splitlines()
  assert len(expected) > 0
  if (anchor_index < 0):
    anchor_index = len(lines) + anchor_index - len(expected) + 1
  lines = "\n".join(lines[anchor_index:anchor_index+len(expected)])+"\n"
  expected = "\n".join(expected)+"\n"
  return show_diff(lines, expected)

def contains_substring(
      actual,
      expected,
      failure_prefix="contains_substring() "):
  assert isinstance(actual, str)
  assert isinstance(expected, str)
  if (actual.find(expected) < 0):
    print "%sFAILURE:" % failure_prefix
    def show(s):
      print "v"*79
      if (s.endswith("\n") or s.endswith(os.linesep)):
        sys.stdout.write(s)
      else:
        print s
      print "^"*79
    show(actual)
    print "  ACTUAL ----^"
    print "EXPECTED ----v"
    show(expected)
    return False
  return True

def contains_lines(lines, expected):
  return contains_substring(
    actual=lines, expected=expected, failure_prefix="contains_lines() ")

class RunCommandError(RuntimeError): pass

def _check_command_output(
    lines=None,
    file_name=None,
    show_command_if_error=None,
    sorry_expected=False):
  assert [lines, file_name].count(None) == 1
  if (lines is None):
    lines = open(file_name).read().splitlines()
  def show_and_raise(detected):
    if (show_command_if_error):
      print show_command_if_error
      print
    print "\n".join(lines)
    msg = detected + " detected in output"
    if (file_name is None):
      msg += "."
    else:
      msg += ": " + show_string(file_name)
    raise RunCommandError(msg)
  have_sorry = False
  for line in lines:
    if (line == "Traceback (most recent call last):"):
      show_and_raise(detected="Traceback")
    if (line.startswith("Sorry:")):
      have_sorry = True
  if (have_sorry and not sorry_expected):
    show_and_raise(detected='"Sorry:"')

def run_command(
      command,
      verbose=0,
      buffered=True,
      log_file_name=None,
      stdout_file_name=None,
      result_file_names=[],
      show_diff_log_stdout=False,
      sorry_expected=False,
      join_stdout_stderr=False):
  """\
This function starts another process to run command, with some
pre-call and post-call processing.
Before running command, the expected output files are removed:

  log_file_name
  stdout_file_name
  result_file_names

After command is finished, log_file_name and stdout_file_name are scanned
for Traceback and Sorry. An exception is raised if there are any
matches. sorry_expected=True suppresses the scanning for Sorry.

With buffered=True easy_run.fully_buffered() is used. If there
is any output to stderr of the child process, an exception is
raised. The run_command() return value is the result of the
easy_run.fully_buffered() call.

With buffered=False easy_run.call() is used. I.e. stdout and stderr
of the command are connected to stdout and stderr of the parent
process. stderr is not checked. The run_command() return value is None.

It is generally best to use buffered=True, and even better not to
use this function at all if command is another Python script. It
is better to organize the command script such that it can be called
directly from within the same Python process running the unit tests.
"""
  assert verbose >= 0
  if (verbose > 0):
    print command
    print
    show_command_if_error = None
  else:
    show_command_if_error = command
  all_file_names = [log_file_name, stdout_file_name] + result_file_names
  for file_name in all_file_names:
    if (file_name is None): continue
    if (os.path.isfile(file_name)): os.remove(file_name)
    if (os.path.exists(file_name)):
      raise RunCommandError(
        "Unable to remove file: %s" % show_string(file_name))
  if (buffered):
    sys.stdout.flush()
    sys.stderr.flush()
    cmd_result = easy_run.fully_buffered(
      command=command,
      join_stdout_stderr=join_stdout_stderr)
    if (len(cmd_result.stderr_lines) != 0):
      if (verbose == 0):
        print command
        print
      print "\n".join(cmd_result.stdout_lines)
      cmd_result.raise_if_errors()
    _check_command_output(
      lines=cmd_result.stdout_lines,
      show_command_if_error=show_command_if_error,
      sorry_expected=sorry_expected)
  else:
    easy_run.call(command=command)
    cmd_result = None
  for file_name in [log_file_name, stdout_file_name]:
    if (file_name is None or not os.path.isfile(file_name)): continue
    _check_command_output(
      file_name=file_name,
      show_command_if_error=show_command_if_error,
      sorry_expected=sorry_expected)
  for file_name in all_file_names:
    if (file_name is None): continue
    if (not os.path.isfile(file_name)):
      raise RunCommandError(
        "Missing output file: %s" % show_string(file_name))
  if (verbose > 1 and cmd_result is not None):
    print "\n".join(cmd_result.stdout_lines)
    print
  if (    show_diff_log_stdout
      and log_file_name is not None
      and stdout_file_name is not None):
    if (verbose > 0):
      print "diff %s %s" % (show_string(log_file_name),
                            show_string(stdout_file_name))
      print
    if (show_diff(open(log_file_name).read(), open(stdout_file_name).read())):
      introspection.show_stack(
        frames_back=1, reverse=True, prefix="INFO_LOG_STDOUT_DIFFERENCE: ")
      print "ERROR_LOG_STDOUT_DIFFERENCE"
  sys.stdout.flush()
  return cmd_result

def exercise():
  from cStringIO import StringIO
  assert approx_equal(1, 1)
  out = StringIO()
  assert not approx_equal(1, 0, out=out)
  assert not show_diff(out.getvalue().replace("1e-006", "1e-06"), """\
approx_equal eps: 1e-06
approx_equal multiplier: 10000000000.0
1 approx_equal ERROR
0 approx_equal ERROR

""")
  out = StringIO()
  assert not approx_equal(1, 2, out=out)
  assert not show_diff(out.getvalue().replace("1e-006", "1e-06"), """\
approx_equal eps: 1e-06
approx_equal multiplier: 10000000000.0
1 approx_equal ERROR
2 approx_equal ERROR

""")
  out = StringIO()
  assert not approx_equal(1, 1+1.e-5, out=out)
  assert approx_equal(1, 1+1.e-6)
  out = StringIO()
  assert not approx_equal(0, 1.e-5, out=out)
  assert approx_equal(0, 1.e-6)
  out = StringIO()
  assert not approx_equal([[0,1],[2j,3]],[[0,1],[-2j,3]], out=out, prefix="$%")
  assert not show_diff(out.getvalue().replace("1e-006", "1e-06"), """\
$%approx_equal eps: 1e-06
$%approx_equal multiplier: 10000000000.0
$%    0
$%    0
$%
$%    1
$%    1
$%
$%    real 0.0
$%    real 0.0
$%    real
$%    imag 2.0 approx_equal ERROR
$%    imag -2.0 approx_equal ERROR
$%    imag
$%    3
$%    3
$%
""")
  assert eps_eq(1, 1)
  out = StringIO()
  assert not eps_eq(1, 0, out=out)
  assert not show_diff(out.getvalue().replace("1e-006", "1e-06"), """\
eps_eq eps: 1e-06
1 eps_eq ERROR
0 eps_eq ERROR

""")
  out = StringIO()
  assert not eps_eq(1, 2, out=out)
  assert not show_diff(out.getvalue().replace("1e-006", "1e-06"), """\
eps_eq eps: 1e-06
1 eps_eq ERROR
2 eps_eq ERROR

""")
  out = StringIO()
  assert not eps_eq(1, 1+1.e-5, out=out)
  assert eps_eq(1, 1+1.e-6)
  out = StringIO()
  assert not eps_eq(0, 1.e-5, out=out)
  assert eps_eq(0, 1.e-6)
  out = StringIO()
  assert not eps_eq([[0,1],[2j,3]],[[0,1],[-2j,3]], out=out, prefix="$%")
  assert not show_diff(out.getvalue().replace("1e-006", "1e-06"), """\
$%eps_eq eps: 1e-06
$%    0
$%    0
$%
$%    1
$%    1
$%
$%    real 0.0
$%    real 0.0
$%    real
$%    imag 2.0 eps_eq ERROR
$%    imag -2.0 eps_eq ERROR
$%    imag
$%    3
$%    3
$%
""")
  assert is_below_limit(value=5, limit=10, eps=2)
  out = StringIO()
  assert is_below_limit(value=5, limit=10, eps=2, info_low_eps=1, out=out)
  assert not show_diff(out.getvalue(), """\
INFO LOW VALUE: is_below_limit(value=5, limit=10, info_low_eps=1)
""", selections=[[-1]], expected_number_of_lines=3)
  out = StringIO()
  assert not is_below_limit(value=15, limit=10, eps=2, out=out)
  assert not show_diff(out.getvalue(), """\
ERROR: is_below_limit(value=15, limit=10, eps=2)
""")
  out = StringIO()
  assert not is_below_limit(value=None, limit=3, eps=1, out=out)
  assert not is_below_limit(value=None, limit=-3, eps=1, out=out)
  assert not show_diff(out.getvalue(), """\
ERROR: is_below_limit(value=None, limit=3, eps=1)
ERROR: is_below_limit(value=None, limit=-3, eps=1)
""")
  assert is_above_limit(value=10, limit=5, eps=2)
  out = StringIO()
  assert is_above_limit(value=10, limit=5, eps=2, info_high_eps=1, out=out)
  assert not show_diff(out.getvalue(), """\
INFO HIGH VALUE: is_above_limit(value=10, limit=5, info_high_eps=1)
""", selections=[[-1]], expected_number_of_lines=3)
  out = StringIO()
  assert not is_above_limit(value=10, limit=15, eps=2, out=out)
  assert not show_diff(out.getvalue(), """\
ERROR: is_above_limit(value=10, limit=15, eps=2)
""")
  out = StringIO()
  assert not is_above_limit(value=None, limit=-3, eps=1, out=out)
  assert not is_above_limit(value=None, limit=3, eps=1, out=out)
  assert not show_diff(out.getvalue(), """\
ERROR: is_above_limit(value=None, limit=-3, eps=1)
ERROR: is_above_limit(value=None, limit=3, eps=1)
""")
  #
  import pickle, cPickle
  for p in [pickle, cPickle]:
    d = pickle_detector()
    assert d.unpickled_counter is None
    assert d.counter == 0
    s = p.dumps(d, 1)
    assert d.unpickled_counter is None
    assert d.counter == 1
    l = p.loads(s)
    assert l.unpickled_counter == 1
    assert l.counter == 0
    p.dumps(d, 1)
    assert d.counter == 2
    assert l.counter == 0
    p.dumps(l, 1)
    assert l.counter == 1
    s = p.dumps(l, 1)
    assert l.counter == 2
    k = p.loads(s)
    assert k.unpickled_counter == 2
    assert k.counter == 0
  #
  assert precision_approx_equal(0.799999,0.800004,precision=17)==True
  assert precision_approx_equal(0.799999,0.800004,precision=18)==False
  print "OK"

if (__name__ == "__main__"):
  exercise()
