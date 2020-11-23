from __future__ import absolute_import, division, print_function
from libtbx.test_utils import show_diff
from libtbx import easy_run
import os
op = os.path

def remove_file(path):
  if (op.exists(path)):
    os.remove(path)
  assert not op.exists(path)

class build_cmds_class(object):
 def __init__(self):
  # prevent race condition by naming each executable uniquely
  self.iteration = 0
 def __call__(self,tst_f, opts, ignore_ifort=False):
  self.iteration += 1
  from fable import simple_compilation
  comp_env = simple_compilation.environment()
  exe = comp_env.exe_suffix
  result = []
  if (opts.ifort and not ignore_ifort):
    remove_file("a.out")
    cmd = "ifort -diag-disable 7951 %s" % tst_f
    if (opts.verbose): print(cmd)
    easy_run.fully_buffered(command=cmd).raise_if_errors()
    assert op.exists("a.out")
    result.append("a.out")
  #remove_file("fable_cout")
  cmd = "fable.cout %s --link --exe_name=fable_cout%02d" %(tst_f,self.iteration)
  if (opts.verbose): print(cmd)
  easy_run.fully_buffered(command=cmd).raise_if_errors()
  assert op.exists("fable_cout%02d"%self.iteration+exe)
  result.append("fable_cout%02d"%self.iteration+exe)
  return [op.join(".", cmd) for cmd in result]

build_cmds = build_cmds_class()

def exercise_open(opts):
  for status in ["old", "new", "unknown", "scratch"]:
    tst_f = "exercise_open_%s.f" % status
    with open(tst_f, "w") as f:
      f.write("""\
        program prog
        open(1, file='exercise_open.tmp', status='%s', iostat=ios)
        write(6, '(l1)') (ios .eq. 0)
        end
""" % status)
    #
    for cmd in build_cmds(tst_f=tst_f, opts=opts):
      remove_file("exercise_open.tmp")
      stdout_1 = easy_run.fully_buffered(
        command=cmd).raise_if_errors().stdout_lines
      exists_1 = op.exists("exercise_open.tmp")
      #
      f = open("exercise_open.tmp", "w")
      f.close()
      assert op.exists("exercise_open.tmp")
      stdout_2 = easy_run.fully_buffered(
        command=cmd).raise_if_errors().stdout_lines
      exists_2 = op.exists("exercise_open.tmp")
      #
      with open("exercise_open.tmp", "w") as f:
        f.write("X")
      stdout_3 = easy_run.fully_buffered(
        command=cmd).raise_if_errors().stdout_lines
      exists_3 = op.exists("exercise_open.tmp")
      if (exists_3):
        size_3 = op.getsize("exercise_open.tmp")
      else:
        size_3 = None
      #
      results = [
        stdout_1, exists_1, stdout_2, exists_2, stdout_3, exists_3, size_3]
      if (opts.verbose): print("%-12s" % cmd, results)
      expected = {
        "old": [['F'], False, ['T'], True, ['T'], True, 1],
        "new": [['T'], True, ['F'], True, ['F'], True, 1],
        "unknown": [['T'], True, ['T'], True, ['T'], True, 1],
        "scratch": [['T'], False, ['T'], True, ['T'], True, 1]}
      assert results == expected[status]

def exercise_mixed_read_write(opts):
  tmp = "exercise_mixed_read_write.tmp"
  tst_f = "exercise_mixed_read_write.f"
  with open(tst_f, "w") as f:
    f.write("""\
      program prog
      open(
     &  unit=1,
     &  file='%s',
     &  status='old')
      read(1, '(i2)') num
      write(6, '(i2)') num*2
      write(1, '(i2)') 78
      end
""" % tmp)
  for cmd in build_cmds(tst_f=tst_f, opts=opts):
    if (opts.verbose): print(cmd)
    with open(tmp, "w") as f:
      f.write("""\
12
34
56
""")
    stdout = easy_run.fully_buffered(
      command=cmd).raise_if_errors().stdout_lines
    assert stdout == ["24"]
    with open(tmp, "rb") as f:
      tmp_text = f.read()
    assert not show_diff(tmp_text, b"""\
12
78
""".replace(b"\n", os.linesep.encode("latin-1")))

def exercise_read_from_non_existing_file(opts):
  tst_f = "exercise_read_from_non_existing_file.f"
  remove_file("fem_io_unit_001")
  with open(tst_f, "w") as f:
    f.write("""\
      program prog
      read(1, *) num
      end
""")
  for cmd in build_cmds(tst_f=tst_f, opts=opts, ignore_ifort=True):
    stdout = easy_run.fully_buffered(
      command=cmd, join_stdout_stderr=True).stdout_lines
    assert not show_diff(stdout, """\
std::exception what(): End of input during read
""")

def run(args):
  from libtbx.option_parser import option_parser
  command_line = (option_parser(
    usage="fable.python %s [options]" % __file__)
    .option(None, "--ifort",
      action="store_true",
      default=False)
    .option(None, "--verbose",
      action="store_true",
      default=False)
  ).process(args=args)
  keys = set(command_line.args)
  exercises = set()
  for key in globals().keys():
    if (key.startswith("exercise_")):
      exercises.add(key[9:])
  assert len(keys) == 0 or keys.issubset(exercises)
  co = command_line.options
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  if (len(keys) == 0 or "open" in keys):
    exercise_open(opts=co)
  if (len(keys) == 0 or "mixed_read_write" in keys):
    exercise_mixed_read_write(opts=co)
  if (len(keys) == 0 or "read_from_non_existing_file" in keys):
    exercise_read_from_non_existing_file(opts=co)
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
