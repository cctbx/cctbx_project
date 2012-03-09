import fable.cout

file_names_disable_warnings = set("""\
add_reals.f
add_real_integer.f
logical_a_or_b.f
add_dp_integer.f
real_array_sum.f
""".splitlines())

file_names_join_stdout_stderr = set("""\
stop_bare.f
stop_integer.f
stop_string.f
""".splitlines())

top_procedures_by_file_name = {}
for line in """\
const_analysis_1.f prog
const_analysis_2.f prog
""".splitlines():
  file_name, top_procedure = line.split()
  top_procedures_by_file_name[file_name] = [top_procedure]

dynamic_parameters_by_file_name = {
  "dynamic_parameters_1.f": [fable.cout.dynamic_parameter_props(
    name="root_size", ctype="int", default="3")],
  "dynamic_parameters_2.f": [fable.cout.dynamic_parameter_props(
    name="nums_size", ctype="int", default="2")],
  "dynamic_parameters_3.f": [fable.cout.dynamic_parameter_props(
    name="base_size", ctype="int", default="3")],
  "dynamic_parameters_4.f": [fable.cout.dynamic_parameter_props(
    name="base_size", ctype="int", default="3")],
  "dynamic_parameters_5.f": [fable.cout.dynamic_parameter_props(
    name="base_size", ctype="int", default="3")]}

common_equivalence_simple_by_file_name = {
  "common_equivalence_simple_1.f": ["info"],
  "common_equivalence_simple_2.f": ["info"],
  "common_equivalence_simple_3.f": ["tab"],
  "common_equivalence_simple_4.f": ["first"],
  "common_equivalence_simple_5.f": ["all"],
  "common_equivalence_simple_6.f": ["com"]}

def check_intrinsics_extra(text):
  import re
  lines = text.splitlines()
  def check():
    if (len(lines) != 6): return False
    if (re.match(r'\d\d-[A-Z][a-z][a-z]-\d\d', lines[0]) is None): return False
    if (re.match(r'\d\d:\d\d:\d\d', lines[1]) is None): return False
    if (len(lines[2]) != 70): return False
    if (len(lines[3]) != 6): return False
    if (lines[4] != "YkD"): return False
    if (lines[5] != "           0"): return False
    return True
  if (not check()):
    print "Unexpected output:"
    print text
    raise AssertionError

class file_name_and_expected_cout_info(object):

  __slots__ = [
    "inp_lines",
    "out_lines",
    "skip_run",
    "ifort_diff_behavior",
    "ifort_diff_floating_point_format"]

  def __init__(O):
    O.inp_lines = []
    O.out_lines = []
    O.skip_run = False
    O.ifort_diff_behavior = False
    O.ifort_diff_floating_point_format = False

def read_file_names_and_expected_cout(test_valid):
  from fable.utils import keyed_lists
  import os.path as op
  text = open(op.join(test_valid, "file_names_and_expected_cout")).read()
  result = keyed_lists()
  file_name, info = None, None
  for line in text.splitlines():
    if (line.startswith("<")):
      assert file_name is not None
      assert line.endswith("<")
      info.inp_lines.append(line[1:-1])
    elif (line.startswith("|")):
      assert file_name is not None
      assert line.endswith("|")
      info.out_lines.append(line[1:-1])
    else:
      if (file_name is not None):
        result.get(file_name).append(info)
      info = file_name_and_expected_cout_info()
      if (line.startswith("!@")):
        file_name = line[2:]
        info.skip_run = True
      elif (line.startswith("!=")):
        file_name = line[2:]
        info.ifort_diff_behavior = True
      elif (line.startswith("!%")):
        file_name = line[2:]
        info.ifort_diff_floating_point_format = True
      else:
        file_name = line
  if (file_name is not None):
    result.get(file_name).append(info)
  return result

def regex_select(keyed_lists, regex_patterns):
  result = []
  for key,list in keyed_lists.items():
    def key_matches_regex():
      if (len(regex_patterns) == 0):
        return True
      from re import search
      for pattern in regex_patterns:
        if (search(pattern, key) is not None):
          return True
      return False
    if (key_matches_regex()):
      result.append((key,list))
  return result

class process_file_info(object):

  __slots__ = ["opts", "comp_env", "test_valid"]

  def __init__(O, opts, comp_env, test_valid):
    O.opts = opts
    O.comp_env = comp_env
    O.test_valid = test_valid

  def __call__(O, file_info):
    from libtbx import easy_run
    from libtbx.str_utils import show_string
    from libtbx.test_utils import show_diff
    from cStringIO import StringIO
    import os.path as op
    import sys
    opts = O.opts
    file_name, io_infos = file_info
    if (opts.verbose):
      print file_name
    file_path = op.join(O.test_valid, file_name)
    top_procedures = top_procedures_by_file_name.get(file_name)
    common_equivalence_simple_list = [set(
      common_equivalence_simple_by_file_name.get(file_name, []))]
    if (len(common_equivalence_simple_list[0]) != 0):
      common_equivalence_simple_list.append([])
    for i_ces,common_equivalence_simple in \
          enumerate(common_equivalence_simple_list):
      common_report_stringio = StringIO()
      try:
        lines = fable.cout.process(
          file_names=[file_path],
          top_procedures=top_procedures,
          dynamic_parameters=dynamic_parameters_by_file_name.get(file_name),
          common_equivalence_simple=common_equivalence_simple,
          common_report_stringio=common_report_stringio)
      except Exception:
        if (not opts.keep_going): raise
        print "\nEXCEPTION: fable.cout.process([%s])\n" % file_name
        return 1
      have_simple_equivalence = (
        "\n".join(lines).find(" // SIMPLE EQUIVALENCE") >= 0)
      if (len(common_equivalence_simple) != 0):
        assert have_simple_equivalence
      else:
        assert not have_simple_equivalence
      assert file_name.endswith(".f")
      base_name = file_name[:-2]
      if (len(common_equivalence_simple_list) != 1):
        base_name += "_alt%d" % i_ces
      fem_cpp = base_name + "_fem.cpp"
      fem_exe_name = fem_cpp[:-4] + O.comp_env.exe_suffix
      print >> open(fem_cpp, "w"), "\n".join(lines)
      if (opts.ifort):
        ifort_exe_name = base_name + "_ifort"
        ifort_cmd = "ifort -diag-disable 7000 -o %s %s" % (
          ifort_exe_name, show_string(file_path))
      else:
        ifort_exe_name = None
        ifort_cmd = None
      if (opts.dry_run):
        return 0
      #
      n_failures = [0]
      def handle_exception(e):
        n_failures[0] += 1
        if (not opts.keep_going): raise
        print
        print str(e)
        print
        sys.stdout.flush()
      #
      class BuildError(RuntimeError): pass
      try:
        O.comp_env.build(
          link=True,
          file_name_cpp=fem_cpp,
          exe_name=fem_exe_name,
          disable_warnings=(file_name in file_names_disable_warnings),
          show_command=opts.verbose,
          Error=BuildError)
      except BuildError, e:
        handle_exception(e)
        fem_exe_name = None
      #
      if (ifort_cmd is not None):
        if (opts.verbose):
          print ifort_cmd
        buffers = easy_run.fully_buffered(command=ifort_cmd)
        try:
          buffers.raise_if_errors_or_output(Error=BuildError)
        except BuildError, e:
          handle_exception(e)
          ifort_exe_name = None
      #
      for info in io_infos:
        if (info.skip_run):
          if (opts.verbose):
            print "Skipping run:", file_name
          continue
        if (len(info.inp_lines) != 0 and opts.verbose):
          print "  number of input lines:", len(info.inp_lines)
        sys.stdout.flush()
        for exe_name in [fem_exe_name, ifort_exe_name]:
          if (exe_name is None): continue
          cmd = cmd0 = op.join(".", exe_name)
          if (opts.valgrind):
            cmd = "valgrind " + cmd
          if (opts.verbose):
            print cmd
            sys.stdout.flush()
          join_stdout_stderr = (
               opts.valgrind
            or (file_name in file_names_join_stdout_stderr))
          buffers = easy_run.fully_buffered(
            command=cmd,
            stdin_lines=info.inp_lines,
            join_stdout_stderr=join_stdout_stderr)
          if (not join_stdout_stderr):
            class ExeError(RuntimeError): pass
            try:
              buffers.raise_if_errors(Error=ExeError)
            except ExeError, e:
              handle_exception(e)
              buffers = None
          if (buffers is not None):
            text = "\n".join(buffers.stdout_lines)
            if (opts.valgrind):
              print text
            else:
              def check(text):
                if (file_name == "intrinsics_extra.f"):
                  check_intrinsics_extra(text)
                  return
                if (file_name == "sf.f"):
                  text = text.replace(" -0.620088", " -0.620087")
                elif (file_name == "unformatted_experiments.f"):
                  if (sys.byteorder == "big"):
                    text = text \
                      .replace(
                        "        1234        5678",
                        "        5678        1234") \
                      .replace(
                        "        18558553691448",
                        "        23330262356193")
                have_diffs = show_diff(text, "\n".join(info.out_lines))
                def assert_not_have_diffs():
                  if (opts.keep_going):
                    print "WARNING: --keep-going after show_diff:", exe_name
                  else:
                    assert not have_diffs
                if (have_diffs):
                  if (exe_name is fem_exe_name):
                    assert_not_have_diffs()
                  elif (exe_name is ifort_exe_name):
                    if (    not info.ifort_diff_behavior
                        and not info.ifort_diff_floating_point_format):
                      assert_not_have_diffs()
                  else:
                    raise AssertionError
              check(text)
          def run_with_args(args):
            cmda = cmd0 + " " + args
            if (opts.verbose):
              print cmda
              sys.stdout.flush()
            result = easy_run.fully_buffered(
              command=cmda, join_stdout_stderr=True)
            if (opts.valgrind):
              cmda = "valgrind " + cmda
              if (opts.verbose):
                print cmda
                sys.stdout.flush()
              buffers = easy_run.fully_buffered(
                command=cmda, join_stdout_stderr=True)
              print "\n".join(buffers.stdout_lines)
            return result
          if (file_name == "read_lines.f"):
            exercise_end_of_line(exe_name=exe_name, verbose=opts.verbose)
          elif (file_name == "dynamic_parameters_1.f"):
            buffers = run_with_args("--fem-dynamic-parameters=5")
            assert not show_diff(buffers.stdout_lines, """\
          14          15          16          17          18          19
          20          21          22          23
""")
            buffers = run_with_args("--fem-dynamic-parameters=5,6")
            assert buffers.stdout_lines[0].endswith(
              "Too many --fem-dynamic-parameters fields"
              " (given: 2, max. expected: 1)")
            buffers = run_with_args("--fem-dynamic-parameters=x")
            assert buffers.stdout_lines[0].endswith(
              'Invalid --fem-dynamic-parameters field (field 1): "x"')
          elif (file_name == "intrinsics_iargc_getarg.f"):
            buffers = run_with_args("D rP uWq")
            assert not show_diff(buffers.stdout_lines, "\n".join([
              "A", "D   ", "rP  ", "uWq ",
              "B", "uWq ", "rP  ", "D   ",
              "C", "rP  ", "uWq ", "D   "]) + "\n")
    #
    return n_failures[0]

def exercise_end_of_line(exe_name, verbose):
  lines = """\
a
bc
def
ghij
klmno
""".splitlines()
  open("unix.txt", "wb").write("\n".join(lines)+"\n")
  open("dos.txt", "wb").write("\r\n".join(lines)+"\r\n")
  open("dos2.txt", "wb").write("\r\r\n".join(lines)+"\r\n")
  open("mac.txt", "wb").write("\r".join(lines)+"\r")
  from libtbx import easy_run
  from libtbx.utils import remove_files
  import os
  op = os.path
  expected_outputs = [
    "a   \nbc  \ndef \nghij\nklmn\n",
    "a   \nbc  \ndef \nghij\nklmn\n",
    "a\r  \nbc\r \ndef\r\nghij\nklmn\n",
    "a\rbc\n"]
  for vers,expected in zip(["unix", "dos", "dos2", "mac"], expected_outputs):
    remove_files(paths=["read_lines_out"])
    cmd = "%s < %s.txt > read_lines_out" % (op.join(".", exe_name), vers)
    if (verbose): print cmd
    easy_run.fully_buffered(command=cmd).raise_if_errors_or_output()
    assert op.isfile("read_lines_out")
    result = open("read_lines_out", "rb").read()
    assert result == expected

def exercise_compile_valid(regex_patterns, opts):
  from fable import cout
  from fable import simple_compilation
  comp_env = simple_compilation.environment()
  if (comp_env.compiler_path is None):
    print "Skipping exercise_compile_valid(): %s not available." % \
      comp_env.compiler
    return
  import libtbx.load_env
  import os.path as op
  fable_dist = libtbx.env.dist_path(module_name="fable")
  test_valid = op.join(fable_dist, "test/valid")
  selected_file_names_and_expected_cout = regex_select(
    keyed_lists=read_file_names_and_expected_cout(test_valid=test_valid),
    regex_patterns=regex_patterns)
  assert len(selected_file_names_and_expected_cout) != 0
  #
  if (opts.pch):
    comp_env.build(
      link=False,
      file_name_cpp=op.join(fable_dist, "fem.hpp"),
      pch_name="fem.hpp",
      show_command=True)
    comp_env.set_have_pch()
    print
  #
  processor = process_file_info(
    opts=opts, comp_env=comp_env, test_valid=test_valid)
  #
  if (not opts.multiprocessing):
    n_failures = 0
    for file_info in selected_file_names_and_expected_cout:
      n_failures += processor(file_info=file_info)
    return n_failures
  #
  n_proc = min(
    len(selected_file_names_and_expected_cout),
    opts.max_proc)
  print "Number of processors:", n_proc
  import multiprocessing
  mp_pool = multiprocessing.Pool(processes=n_proc)
  return sum(mp_pool.map(processor, selected_file_names_and_expected_cout))

def run(args):
  from libtbx.option_parser import option_parser
  command_line = (option_parser(
    usage="fable.python %s [options] regex_pattern ..." % __file__)
    .enable_multiprocessing()
    .option(None, "--dry_run",
      action="store_true",
      default=False)
    .option(None, "--valgrind",
      action="store_true",
      default=False)
    .option(None, "--ifort",
      action="store_true",
      default=False)
    .option(None, "--keep_going",
      action="store_true",
      default=False)
    .option(None, "--pch",
      action="store_true",
      default=False)
    .option(None, "--verbose",
      action="store_true",
      default=False)
  ).process(args=args)
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  n_failures = exercise_compile_valid(
    regex_patterns=command_line.args,
    opts=command_line.options)
  if (n_failures != 0):
    print "Done."
  else:
    print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
