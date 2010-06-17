import os
op = os.path

def remove_file_if_necessary(file_name):
  if (op.isfile(file_name)): os.remove(file_name)
  if (op.exists(file_name)):
    from libtbx.str_utils import show_string
    raise RuntimeError(
      "Unable to remove file: %s" % show_string(file_name))

def exercise(
      verbose,
      file_names_cpp,
      number_of_function_files=None,
      separate_files_main_namespace={},
      separate_files_separate_namespace={}):
  if (verbose): print "next exercise"
  import libtbx.load_env
  test_valid = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  import fable.cout
  top_cpp = fable.cout.process(
    file_names=[op.join(test_valid, "subroutine_3.f")],
    top_unit_name="prog",
    namespace="tst_separate_files",
    top_cpp_file_name=file_names_cpp[0],
    number_of_function_files=number_of_function_files,
    separate_files_main_namespace=separate_files_main_namespace,
    separate_files_separate_namespace=separate_files_separate_namespace)
  from fable import simple_compilation
  comp_env = simple_compilation.environment()
  from libtbx import easy_run
  file_names_obj = []
  for file_name_cpp in file_names_cpp:
    obj = comp_env.file_name_obj(file_name_cpp=file_name_cpp)
    remove_file_if_necessary(file_name=obj)
    cmd = comp_env.compilation_command(file_name_cpp=file_name_cpp)
    if (verbose): print cmd
    easy_run.call(command=cmd)
    assert op.exists(obj)
    file_names_obj.append(obj)
  exe_root = "tst_separate_files"
  exe = comp_env.file_name_exe(exe_root=exe_root)
  remove_file_if_necessary(file_name=exe)
  cmd = comp_env.link_command(file_names_obj=file_names_obj, exe_root=exe_root)
  if (verbose): print cmd
  easy_run.call(command=cmd)
  cmd = op.join(".", exe)
  if (verbose): print cmd
  assert op.exists(cmd)
  stdout = easy_run.fully_buffered(command=cmd).raise_if_errors().stdout_lines
  text = "\n".join(stdout)
  if (verbose):
    print text
  from fable.tst_cout_compile import read_file_names_and_expected_cout
  info = read_file_names_and_expected_cout(test_valid=test_valid).get(
    "subroutine_3.f")[0]
  from libtbx.test_utils import show_diff
  assert not show_diff(text, "\n".join(info.out_lines))
  if (verbose): print

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = (args == ["--verbose"])
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  all = True
  if (0 or all):
    exercise(verbose,
      file_names_cpp=["top.cpp", "functions.cpp"],
      number_of_function_files=1)
  if (0 or all):
    exercise(verbose,
      file_names_cpp=["top.cpp", "subs.cpp"],
      separate_files_separate_namespace={"subs": ["sub1", "sub2"]})
  if (0 or all):
    exercise(verbose,
      file_names_cpp=["top.cpp", "subs.cpp", "functions.cpp"],
      number_of_function_files=1,
      separate_files_separate_namespace={"subs": ["sub1", "sub2"]})
  if (0 or all):
    exercise(verbose,
      file_names_cpp=["top.cpp", "subs.cpp"],
      separate_files_main_namespace={"subs": ["sub1", "sub2"]})
  if (0 or all):
    exercise(verbose,
      file_names_cpp=["top.cpp", "subs.cpp", "functions.cpp"],
      number_of_function_files=1,
      separate_files_main_namespace={"subs": ["sub1", "sub2"]})
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
