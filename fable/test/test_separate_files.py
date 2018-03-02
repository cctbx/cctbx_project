from __future__ import absolute_import, division, print_function

import os
import pytest
import fable.cout
import fable.simple_compilation
from libtbx import easy_run

@pytest.mark.parametrize("""
    file_names_cpp,
    number_of_function_files,
    separate_files_main_namespace,
    separate_files_separate_namespace
    """, [
      (["top.cpp", "functions.cpp"], 1, {}, {}),
      (["top.cpp", "subs.cpp"], None, {}, {"subs": ["sub1", "sub2"]}),
      (["top.cpp", "subs.cpp", "functions.cpp"], 1, {}, {"subs": ["sub1", "sub2"]}),
      (["top.cpp", "subs.cpp"], None, {"subs": ["sub1", "sub2"]}, {}),
      (["top.cpp", "subs.cpp", "functions.cpp"], 1, {"subs": ["sub1", "sub2"]}, {})
])
def test_exercise(
    tmpdir, testsdir,
    file_names_cpp,
    number_of_function_files,
    separate_files_main_namespace,
    separate_files_separate_namespace):
  tmpdir.chdir()
  test_valid = os.path.join(testsdir, 'valid')
  top_cpp = fable.cout.process(
    file_names=[os.path.join(test_valid, "subroutine_3.f")],
    top_procedures=["prog"],
    namespace="tst_separate_files",
    top_cpp_file_name=file_names_cpp[0],
    number_of_function_files=number_of_function_files,
    separate_files_main_namespace=separate_files_main_namespace,
    separate_files_separate_namespace=separate_files_separate_namespace)
  comp_env = fable.simple_compilation.environment()
  file_names_obj = []
  for file_name_cpp in file_names_cpp:
    obj = comp_env.file_name_obj(file_name_cpp=file_name_cpp)
    assert not os.path.exists(obj)
    cmd = comp_env.compilation_command(file_name_cpp=file_name_cpp)
    print(cmd)
    easy_run.call(command=cmd)
    assert os.path.exists(obj)
    file_names_obj.append(obj)
  exe_root = "tst_separate_files"
  exe = comp_env.file_name_exe(exe_root=exe_root)
  assert not os.path.exists(exe)
  cmd = comp_env.link_command(file_names_obj=file_names_obj, exe_root=exe_root)
  print(cmd)
  easy_run.call(command=cmd)
  cmd = os.path.join(".", exe)
  print(cmd)
  assert os.path.exists(cmd)
  stdout = easy_run.fully_buffered(command=cmd).raise_if_errors().stdout_lines
  from fable.tst_cout_compile import read_file_names_and_expected_cout
  info = read_file_names_and_expected_cout(test_valid=test_valid).get("subroutine_3.f")[0]
  assert stdout == info.out_lines
