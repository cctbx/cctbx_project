from __future__ import absolute_import, division, print_function

import os

from libtbx import easy_run
import pytest

command_expectation_list = [
    ('fable.split "%ssubroutine_1.f"', 'program_prog.f'),
    ('fable.read --each "%swrite_star.f"', 'Success: 1'),
    ('fable.read --warnings "%sequivalence_mixed.f"',
      'Warning: EQUIVALENCE cluster with mixed data types: integer, real:'),
    ('fable.show_calls --write-graphviz-dot=tmp.dot'
      ' --top-procedure=sub "%sexternal_arg_layers.f"',
        'exch->exch_imp'),
    ('fable.show_calls "%sdependency_cycle.f"', 'sub1 sub2'),
    ('fable.fem_include_search_paths --with-quotes', 'fable'),
    ('fable.cout --each "%swrite_star.f"', 'return fem::main_with_catch'),
    ('fable.cout "%scommon_variants.f"',
      'Writing file: "fable_cout_common_report"'),
    ('fable.cout "%ssubroutine_3.f" --top-procedure=sub3'
     ' --fortran-file-comments',
      'nums(i) = i * 20;'),
    ('fable.cout --compile "%swrite_star.f"',
      'placeholder_please_replace::program_prog);'),
    ('fable.cout --namespace=test --run "%swrite_star.f"',
      'test::program_prog);'),
    ('fable.cout %ssf.f --namespace example --run', '  -3   1  -5'),
    ('fable.cout "%sdynamic_parameters_1.f"'
      ' --dynamic-parameter="int root_size=1"',
        "const int root_size = cmn.dynamic_params.root_size;"),
]

@pytest.mark.parametrize('test_number', range(len(command_expectation_list)))
def test_run_fable_command(tmpdir, testsdir, test_number):
  command, expected_output_fragment = command_expectation_list[test_number]
  join_stdout_stderr = 'fable_cout_common_report' in expected_output_fragment
  if "%s" in command:
    command = command % (os.path.join(testsdir, 'valid') + os.sep)
  print(command)
  with tmpdir.as_cwd():
    run_buffers = easy_run.fully_buffered(
      command=command,
      join_stdout_stderr=join_stdout_stderr)
  if not join_stdout_stderr:
    run_buffers.raise_if_errors()
  stdout_buffer = "\n".join(run_buffers.stdout_lines)
  if expected_output_fragment:
    assert expected_output_fragment in stdout_buffer, \
           "Expected fragment missing from output."
  else:
    assert stdout_buffer == '', \
           "Expected no output from command."
