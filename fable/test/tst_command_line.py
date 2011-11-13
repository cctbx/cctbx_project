def run(args):
  assert len(args) == 0
  from libtbx import easy_run
  import libtbx.load_env
  import os
  op = os.path
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  assert t_dir.find('"') < 0
  n_errors = 0
  for command,expected_output_fragment in [
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
         ' --fortran-file_comments',
          'nums(i) = i * 20;'),
        ('fable.cout --compile "%swrite_star.f"',
          'placeholder_please_replace::program_prog);'),
        ('fable.cout --namespace=test --run "%swrite_star.f"',
          'test::program_prog);'),
        ('fable.cout --example', '  -3   1  -5'),
        ('fable.cout "%sdynamic_parameters_1.f"'
          ' --dynamic-parameter="int root_size=1"',
            "const int root_size = cmn.dynamic_params.root_size;")]:
    if (expected_output_fragment.find("fable_cout_common_report") >= 0):
      join_stdout_stderr = True
    else:
      join_stdout_stderr = False
    if (command.find("%s") >= 0):
      command = command % (t_dir + os.sep)
    print command
    run_buffers = easy_run.fully_buffered(
      command=command,
      join_stdout_stderr=join_stdout_stderr)
    class SpecificError(RuntimeError): pass
    try:
      if (not join_stdout_stderr):
        run_buffers.raise_if_errors(Error=SpecificError)
    except SpecificError, e:
      n_errors += 1
      print "ERROR:"
      print str(e)
      print
    else:
      stdout_buffer = "\n".join(run_buffers.stdout_lines)
      if (expected_output_fragment is None):
        if (len(stdout_buffer) != 0):
          n_errors += 1
          print stdout_buffer
          print "ERROR: unexpected output above."
          print
      elif (stdout_buffer.find(expected_output_fragment) < 0):
        n_errors += 1
        print stdout_buffer
        print "ERROR: not found in output above:"
        print [expected_output_fragment]
        print
  if (n_errors != 0):
    print "Number of errors:", n_errors
    print "Done."
  else:
    print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
