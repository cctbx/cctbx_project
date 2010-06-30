from fable import read
from libtbx.test_utils import Exception_expected, show_diff
import libtbx.load_env
import os
op = os.path

def exercise_strip_spaces_separate_strings():
  from fable.read import Error, source_line, strip_spaces_separate_strings
  import itertools
  global_line_index_generator = itertools.count()
  def slc(cmbnd):
    return [
      source_line(
        global_line_index_generator=global_line_index_generator,
        file_name="str",
        line_number=i+1,
        text="      "+line)
          for i,line in enumerate(cmbnd.splitlines())]
  def check(cmbnd, expected_code, expected_strings):
    for cmbnd_work,expected_strings_work in [
          (cmbnd, expected_strings),
          (cmbnd.replace("'", '"'), [s.replace("'", '"')
            for s in expected_strings])]:
      ssl = strip_spaces_separate_strings(source_line_cluster=slc(cmbnd_work))
      assert not show_diff(ssl.code, expected_code)
      assert len(ssl.strings) == len(ssl.string_indices)
      assert ssl.strings == expected_strings_work
      expected_string_indices = []
      for i,c in enumerate(expected_code):
        if (c == "'" or c == '"'):
          expected_string_indices.append(i)
      assert ssl.string_indices == expected_string_indices
      if (cmbnd.find('"') > 0):
        break
  check("a = 0", "a=0", [])
  check("a = ''", "a='", [""])
  check("a = '\"'", "a='", ["\""])
  check("a = \"'\"", "a='", ["'"])
  check("a = 'x'", "a='", ["x"])
  check("a = ' x ' ", "a='", [" x "])
  check("a = ''''", "a='", ["'"])
  check("call foo('')", "callfoo(')", [""])
  check("call foo('''')", "callfoo(')", ["'"])
  check("c a l l f o   o ( ''''  )  ", "callfoo(')", ["'"])
  check("   c   a   l   l   f   o   o ( ''''  )  ", "callfoo(')", ["'"])
  check("   C  A   L  LF   O O ( ' abc '''' def '''  , ' g''h''i''' , X )  ",
    "callfoo(',',x)", [" abc '' def '", " g'h'i'"])
  check("a = '\n'", "a='", [""])
  check("a = 'x\n'", "a='", ["x"])
  check("a = '\ny'", "a='", ["y"])
  check("a = 'x\ny'", "a='", ["xy"])
  check("a = '''\n'", "a='", ["'"])
  check("a = '\n'''", "a='", ["'"])
  check("a = '''\n'''", "a='", ["''"])
  #
  for cmbnd,q,nd in [("'abc", "'", 9), ('x="', '"', 11)]:
    try:
      strip_spaces_separate_strings(source_line_cluster=slc(cmbnd))
    except Error, e:
      assert not show_diff(str(e), """\
Missing terminating %s character:
  at str(1):
  |      %s|
%s^""" % (q, cmbnd, "-"*nd))
    else: raise Exception_expected

def exercise_valid(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  #
  read_already = set()
  def get_units(file_name):
    if (verbose):
      print "exercise_valid:", file_name
    read_already.add(file_name)
    return read.process(file_names=[op.join(t_dir, file_name)])
  #
  def get_program(file_name):
    units = get_units(file_name)
    assert len(units.program) == 1
    assert len(units.subroutine) == 0
    assert len(units.function) == 0
    assert len(units.blockdata) == 0
    return units.program[0]
  #
  prog = get_program("goto_forms.f")
  keys = [ei.key for ei in prog.executable]
  assert not show_diff("\n".join(keys), """\
goto_computed
goto_computed
goto_computed
goto_computed
continue
continue""")
  #
  get_program("string_spanning_continuation_lines.f")
  #
  units = get_units("sf.f")
  assert len(units.program) == 1
  assert len(units.subroutine) == 1
  #
  for file_name in sorted(os.listdir(t_dir)):
    if (file_name.endswith(".f") and not file_name in read_already):
      get_units(file_name)

def exercise_lenient(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/lenient", test=op.isdir)
  #
  def get(file_name):
    if (verbose):
      print "exercise_lenient:", file_name
    return read.process(file_names=[op.join(t_dir, file_name)])
  #
  get("str_blank_str.f")
  get("str_cont_line_str.f")

def exercise_syntax_error(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/syntax_error", test=op.isdir)
  def fail(file_name):
    if (verbose):
      print "exercise_syntax_error:", file_name
    read.process(file_names=[op.join(t_dir, file_name)])
  from fable.read import Error
  try:
    fail("label_cont_char.f")
  except Error, e:
    assert str(e).startswith(
      "A continuation character is illegal on a line with a statement label:")
    assert str(e).endswith("""\
  |    3xk=6|
--------^""")
  else: raise Exception_expected
  try:
    fail("label_empty.f")
  except Error, e:
    assert str(e).startswith("Labelled statement is empty:")
    assert str(e).endswith("""\
  |    1  |
---------^""")
  else: raise Exception_expected
  try:
    fail("not_an_identifier.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      external a, x%z|
----------------------^""")
  else: raise Exception_expected
  try:
    fail("bare_integer_i.f")
  except Error, e:
    assert str(e).startswith("Missing END for PROGRAM:")
    assert str(e).endswith("  |      integer i|")
  else: raise Exception_expected
  try:
    fail("closing_parenthesis_without_matching.f")
  except Error, e:
    assert str(e).startswith(
      'Closing ")" without a matching opening parenthesis:')
    assert str(e).endswith("""\
  |      a = 0)|
--------------^""")
  else: raise Exception_expected
  try:
    fail("opening_parenthesis_without_matching.f")
  except Error, e:
    assert str(e).startswith(
      'Missing a closing ")":')
    assert str(e).endswith("""\
  |      a = (0|
-------------^""")
  else: raise Exception_expected
  try:
    fail("dot_e.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |        x       =       .e0|
----------------------------^""")
  else: raise Exception_expected
  try:
    fail("x_assign_dot.f")
  except Error, e:
    assert str(e).startswith("Expression unexpectedly ends with a dot:")
    assert str(e).endswith("""\
  |      x = .|
-------------^""")
  try:
    fail("x_assign_dot_2.f")
  except Error, e:
    assert str(e).startswith("Expression unexpectedly ends with a dot:")
    assert str(e).endswith("""\
  |        x = .|
---------------^""")
  else: raise Exception_expected
  try:
    fail("x_assign_1ex.f")
  except Error, e:
    assert str(e).startswith("Invalid floating-point literal:")
    assert str(e).endswith("""\
  |      x = 1ex|
---------------^""")
  else: raise Exception_expected
  try:
    fail("x_assign_1dotd.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      x = 1.d|
---------------^""")
  else: raise Exception_expected
  try:
    fail("bad_false.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = .fals.|
--------------^""")
  else: raise Exception_expected
  try:
    fail("bad_true.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = .true|
--------------^""")
  else: raise Exception_expected
  try:
    fail("bad_not.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = .not b|
--------------^""")
  else: raise Exception_expected
  try:
    fail("bad_not_2.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = 1 .not. b|
----------------^""")
  else: raise Exception_expected
  try:
    fail("bad_gt.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = 1 .gt 2|
----------------^""")
  else: raise Exception_expected
  try:
    fail("bad_after_dot.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = 1 .2et. 2|
----------------^""")
  else: raise Exception_expected
  try:
    fail("j_assign_i_percent_5.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      j = i % 5|
---------------^""")
  else: raise Exception_expected
  try:
    fail("bad_and.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = b .ad. c|
-----------------^""")
  else: raise Exception_expected
  try:
    fail("bad_or.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      a = b .or c|
-----------------^""")
  else: raise Exception_expected
  try:
    fail("exclamation_mark_syndrome.f")
  except Error, e:
    assert str(e).startswith("Missing terminating ' character:")
    assert str(e).endswith("""\
  |     !ef'|
-----------^""")
  else: raise Exception_expected
  try:
    fail("common_with_data_size.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      common /com/ vals(2), nums*4(2)|
-----------------------------------^""")
  else: raise Exception_expected
  try:
    fail("dimension_with_data_size.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("  |      dimension strings*4(2)|")
  else: raise Exception_expected
  try:
    fail("save_with_dims.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("  |      save nums(2)|")
  else: raise Exception_expected
  try:
    fail("sub_no_name.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      subroutine|
------------------^""")
  else: raise Exception_expected
  try:
    fail("sub_percent_3.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      subroutine sub % 3|
------------------------^""")
  else: raise Exception_expected
  try:
    fail("sub_open_parenthesis.f")
  except Error, e:
    assert str(e).startswith('Missing a closing ")":')
    assert str(e).endswith("""\
  |      subroutine sub(|
-----------------------^""")
  else: raise Exception_expected
  try:
    fail("sub_bad_comma.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      subroutine sub(a,)|
--------------------------^""")
  else: raise Exception_expected
  try:
    fail("fun_star.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      function fun(a,*)|
------------------------^""")
  else: raise Exception_expected
  try:
    fail("sub_bad_trailing.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      subroutine sub(a,b) x|
-----------------------------^""")
  else: raise Exception_expected
  try:
    fail("save_bad_comma.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      save num,|
-----------------^""")
  else: raise Exception_expected
  try:
    fail("save_num_comma_colon.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      save num, :|
-------------------^""")
  else: raise Exception_expected
  try:
    fail("save_num_val_colon.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      save num, val :|
-----------------------^""")
  else: raise Exception_expected
  try:
    fail("bare_external.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      external|
----------------^""")
  else: raise Exception_expected
  try:
    fail("save_slash_slash.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      save //|
--------------^""")
  else: raise Exception_expected
  try:
    fail("bare_data.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      data|""")
  else: raise Exception_expected
  try:
    fail("data_plus_repetition.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      data nums /+2*3/|
----------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_format_1.f")
  except Error, e:
    assert str(e).startswith('Format string must start with "("')
    assert str(e).endswith("""\
  |      write(6, '') num|
------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_format_2.f")
  except Error, e:
    assert str(e).startswith('Format string must end with ")"')
    assert str(e).endswith("""\
  |      write(6, '(') num|
------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_format_3.f")
  except Error, e:
    assert str(e).startswith("""\
Missing terminating ' within character format specifier "(')":""")
    assert str(e).endswith("""\
  |      write(6, '('')')|
--------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_format_4.f")
  except Error, e:
    assert str(e).startswith("Invalid FORMAT specification:")
    assert str(e).endswith("""\
  |      write(6, '(+2x)')|
--------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_format_5.f")
  except Error, e:
    assert str(e).startswith("Invalid FORMAT specification:")
    assert str(e).endswith("""\
  |      write(6, '(i2.)')|
----------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_format_6.f")
  except Error, e:
    assert str(e).startswith("Invalid FORMAT specification:")
    assert str(e).endswith("""\
  |      write(6, '(tx)')|
---------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_format_7.f")
  except Error, e:
    assert str(e).startswith("Invalid FORMAT specification:")
    assert str(e).endswith("""\
  |      write(6, '(tl)')|
---------------------^""")
  else: raise Exception_expected
  try:
    fail("format_without_label.f")
  except Error, e:
    assert str(e).startswith(
      "FORMAT without a statement label in columns 1-5:")
    assert str(e).endswith("""\
  |        format(1x)|
-----------^""")
  try:
    fail("duplicate_format_label.f")
  except Error, e:
    assert str(e).startswith(
      "Duplicate statement label in columns 1-5:")
    assert str(e).endswith("""\
  |   10 format(i2)|
---^""")
  else: raise Exception_expected
  try:
    fail("bad_implied_do_1.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      write(6, *) (i,i=1)|
---------------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_implied_do_2.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      write(6, *) (i=1,2)|
-----------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_implied_do_3.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      write(6, *) i,i=1,2|
------------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_implied_do_4.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      write(6, *) (i,i=1,j=2)|
-----------------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_implied_do_5.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      write(6, *) (i,0=1,2)|
-------------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_implied_do_6.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      write(6, *) (i,i+j=1,2)|
---------------------------^""")
  else: raise Exception_expected
  try:
    fail("bad_data.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      data (nums(i),i) /1,2/|
------------------------^""")
  else: raise Exception_expected
  try:
    fail("read_star_comma_empty.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      read *,|
---------------^""")
  else: raise Exception_expected
  try:
    fail("read_star_name.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      read * name|
----------------^""")
  else: raise Exception_expected
  try:
    fail("read_plus_name.f")
  except Error, e:
    assert str(e).startswith("Syntax error:")
    assert str(e).endswith("""\
  |      read + name|
--------------^""")
  else: raise Exception_expected

def exercise_semantic_error(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/semantic_error", test=op.isdir)
  from fable import SemanticError
  def fail(file_name):
    if (verbose):
      print "exercise_semantic_error:", file_name
    read.process(file_names=[op.join(t_dir, file_name)])
  try:
    fail("missing_include.f")
  except SemanticError, e:
    assert str(e).startswith("Missing include file:")
    assert str(e).endswith("""\
  |      include '/bin/sh/should/never/exist'|
-----------------^""")
  else: raise Exception_expected
  try:
    fail("real_declared_twice.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: val:")
    assert str(e).endswith("""\
  |      real val|
--------------^""")
  else: raise Exception_expected
  try:
    fail("external_array.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: f2:")
    assert str(e).endswith("""\
  |      external f2|
------------------^""")
  else: raise Exception_expected
  try:
    fail("dimension_unknown_data_type.f")
  except SemanticError, e:
    assert str(e).startswith("Unknown data type: nums:")
    assert str(e).endswith("""\
  |      write(6, *) nums(1)|
---------------------^""")
  else: raise Exception_expected
  try:
    fail("dims_repeated_in_dimension.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated dimension: num:")
    assert str(e).endswith("""\
  |      dimension num(2)|
-------------------^""")
  else: raise Exception_expected
  try:
    fail("dims_repeated_in_common.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated dimension: num:")
    assert str(e).endswith("""\
  |      common /com/ num(2)|
----------------------^""")
  else: raise Exception_expected
  try:
    fail("parameter_array.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: nums:")
    assert str(e).endswith("""\
  |      integer nums(2)|
-----------------^""")
  else: raise Exception_expected
  try:
    fail("parameter_in_common.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: num:")
    assert str(e).endswith("""\
  |      common /com/ num|
----------------------^""")
  else: raise Exception_expected
  try:
    fail("parameter_save.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: num:")
    assert str(e).endswith("""\
  |      save num|
--------------^""")
  else: raise Exception_expected
  try:
    fail("parameter_in_sub_args.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: num:")
    assert str(e).endswith("""\
  |      parameter(num=0)|
-------------------^""")
  else: raise Exception_expected
  try:
    fail("twice_in_sub_args.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: val:")
    assert str(e).endswith("""\
  |      subroutine sub(val, val)|
-----------------------------^""")
  else: raise Exception_expected
  try:
    fail("subroutine_name_is_also_arg.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: sub:")
    assert str(e).endswith("""\
  |      subroutine sub(sub)|
------------------------^""")
  else: raise Exception_expected
  try:
    fail("function_name_is_also_arg.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: fun:")
    assert str(e).endswith("""\
  |      function fun(fun)|
----------------------^""")
  else: raise Exception_expected
  try:
    fail("unknown_intrinsic.f")
  except SemanticError, e:
    assert str(e).startswith("Unknown intrinsic: unk:")
    assert str(e).endswith("""\
  |      write(6, *) unk(0)|
---------------------^""")
  else: raise Exception_expected
  try:
    fail("intrinsic_common.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: abs:")
    assert str(e).endswith("""\
  |      intrinsic abs|
-------------------^""")
  else: raise Exception_expected
  try:
    fail("external_common.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: nums:")
    assert str(e).endswith("""\
  |      external nums|
------------------^""")
  else: raise Exception_expected
  try:
    fail("calling_array.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting declaration: nums:")
    assert str(e).endswith("""\
  |      call nums(3)|
--------------^""")
  else: raise Exception_expected
  try:
    fail("calling_dimension.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting declaration: nums:")
    assert str(e).endswith("""\
  |      call nums(2)|
--------------^""")
  else: raise Exception_expected
  try:
    fail("function_data_type_decl_twice.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: fun:")
    assert str(e).endswith("""\
  |      integer fun|
-----------------^""")
  else: raise Exception_expected
  try:
    fail("intrinsic_dimension.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting or repeated declaration: nint:")
    assert str(e).endswith("""\
  |      intrinsic nint|
-------------------^""")
  else: raise Exception_expected
  try:
    fail("sub_fun_2.f")
  except SemanticError, e:
    assert str(e).startswith("Conflicting declaration: fun:")
    assert str(e).endswith("""\
  |      y = fun(x)|
-------------^""")
  else: raise Exception_expected
  try:
    fail("write_without_unit.f")
  except SemanticError, e:
    assert str(e).startswith("Required UNIT information is not defined:")
    assert str(e).endswith("""\
  |      write(fmt='(i3)') num|
--------------^""")
  else: raise Exception_expected
  try:
    fail("write_end.f")
  except SemanticError, e:
    assert str(e).startswith("END is invalid for WRITE statements:")
    assert str(e).endswith("""\
  |      write(10, end=20) num|
-----------------------^""")
  else: raise Exception_expected
  try:
    fail("equivalence_external.f")
  except SemanticError, e:
    assert str(e).startswith("Invalid EQUIVALENCE:")
    assert str(e).endswith("""\
  |      equivalence (ne, nl)|
----------------------^""")
  else: raise Exception_expected

def exercise_unsupported(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/unsupported", test=op.isdir)
  def fail(file_name):
    if (verbose):
      print "exercise_unsupported:", file_name
    read.process(file_names=[op.join(t_dir, file_name)])
  try:
    fail("hollerith_cont_lines.f")
  except RuntimeError, e:
    assert str(e).startswith(
      "FATAL: Not supported:"
      " FORMAT Hollerith edit descriptor spanning continuation lines:")
    assert str(e).endswith("""\
  |      write(6, '(4h|
--------------------^""")
  else: raise Exception_expected
  try:
    fail("hollerith_with_quotes.f")
  except RuntimeError, e:
    assert str(e).startswith(
      "FATAL: Not supported:"
      " FORMAT Hollerith edit descriptor with quotes:")
    assert str(e).endswith("""\
  |      write(6, '(2h'''')')|
--------------------^""")
  else: raise Exception_expected

def exercise_tokens_as_string(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  from tokenization import tokens_as_string
  for file_name in sorted(os.listdir(t_dir)):
    if (not file_name.endswith(".f")): continue
    if (verbose):
      print "exercise_tokens_as_string:", file_name
    units = read.process(file_names=[op.join(t_dir, file_name)])
    for unit in units.all_in_input_order:
      for ei in unit.executable:
        if (ei.key == "write"):
          s = tokens_as_string(tokens=ei.iolist)
          if (verbose):
            print s
      for tokens in unit.format.values():
        s = tokens_as_string(tokens=tokens)
        if (verbose):
          print s
      if (verbose):
        print

def exercise_show():
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  all_units = read.process(file_names=[op.join(t_dir, "subroutine_3.f")])
  from cStringIO import StringIO
  cio = StringIO()
  all_units.show_counts_by_type(out=cio, prefix="$ ")
  assert not show_diff(cio.getvalue(), """\
$ Counts by Fortran unit type:
$   program: 1
$   subroutine: 3
$   function: 0
$   blockdata: 0
""")

def exercise_build_units_by_name():
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  for pair in [
        ("subroutine_3.f", "subroutine_4.f"),
        ("implied_program.f", "implied_program.f")]:
    file_names = [op.join(t_dir, file_name) for file_name in pair]
    all_units = read.process(file_names=file_names)
    from libtbx.utils import Sorry
    try:
      all_units.units_by_name()
    except Sorry, e:
      if (pair[0] == "subroutine_3.f"):
        assert str(e).startswith("Fortran unit name conflict:")
        assert str(e).endswith("  -----------------^")
      else:
        assert str(e).startswith("""\
Fortran unit name conflict:
  1. definition: program_unnamed (implied)
    before """)
        assert str(e).endswith("implied_program.f(2)")
    else: raise Exception_expected

def exercise_eval_const_expression_simple(verbose):
  t_dir = libtbx.env.under_dist(
    module_name="fable", path="test/valid", test=op.isdir)
  file_name = "const_expressions.f"
  all_units = read.process(file_names=[op.join(t_dir, file_name)])
  assert len(all_units.all_in_input_order) == 2
  unit = all_units.all_in_input_order[0]
  val = unit.eval_const_expression_simple(identifier="n5")
  assert val == 296356
  for identifier,expected_vals in [
        ("nums1", [3,2]),
        ("nums2", [1,3]),
        ("nums3", [4])]:
    vals = unit.eval_dimensions_simple(
      dim_tokens=unit.fdecl_by_identifier[identifier].dim_tokens)
    assert vals == expected_vals
  vals = unit.eval_dimensions_simple(
    dim_tokens=unit.fdecl_by_identifier["nums3"].dim_tokens,
    allow_power=False)
  assert vals == [None]
  unit = all_units.all_in_input_order[1]
  vals = unit.eval_dimensions_simple(
    dim_tokens=unit.fdecl_by_identifier["nums"].dim_tokens)
  assert vals == [None, None]

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = (len(args) != 0)
  exercise_strip_spaces_separate_strings()
  exercise_valid(verbose=verbose)
  exercise_lenient(verbose=verbose)
  exercise_syntax_error(verbose=verbose)
  exercise_semantic_error(verbose=verbose)
  exercise_unsupported(verbose=verbose)
  exercise_tokens_as_string(verbose=verbose)
  exercise_show()
  exercise_build_units_by_name()
  exercise_eval_const_expression_simple(verbose=verbose)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
