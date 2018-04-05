from __future__ import absolute_import, division, print_function

from cStringIO import StringIO
import itertools
import os

import fable.read
import pytest

def test_strip_spaces_separate_strings():
  from fable.read import source_line, strip_spaces_separate_strings
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
      assert ssl.code == expected_code
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
    with pytest.raises(fable.read.Error) as e:
      strip_spaces_separate_strings(source_line_cluster=slc(cmbnd))
    assert str(e.value) == """\
Missing terminating %s character:
  at str(1):
  |      %s|
%s^""" % (q, cmbnd, "-"*nd)

def test_valid(testsdir):
  t_dir = os.path.join(testsdir, 'valid')
  #
  read_already = set()
  def get_fprocs(file_name):
    print(file_name)
    read_already.add(file_name)
    return fable.read.process(file_names=[os.path.join(t_dir, file_name)])
  #
  def get_program(file_name):
    fprocs = get_fprocs(file_name)
    assert len(fprocs.program) == 1
    assert len(fprocs.subroutine) == 0
    assert len(fprocs.function) == 0
    assert len(fprocs.blockdata) == 0
    return fprocs.program[0]
  #
  prog = get_program("goto_forms.f")
  keys = [ei.key for ei in prog.executable]
  assert keys == ['goto_computed'] * 4 + ['continue'] * 2

  #
  get_program("string_spanning_continuation_lines.f")
  #
  fprocs = get_fprocs("sf.f")
  assert len(fprocs.program) == 1
  assert len(fprocs.subroutine) == 1
  #
  for file_name in sorted(os.listdir(t_dir)):
    if (file_name.endswith(".f") and not file_name in read_already):
      get_fprocs(file_name)

def test_lenient(testsdir):
  t_dir = os.path.join(testsdir, 'lenient')
  #
  def get(file_name):
    print(file_name)
    return fable.read.process(file_names=[os.path.join(t_dir, file_name)])
  #
  get("str_blank_str.f")
  get("str_cont_line_str.f")

def test_syntax_error(testsdir):
  t_dir = os.path.join(testsdir, 'syntax_error')
  def fail(file_name):
    print(file_name)
    fable.read.process(file_names=[os.path.join(t_dir, file_name)])
  with pytest.raises(fable.read.Error) as e:
    fail("label_cont_char.f")
  assert str(e.value).startswith(
      "A continuation character is illegal on a line with a statement label:")
  assert str(e.value).endswith("""\
  |    3xk=6|
--------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("label_empty.f")
  assert str(e.value).startswith("Labelled statement is empty:")
  assert str(e.value).endswith("""\
  |    1  |
---------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("not_an_identifier.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      external a, x%z|
----------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bare_integer_i.f")
  assert str(e.value).startswith("Missing END for PROGRAM:")
  assert str(e.value).endswith("  |      integer i|")
  with pytest.raises(fable.read.Error) as e:
    fail("closing_parenthesis_without_matching.f")
  assert str(e.value).startswith(
      'Closing ")" without a matching opening parenthesis:')
  assert str(e.value).endswith("""\
  |      a = 0)|
--------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("opening_parenthesis_without_matching.f")
  assert str(e.value).startswith(
      'Missing a closing ")":')
  assert str(e.value).endswith("""\
  |      a = (0|
-------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("dot_e.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |        x       =       .e0|
----------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("x_assign_dot.f")
  assert str(e.value).startswith("Expression unexpectedly ends with a dot:")
  assert str(e.value).endswith("""\
  |      x = .|
-------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("x_assign_dot_2.f")
  assert str(e.value).startswith("Expression unexpectedly ends with a dot:")
  assert str(e.value).endswith("""\
  |        x = .|
---------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("x_assign_1ex.f")
  assert str(e.value).startswith("Invalid floating-point literal:")
  assert str(e.value).endswith("""\
  |      x = 1ex|
---------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("x_assign_1dotd.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      x = 1.d|
---------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_false.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = .fals.|
--------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_true.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = .true|
--------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_not.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = .not b|
--------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_not_2.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = 1 .not. b|
----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_gt.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = 1 .gt 2|
----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_after_dot.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = 1 .2et. 2|
----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("j_assign_i_percent_5.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      j = i % 5|
---------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_and.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = b .ad. c|
-----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_or.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      a = b .or c|
-----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("exclamation_mark_syndrome.f")
  assert str(e.value).startswith("Missing terminating ' character:")
  assert str(e.value).endswith("""\
  |     !ef'|
-----------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("common_with_data_size.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      common /com/ vals(2), nums*4(2)|
-----------------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("dimension_with_data_size.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("  |      dimension strings*4(2)|")
  with pytest.raises(fable.read.Error) as e:
    fail("save_with_dims.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("  |      save nums(2)|")
  with pytest.raises(fable.read.Error) as e:
    fail("sub_no_name.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      subroutine|
------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("sub_percent_3.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      subroutine sub % 3|
------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("sub_open_parenthesis.f")
  assert str(e.value).startswith('Missing a closing ")":')
  assert str(e.value).endswith("""\
  |      subroutine sub(|
-----------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("sub_bad_comma.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      subroutine sub(a,)|
--------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("fun_star.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      function fun(a,*)|
------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("sub_bad_trailing.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      subroutine sub(a,b) x|
-----------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("save_bad_comma.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      save num,|
-----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("save_num_comma_colon.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      save num, :|
-------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("save_num_val_colon.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      save num, val :|
-----------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bare_external.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      external|
----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("save_slash_slash.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      save //|
--------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bare_data.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      data|""")
  with pytest.raises(fable.read.Error) as e:
    fail("data_plus_repetition.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      data nums /+2*3/|
----------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_format_1.f")
  assert str(e.value).startswith('Format string must start with "("')
  assert str(e.value).endswith("""\
  |      write(6, '') num|
------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_format_2.f")
  assert str(e.value).startswith('Format string must end with ")"')
  assert str(e.value).endswith("""\
  |      write(6, '(') num|
------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_format_3.f")
  assert str(e.value).startswith("""\
Missing terminating ' within character format specifier "(')":""")
  assert str(e.value).endswith("""\
  |      write(6, '('')')|
--------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_format_4.f")
  assert str(e.value).startswith("Invalid FORMAT specification:")
  assert str(e.value).endswith("""\
  |      write(6, '(+2x)')|
--------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_format_5.f")
  assert str(e.value).startswith("Invalid FORMAT specification:")
  assert str(e.value).endswith("""\
  |      write(6, '(i2.)')|
----------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_format_6.f")
  assert str(e.value).startswith("Invalid FORMAT specification:")
  assert str(e.value).endswith("""\
  |      write(6, '(tx)')|
---------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_format_7.f")
  assert str(e.value).startswith("Invalid FORMAT specification:")
  assert str(e.value).endswith("""\
  |      write(6, '(tl)')|
---------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("format_without_label.f")
  assert str(e.value).startswith(
      "FORMAT without a statement label in columns 1-5:")
  assert str(e.value).endswith("""\
  |        format(1x)|
-----------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("duplicate_format_label.f")
  assert str(e.value).startswith(
      "Duplicate statement label in columns 1-5:")
  assert str(e.value).endswith("""\
  |   10 format(i2)|
---^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_implied_do_1.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      write(6, *) (i,i=1)|
---------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_implied_do_2.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      write(6, *) (i=1,2)|
-----------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_implied_do_3.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      write(6, *) i,i=1,2|
------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_implied_do_4.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      write(6, *) (i,i=1,j=2)|
-----------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_implied_do_5.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      write(6, *) (i,0=1,2)|
-------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_implied_do_6.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      write(6, *) (i,i+j=1,2)|
---------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bad_data.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      data (nums(i),i) /1,2/|
------------------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("read_star_comma_empty.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      read *,|
---------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("read_star_name.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      read * name|
----------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("read_plus_name.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      read + name|
--------------^""")
  with pytest.raises(fable.read.Error) as e:
    fail("bare_print.f")
  assert str(e.value).startswith("Syntax error:")
  assert str(e.value).endswith("""\
  |      print|
-------------^""")

def test_semantic_error(testsdir):
  t_dir = os.path.join(testsdir, 'semantic_error')
  def fail(file_name):
    print(file_name)
    fable.read.process(file_names=[os.path.join(t_dir, file_name)])
  with pytest.raises(fable.SemanticError) as e:
    fail("missing_include.f")
  assert str(e.value).startswith("Missing include file:")
  assert str(e.value).endswith("""\
  |      include '/bin/sh/should/never/exist'|
-----------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("real_declared_twice.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: val:")
  assert str(e.value).endswith("""\
  |      real val|
--------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("external_array.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: f2:")
  assert str(e.value).endswith("""\
  |      external f2|
------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("dimension_unknown_data_type.f")
  assert str(e.value).startswith("Unknown data type: nums:")
  assert str(e.value).endswith("""\
  |      write(6, *) nums(1)|
---------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("dims_repeated_in_dimension.f")
  assert str(e.value).startswith("Conflicting or repeated dimension: num:")
  assert str(e.value).endswith("""\
  |      dimension num(2)|
-------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("dims_repeated_in_common.f")
  assert str(e.value).startswith("Conflicting or repeated dimension: num:")
  assert str(e.value).endswith("""\
  |      common /com/ num(2)|
----------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("parameter_array.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: nums:")
  assert str(e.value).endswith("""\
  |      integer nums(2)|
-----------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("parameter_in_common.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: num:")
  assert str(e.value).endswith("""\
  |      common /com/ num|
----------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("parameter_save.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: num:")
  assert str(e.value).endswith("""\
  |      save num|
--------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("parameter_in_sub_args.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: num:")
  assert str(e.value).endswith("""\
  |      parameter(num=0)|
-------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("twice_in_sub_args.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: val:")
  assert str(e.value).endswith("""\
  |      subroutine sub(val, val)|
-----------------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("subroutine_name_is_also_arg.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: sub:")
  assert str(e.value).endswith("""\
  |      subroutine sub(sub)|
------------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("function_name_is_also_arg.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: fun:")
  assert str(e.value).endswith("""\
  |      function fun(fun)|
----------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("unknown_intrinsic.f")
  assert str(e.value).startswith("Unknown intrinsic: unk:")
  assert str(e.value).endswith("""\
  |      write(6, *) unk(0)|
---------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("intrinsic_common.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: abs:")
  assert str(e.value).endswith("""\
  |      intrinsic abs|
-------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("external_common.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: nums:")
  assert str(e.value).endswith("""\
  |      external nums|
------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("calling_array.f")
  assert str(e.value).startswith("Conflicting declaration: nums:")
  assert str(e.value).endswith("""\
  |      call nums(3)|
--------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("calling_dimension.f")
  assert str(e.value).startswith("Conflicting declaration: nums:")
  assert str(e.value).endswith("""\
  |      call nums(2)|
--------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("function_data_type_decl_twice.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: fun:")
  assert str(e.value).endswith("""\
  |      integer fun|
-----------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("intrinsic_dimension.f")
  assert str(e.value).startswith("Conflicting or repeated declaration: nint:")
  assert str(e.value).endswith("""\
  |      intrinsic nint|
-------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("sub_fun_2.f")
  assert str(e.value).startswith("Conflicting declaration: fun:")
  assert str(e.value).endswith("""\
  |      y = fun(x)|
-------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("write_without_unit.f")
  assert str(e.value).startswith("Required UNIT information is not defined:")
  assert str(e.value).endswith("""\
  |      write(fmt='(i3)') num|
--------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("write_end.f")
  assert str(e.value).startswith("END is invalid for WRITE statements:")
  assert str(e.value).endswith("""\
  |      write(10, end=20) num|
-----------------------^""")
  with pytest.raises(fable.SemanticError) as e:
    fail("equivalence_external.f")
  assert str(e.value).startswith("Invalid EQUIVALENCE:")
  assert str(e.value).endswith("""\
  |      equivalence (ne, nl)|
----------------------^""")

def test_unsupported(testsdir):
  t_dir = os.path.join(testsdir, 'unsupported')
  def fail(file_name):
    print(file_name)
    fable.read.process(file_names=[os.path.join(t_dir, file_name)])
  with pytest.raises(RuntimeError) as e:
    fail("hollerith_cont_lines.f")
  assert str(e.value).startswith(
      "FATAL: Not supported:"
      " FORMAT Hollerith edit descriptor spanning continuation lines:")
  assert str(e.value).endswith("""\
  |      write(6, '(4h|
--------------------^""")
  with pytest.raises(RuntimeError) as e:
    fail("hollerith_with_quotes.f")
  assert str(e.value).startswith(
      "FATAL: Not supported:"
      " FORMAT Hollerith edit descriptor with quotes:")
  assert str(e.value).endswith("""\
  |      write(6, '(2h'''')')|
--------------------^""")

def test_tokens_as_string(testsdir):
  t_dir = os.path.join(testsdir, 'valid')
  verbose = False
  from fable.tokenization import tokens_as_string
  for file_name in sorted(os.listdir(t_dir)):
    if (not file_name.endswith(".f")): continue
    print(file_name)
    all_fprocs = fable.read.process(file_names=[os.path.join(t_dir, file_name)])
    for fproc in all_fprocs.all_in_input_order:
      for ei in fproc.executable:
        if (ei.key == "write"):
          s = tokens_as_string(tokens=ei.iolist)
          if (verbose):
            print(s)
      for tokens in fproc.format.values():
        s = tokens_as_string(tokens=tokens)
        if (verbose):
          print(s)
      if (verbose):
        print()

def test_show(testsdir):
  t_dir = os.path.join(testsdir, 'valid')
  all_fprocs = fable.read.process(file_names=[os.path.join(t_dir, "subroutine_3.f")])
  cio = StringIO()
  all_fprocs.show_counts_by_type(out=cio, prefix="$ ")
  assert cio.getvalue() == """\
$ Counts by Fortran procedure type:
$   program: 1
$   subroutine: 3
$   function: 0
$   blockdata: 0
"""

def test_build_fprocs_by_name(testsdir):
  t_dir = os.path.join(testsdir, 'valid')
  for pair in [
        ("subroutine_3.f", "subroutine_4.f"),
        ("implied_program.f", "implied_program.f")]:
    file_names = [os.path.join(t_dir, file_name) for file_name in pair]
    all_fprocs = fable.read.process(file_names=file_names)
    with pytest.raises(SystemExit) as e:
      all_fprocs.fprocs_by_name()
    if (pair[0] == "subroutine_3.f"):
      assert str(e.value.code).startswith("Fortran procedure name conflict:")
      assert str(e.value.code).endswith("  -----------------^")
    else:
      assert str(e.value.code).startswith("""\
Fortran procedure name conflict:
  1. definition: program_unnamed (implied)
    before """)
      assert str(e.value.code).endswith("implied_program.f(2)")

def test_eval_const_expression_simple(testsdir):
  t_dir = os.path.join(testsdir, 'valid')
  file_name = "const_expressions.f"
  all_fprocs = fable.read.process(file_names=[os.path.join(t_dir, file_name)])
  assert len(all_fprocs.all_in_input_order) == 2
  fproc = all_fprocs.all_in_input_order[0]
  val = fproc.eval_const_expression_simple(identifier="n5")
  assert val == 296356
  for identifier,expected_vals in [
        ("nums1", [3,2]),
        ("nums2", [1,3]),
        ("nums3", [4])]:
    vals = fproc.eval_dimensions_simple(
      dim_tokens=fproc.fdecl_by_identifier[identifier].dim_tokens)
    assert vals == expected_vals
  vals = fproc.eval_dimensions_simple(
    dim_tokens=fproc.fdecl_by_identifier["nums3"].dim_tokens,
    allow_power=False)
  assert vals == [None]
  fproc = all_fprocs.all_in_input_order[1]
  vals = fproc.eval_dimensions_simple(
    dim_tokens=fproc.fdecl_by_identifier["nums"].dim_tokens)
  assert vals == [None, None]
  #
  file_name = "const_expressions_2.f"
  all_fprocs = fable.read.process(file_names=[os.path.join(t_dir, file_name)])
  assert len(all_fprocs.all_in_input_order) == 1
  fproc = all_fprocs.all_in_input_order[0]
  val = fproc.eval_const_expression_simple(identifier="n5f")
  assert val == pytest.approx(297845.226131, abs=1e-6)
