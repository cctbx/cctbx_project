from __future__ import absolute_import, division, print_function
from six.moves import range
def exercise():
  from libtbx.test_utils import show_diff, Exception_expected
  from six.moves import cPickle as pickle
  #
  from libtbx.str_utils import split_keeping_spaces
  assert split_keeping_spaces(s="") == []
  assert split_keeping_spaces(s=" ") == [" "]
  assert split_keeping_spaces(s="a") == ["a"]
  assert split_keeping_spaces(s="abc") == ["abc"]
  assert split_keeping_spaces(s=" a") == [" ", "a"]
  assert split_keeping_spaces(s="  a") == ["  ", "a"]
  assert split_keeping_spaces(s="  abc") == ["  ", "abc"]
  assert split_keeping_spaces(s="  abc ") == ["  ", "abc", " "]
  assert split_keeping_spaces(s="  abc  ") == ["  ", "abc", "  "]
  assert split_keeping_spaces(s="a ") == ["a", " "]
  assert split_keeping_spaces(s="a  ") == ["a", "  "]
  assert split_keeping_spaces(s="abc  ") == ["abc", "  "]
  assert split_keeping_spaces(s="a b") == ["a", " ", "b"]
  assert split_keeping_spaces(s="a  b") == ["a", "  ", "b"]
  assert split_keeping_spaces(s="  a  b c   d ") == [
    "  ", "a", "  ", "b", " ", "c", "   ", "d", " "]
  #
  from libtbx.str_utils import size_as_string_with_commas
  assert size_as_string_with_commas(0) == "0"
  assert size_as_string_with_commas(1) == "1"
  assert size_as_string_with_commas(-1) == "-1"
  assert size_as_string_with_commas(10) == "10"
  assert size_as_string_with_commas(100) == "100"
  assert size_as_string_with_commas(1000) == "1,000"
  assert size_as_string_with_commas(12345) == "12,345"
  assert size_as_string_with_commas(12345678) == "12,345,678"
  assert size_as_string_with_commas(-12345678) == "-12,345,678"
  #
  from libtbx.str_utils import show_string
  assert show_string("abc") == '"abc"'
  assert show_string("a'c") == '"a\'c"'
  assert show_string('a"c') == "'a\"c'"
  assert show_string('\'"c') == '"\'\\"c"'
  #
  from libtbx.str_utils import prefix_each_line
  assert prefix_each_line(prefix="^", lines_as_one_string="""\
hello
world""") == """\
^hello
^world"""
  #
  from libtbx.str_utils import prefix_each_line_suffix
  assert prefix_each_line_suffix(prefix="^", lines_as_one_string="""\
hello
world""", suffix=" ") == """\
^hello
^world"""
  assert prefix_each_line_suffix(prefix="^", lines_as_one_string="""\
hello
world""", suffix=" ", rstrip=False) == """\
^hello%s
^world """ % " "
  #
  from libtbx.str_utils import show_sorted_by_counts
  from six.moves import cStringIO
  out = cStringIO()
  assert show_sorted_by_counts(
    label_count_pairs=[("b", 3), ("a", 3), ("c", -2)],
    out=out, prefix="%")
  assert not show_diff(out.getvalue(), """\
%"a"  3
%"b"  3
%"c" -2
""")
  out = cStringIO()
  assert show_sorted_by_counts(
    label_count_pairs=[("b", -3), ("a", -3), ("c", 2)], reverse=False,
     out=out, prefix="%", annotations=[None, "", "x"])
  assert not show_diff(out.getvalue(), """\
%"a" -3
%"b" -3
%"c"  2 x
""")
  #
  from libtbx.str_utils import line_breaker
  for string, expected_result in [
    ("", [""]),
    ("this is", ["this is"]),
    ("this is a", ["this is", "a"]),
    ("this is a sentence", ["this is", "a", "sentence"]),
    ("this is a longer sentence", ["this is", "a", "longer", "sentence"]),
    ("this is a very long sentence indeed",
      ["this is", "a very", "long", "sentence", "indeed"])]:
    assert [block for block in line_breaker(string, width=7)]==expected_result
  #
  from libtbx.str_utils import StringIO
  out1 = cStringIO()
  out2 = StringIO()
  out3 = StringIO("Hello world!\n")
  print("Hello world!", file=out1)
  print("Hello world!", file=out2)
  try :
      print("Hello world!", file=out3)
  except AttributeError :
    pass
  else :
    raise Exception_expected
  out4 = pickle.loads(pickle.dumps(out2))
  out5 = pickle.loads(pickle.dumps(out3))
  assert out4.getvalue()==out1.getvalue()==out2.getvalue()==out5.getvalue()
  #
  from libtbx.str_utils import reformat_terminal_text
  txt1 = """
This is some
terminal-formatted
text which needs
to be reset.
"""
  assert (reformat_terminal_text(txt1) ==
          "This is some terminal-formatted text which needs to be reset.")
  txt2 = """
  This is more
  terminal-formatted
  text which needs
  to be reset.
"""
  #
  from libtbx.str_utils import strip_lines, rstrip_lines
  lines = ["  This is more ", "  terminal-formatted ", "  text "]
  assert (strip_lines(txt2) ==
    "\nThis is more\nterminal-formatted\ntext which needs\nto be reset.")
  assert (rstrip_lines(txt2) ==
    "\n  This is more\n  terminal-formatted\n  text which needs\n  to be reset."
  )
  #
  from libtbx.str_utils import expandtabs_track_columns
  def check(s):
    es,js = expandtabs_track_columns(s=s)
    assert len(js) == len(s)
    assert es == s.expandtabs()
    sr = "".join([es[j] for j in js])
    assert sr == s.replace("\t", " ")
  check("")
  check("\t")
  check("\t\t")
  check("\ty")
  check("x\ty")
  check("x\ty\tz")
  check("\txy\t\tz")
  check("abcdefg\txy\t\tz")
  check("ab defgh\txyz\t\tu")
  #
  from libtbx.str_utils import format_value
  assert format_value("%.4f", 1.2345678) == "1.2346"
  assert format_value("%.4f", None) == "  None"
  assert format_value("%.4f", None, replace_none_with="---") == "   ---"
  #
  from libtbx.str_utils import make_header
  out = StringIO()
  make_header("Header 1", out=out)
  assert (out.getvalue() == """
=================================== Header 1 ==================================

""")
  out = StringIO()
  make_header("Header 2", out=out)
  assert (out.getvalue() == """
=================================== Header 2 ==================================

""")
  #
  import sys
  from libtbx.str_utils import string_representation
  iset = list(range(130)) + list(range(250,256))
  for i in iset:
    s = chr(i)
    for j in iset:
      ss = s + chr(j)
      sr = string_representation(
        string=ss, preferred_quote="'", alternative_quote='"')
      if sys.hexversion < 0x03000000:
        assert sr == repr(ss)
      else:
        assert eval(sr) == ss
  from libtbx.str_utils import framed_output
  out = StringIO()
  box = framed_output(out, frame='#')
  print("Hello, world!", file=box)
  box.close()
  assert (out.getvalue() == """
#################
# Hello, world! #
#################
""")
  out = StringIO()
  box = framed_output(out, frame='-', width=80, center=True,
    title="Refinement stats")
  box.write("r_free = 0.1234")
  box.write("  ")
  box.write("r_work = 0.1567")
  box.close()
  assert (out.getvalue() == """
|--------------------------------Refinement stats------------------------------|
|                       r_free = 0.1234  r_work = 0.1567                       |
|------------------------------------------------------------------------------|
""")
  out = StringIO()
  box = framed_output(out, frame='-', width=72, prefix="    ",
    title="Validation summary")
  print("Overall MolProbity score: 2.56", file=box)
  box.add_separator()
  print("""\
Ramachandran favored:  97.5 %
             outliers:  2.5 %
Rotamer outliers:       5.9 %
Clashscore:            10.9""", file=box)
  assert (out.getvalue() == "")
  del box
  assert (out.getvalue() == """
    |-Validation summary---------------------------------------------------|
    | Overall MolProbity score: 2.56                                       |
    |----------------------------------------------------------------------|
    | Ramachandran favored:  97.5 %                                        |
    |              outliers:  2.5 %                                        |
    | Rotamer outliers:       5.9 %                                        |
    | Clashscore:            10.9                                          |
    |----------------------------------------------------------------------|
""")
  from libtbx.str_utils import print_message_in_box
  out = StringIO()
  print_message_in_box(
    message="This is some terminal-formatted text which needs to be reset.",
    out=out,
    width=32,
    center=True,
    prefix="  ",
    frame='*')
  assert (out.getvalue() == """
  ********************************
  *         This is some         *
  *   terminal-formatted text    *
  *   which needs to be reset.   *
  ********************************
""")
  from libtbx.str_utils import make_big_header
  out = StringIO()
  make_big_header("Section title", out=out)
  assert (out.getvalue() == """
################################################################################
#                                Section title                                 #
################################################################################
""")

def exercise_matching_nested_pairs():
  from libtbx.str_utils import find_matching_closing_symbol
  f = find_matching_closing_symbol('(',')')

  s = '... ( ... ( ... ) ...'
  p = s.find('(')
  q = f(s, p)
  assert q == -1
  p1 = s.find('(', p+1)
  q1 = f(s, p1)
  assert q1 == s.rfind(')')

  s = ' ( ( ) ( ( ( ) ) ) ) ) '
  p = f(s, 1)
  assert p == 19

def run(args):
  assert len(args) == 0
  exercise()
  exercise_matching_nested_pairs()
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
