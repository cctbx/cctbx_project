import iotbx.parameters
from cStringIO import StringIO
import sys

class recycle:

  def __init__(self,
        input_string,
        prefix="",
        attributes_level=0,
        print_width=None,
        expected_out=None):
    self.parameters = iotbx.parameters.parse(input_string=input_string)
    self.out = StringIO()
    self.parameters.show(
      out=self.out,
      prefix=prefix,
      attributes_level=attributes_level,
      print_width=print_width)
    self.out = self.out.getvalue()
    if (expected_out is None or self.out != expected_out):
      sys.stdout.write(self.out)
    if (expected_out is not None and self.out != expected_out):
      raise RuntimeError("out != expected_out")
    if (prefix == ""):
      out_parameters = iotbx.parameters.parse(input_string=self.out)
      out_out = StringIO()
      out_parameters.show(
        out=out_out,
        prefix=prefix,
        attributes_level=attributes_level,
        print_width=print_width)
      out_out = out_out.getvalue()
      if (out_out != self.out):
        print "self.out:"
        sys.stdout.write(self.out)
        print "out_out:"
        sys.stdout.write(out_out)
        raise RuntimeError("out_out != self.out")

def exercise_parse_and_show():
  for input_string in ["", "\n", "   \n", "   \t \n \t "]:
    recycle(input_string=input_string, expected_out="")
  recycle(
    input_string="   name\tvalue\n\n",
    expected_out="name value\n")
  r = recycle(
    input_string="   name\tvalue\n\n",
    attributes_level=3,
    expected_out="""\
name value
  .help None
  .caption None
  .short_caption None
  .required None
  .type None
  .input_size None
  .expert_level None
""")
  input_string = """
name value
.help help message with detailed information
.required True
.type path
"""
  recycle(input_string=input_string, attributes_level=3, expected_out="""\
name value
  .help "help message with detailed information"
  .caption None
  .short_caption None
  .required %s
  .type "path"
  .input_size None
  .expert_level None
""" % str(True))
  recycle(input_string=input_string, attributes_level=2, expected_out="""\
name value
  .help "help message with detailed information"
  .required %s
  .type "path"
""" % str(True))
  recycle(input_string=input_string, attributes_level=1, expected_out="""\
name value
  .help "help message with detailed information"
""")
  recycle(input_string=input_string, attributes_level=0, expected_out="""\
name value
""")
  recycle(input_string=input_string, attributes_level=1, print_width=25,
    expected_out="""\
name value
  .help "help message"
        "with detailed"
        "information"
""")
  recycle(
    input_string="name { }", expected_out="""\
name {
}
""")
  recycle(
    input_string="name { a b \n c d }", expected_out="""\
name {
  a b
  c d
}
""")
  recycle(
    input_string="table name { }", expected_out="""\
table name
{
}
""")
  recycle(
    input_string="table name { }", attributes_level=3, expected_out="""\
table name
  .style None
  .help None
  .caption None
  .short_caption None
  .sequential_format None
  .disable_add None
  .disable_delete None
  .expert_level None
{
}
""")
  input_string = """\
table name
  .style column
  .help "help message with detailed information"
  .caption None
  .short_caption None
  .sequential_format None
  .disable_add None
  .disable_delete None
  .expert_level 2
{}
"""
  recycle(input_string=input_string, attributes_level=3,
    expected_out="""\
table name
  .style "column"
  .help "help message with detailed information"
  .caption None
  .short_caption None
  .sequential_format None
  .disable_add None
  .disable_delete None
  .expert_level 2
{
}
""")
  recycle(input_string=input_string, attributes_level=2,
    expected_out="""\
table name
  .style "column"
  .help "help message with detailed information"
  .expert_level 2
{
}
""")
  recycle(input_string=input_string, attributes_level=1,
    expected_out="""\
table name
  .help "help message with detailed information"
{
}
""")
  recycle(input_string=input_string, attributes_level=0,
    expected_out="""\
table name
{
}
""")
  recycle(
    input_string="table name { }", expected_out="""\
table name
{
}
""")
  recycle(
    input_string="table name { {} }", expected_out="""\
table name
{
  {
  }
}
""")
  recycle(
    input_string="table name { row_1 {} }", expected_out="""\
table name
{
  row_1 {
  }
}
""")
  recycle(
    input_string="table name { {} {} }", expected_out="""\
table name
{
  {
  }
  {
  }
}
""")
  recycle(
    input_string="table name { {} row_2 {} }", expected_out="""\
table name
{
  {
  }
  row_2 {
  }
}
""")
  recycle(
    input_string="table name { row_1 {} row_2 {} }", expected_out="""\
table name
{
  row_1 {
  }
  row_2 {
  }
}
""")
  recycle(
    input_string="a b\nc d\n e {} table name { row_1 {} row_2 {} } f g",
    prefix=" prefix ",
    expected_out="""\
 prefix a b
 prefix c d
 prefix
 prefix e {
 prefix }
 prefix
 prefix table name
 prefix {
 prefix   row_1 {
 prefix   }
 prefix   row_2 {
 prefix   }
 prefix }
 prefix
 prefix f g
""")

def test_exception(input_string, exception_string):
  try: iotbx.parameters.parse(input_string=input_string)
  except Exception, e:
    if (exception_string is None or str(e) != exception_string):
      print (e)
    if (exception_string is not None):
      assert str(e) == exception_string
  else: raise RuntimeError("Exception expected.")

def exercise_syntax_errors():
  test_exception("a",
    'Unexpected end of input.')
  test_exception("{}",
    'Syntax error: unexpected "{" (input line 1)')
  test_exception("a {",
    'Syntax error: missing "}".')
  test_exception("table name\n.foo none {}",
    'Unexpected table attribute: .foo (input line 2)')
  test_exception("table name foo {",
    'Syntax error: expected "{", found "foo" (input line 1)')
  test_exception("table name { row_1 bar",
    'Syntax error: expected "{", found "bar" (input line 1)')
  test_exception("a b\n.foo none",
    'Unexpected attribute: .foo (input line 2)')
  test_exception('a b\nc "abc',
    'Syntax error: missing closing quote (input line 2).')

def exercise():
  exercise_parse_and_show()
  exercise_syntax_errors()
  print "OK"

if (__name__ == "__main__"):
  exercise()
