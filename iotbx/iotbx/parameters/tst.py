import iotbx.parameters
from cStringIO import StringIO
import sys, os

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
  recycle(
    input_string='name "value\\\\"',
    expected_out='name "value\\\\"\n')
  recycle(
    input_string='name "value\\\\\\""',
    expected_out='name "value\\\\\\""\n')
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
name
{
}
""")
  recycle(
    input_string="name\n.help message\n.expert_level 3\n{ }",
    attributes_level=3,
    expected_out="""\
name
  .help "message"
  .caption None
  .short_caption None
  .expert_level 3
{
}
""")
  recycle(
    input_string="name{a b\nc d}", expected_out="""\
name
{
  a b
  c d
}
""")
  recycle(
    input_string="table name{}", expected_out="""\
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
    input_string="table name{{}}", expected_out="""\
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
 prefix e
 prefix {
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

def test_exception(input_string, exception_string=None):
  try: iotbx.parameters.parse(input_string=input_string)
  except Exception, e:
    if (exception_string is None or str(e) != exception_string):
      print str(e)
      if (exception_string is not None):
        print exception_string
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
  test_exception('table 1',
    'Syntax error: improper table name "1" (input line 1)')
  test_exception('table name { 1',
    'Syntax error: improper table row name "1" (input line 1)')
  test_exception('1 {',
    'Syntax error: improper scope name "1" (input line 1)')
  test_exception('scope\n.junk',
    'Unexpected scope attribute: .junk (input line 2)')
  test_exception('a. 2',
    'Syntax error: improper definition name "a." (input line 1)')

def check_get(parameters, path, expected_out=None, with_substitution=False):
  out = StringIO()
  if (not with_substitution):
    parameters.get(path=path).show(out=out)
  else:
    parameters.get_with_variable_substitution(path=path).show(out=out)
  out = out.getvalue()
  if (expected_out is None or out != expected_out):
    sys.stdout.write(out)
  if (expected_out is not None and out != expected_out):
    raise RuntimeError("out != expected_out")

def check_get_sub(parameters, path, expected_out=None):
  check_get(parameters, path, expected_out, with_substitution=True)

def exercise_get():
  parameters = iotbx.parameters.parse(input_string="""\
a b
c d
e {
  a 1
  b x
}
table name {
  row_1 {
    a 1
    b x
  }
  row_2 {
    a 2
    b y
  }
  {
    a 3
    b z
  }
}
e g""")
  check_get(parameters, path="a", expected_out="""\
a b
""")
  check_get(parameters, path="e", expected_out="""\
e
{
  a 1
  b x
}

e g
""")
  check_get(parameters, path="e.a", expected_out="""\
a 1
""")
  check_get(parameters, path="e.b", expected_out="""\
b x
""")
  check_get(parameters, path="name", expected_out="""\
table name
{
  row_1 {
    a 1
    b x
  }
  row_2 {
    a 2
    b y
  }
  {
    a 3
    b z
  }
}
""")
  check_get(parameters, path="name.1", expected_out="""\
a 1
b x
""")
  check_get(parameters, path="name.row_1", expected_out="""\
a 1
b x
""")
  check_get(parameters, path="name.2", expected_out="""\
a 2
b y
""")
  check_get(parameters, path="name.row_2", expected_out="""\
a 2
b y
""")
  check_get(parameters, path="name.3", expected_out="""\
a 3
b z
""")
  check_get(parameters, path="name.row_3", expected_out="")

def exercise_nested():
  parameters = recycle(
    input_string="""\
a0 {
  d1 a b c
  a1 {
    table t0 {
      row0
      {
        c yes
        table t1 {
          {
            x 0
            y 1.
          }
        }
      }
    }
  }
  d2 e f 0g
}
""",
    expected_out="""\
a0
{
  d1 a b c

  a1
  {
    table t0
    {
      row0 {
        c yes

        table t1
        {
          {
            x 0
            y 1.
          }
        }
      }
    }
  }

  d2 e f 0g
}
""").parameters
  check_get(parameters, path="a0.d1", expected_out="d1 a b c\n")
  check_get(parameters, path="a0.a1.t0.row0.c", expected_out="c yes\n")
  check_get(parameters, path="a0.a1.t0.row0.t1.1.x", expected_out="x 0\n")
  check_get(parameters, path="a0.a1.t0.row0.t1.1.y", expected_out="y 1.\n")
  parameters.automatic_type_assignment(assignment_if_unknown="UNKNOWN")
  out = StringIO()
  parameters.show(out=out, attributes_level=2)
  assert out.getvalue() == """\
a0
{
  d1 a b c
    .type "str"

  a1
  {
    table t0
    {
      row0 {
        c yes
          .type "bool"

        table t1
        {
          {
            x 0
              .type "int"
            y 1.
              .type "float"
          }
        }
      }
    }
  }

  d2 e f 0g
    .type "UNKNOWN"
}
"""

def exercise_get_with_variable_substitution():
  parameters = iotbx.parameters.parse(input_string="""\
a b
c   d   e   2
""")
  check_get_sub(parameters, path="a", expected_out="a b\n")
  check_get_sub(parameters, path="c", expected_out="c d e 2\n")
  parameters = iotbx.parameters.parse(input_string="""\
a $a
b $c
c $b
d $e
e $f
f $d
g $d
h $_X_Y_Z_
i 0
i 1
j $i
k $l
l x
l y
""")
  assert len(parameters.get(path="a").objects) == 1
  for path in "abcdefg":
    if (path == "g"):
      diag = "d"
    else:
      diag = path
  try: parameters.get_with_variable_substitution(path=path)
  except RuntimeError, e:
    assert str(e) == "Dependency cycle in variable substitution: $%s" % diag
  else: raise RuntimeError("Exception expected.")
  try: parameters.get_with_variable_substitution(path="h")
  except RuntimeError, e:
    assert str(e) == "Undefined variable: $_X_Y_Z_ (input line 8)"
  else: raise RuntimeError("Exception expected.")
  os.environ["_X_Y_Z_"] = "xyz"
  check_get_sub(parameters, path="h", expected_out='h "xyz"\n')
  check_get_sub(parameters, path="j", expected_out='j 1\n')
  check_get_sub(parameters, path="k", expected_out='k y\n')
  parameters = iotbx.parameters.parse(input_string="""\
a x
b $a
c $a $b.d
d $a \$b
answer yes no
e "$a"
f $(a)
g abc$(answer)bc
h abc$(answer)bc 12$answer$(a)56
i $
j $(abc
k $(1bc)
l $()
m $@
n '$a'
""")
  check_get_sub(parameters, path="a", expected_out="a x\n")
  check_get_sub(parameters, path="b", expected_out="b x\n")
  check_get_sub(parameters, path="c", expected_out='c x "x.d"\n')
  check_get_sub(parameters, path="d", expected_out="d x \\$b\n")
  check_get_sub(parameters, path="answer", expected_out="answer yes no\n")
  check_get_sub(parameters, path="e", expected_out='e "x"\n')
  check_get_sub(parameters, path="f", expected_out='f "x"\n')
  check_get_sub(parameters, path="g", expected_out='g "abcyes nobc"\n')
  check_get_sub(parameters, path="h",
    expected_out='h "abcyes nobc" "12yes nox56"\n')
  try: parameters.get_with_variable_substitution(path="i")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: $ must be followed by an identifier:' \
                   + ' "$" (input line 10)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get_with_variable_substitution(path="j")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: missing ")": "$(abc" (input line 11)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get_with_variable_substitution(path="k")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: improper variable name "$(1bc)"' \
                   + ' (input line 12)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get_with_variable_substitution(path="l")
  except RuntimeError, e:
    assert str(e)=='Syntax error: improper variable name "$()" (input line 13)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get_with_variable_substitution(path="m")
  except RuntimeError, e:
    assert str(e)=='Syntax error: improper variable name "$@" (input line 14)'
  else: raise RuntimeError("Exception expected.")
  check_get_sub(parameters, path="n", expected_out="n '$a'\n")

def exercise_include():
  print >> open("tmp1.params", "w"), """\
include tmp1.params
a x
"""
  print >> open("tmp2.params", "w"), """\
b y
"""
  print >> open("tmp3.params", "w"), """\
c z
include tmp2.params
d $z
"""
  parameters = iotbx.parameters.parse(
    input_string="""\
tmp2 tmp2.params
include tmp1.params $tmp2
r 0
include tmp3.params
s 1
""",
    process_includes=True)
  out = StringIO()
  parameters.show(out=out)
  assert out.getvalue() == """\
tmp2 tmp2.params
a x
b y
r 0
c z
d $z
s 1
"""
  try: parameters.get_with_variable_substitution(path="d")
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $z (file "tmp3.params", line 3)'
  else: raise RuntimeError("Exception expected.")
  try: os.makedirs("tmp")
  except OSError: pass
  print >> open("tmp1.params", "w"), """\
a 0
include tmp/tmp1.params
x 1
"""
  print >> open("tmp/tmp1.params", "w"), """\
b 1
include tmp2.params
y 2
"""
  print >> open("tmp/tmp2.params", "w"), """\
c 2
include tmp1.params
z 3
"""
  parameters = iotbx.parameters.parse(file_name="tmp1.params", process_includes=True)
  out = StringIO()
  parameters.show(out=out)
  assert out.getvalue() == """\
a 0
b 1
c 2
z 3
y 2
x 1
"""

def exercise_extract():
  parameters = iotbx.parameters.parse(input_string="""\
group {
  a yes
    .type bool
  a yes
  b 13
    .type int
  c 1.3
    .type float
  d abc def ghi
    .type str
  e a *b c
    .type choice
  e a *b *c
    .type choice
  e a b c
    .type choice
  f a *b c
    .type multi_choice
  f a *b *c
    .type multi_choice
  f a b c
    .type multi_choice
  g "/var/tmp/junk"
    .type path
  h "var.tmp.junk"
    .type key
  u 10,12 13 80,90 100
    .type unit_cell
  s 19
    .type space_group
}
""")
  assert parameters.get(path="group.a").objects[0].extract() is True
  try: parameters.get(path="group.a").objects[1].extract()
  except RuntimeError, e:
    assert str(e) == 'Undefined parameter definition type for' \
                   + ' values of "a" (input line 4).'
  else: raise RuntimeError("Exception expected.")
  assert parameters.get(path="group.b").objects[0].extract() == 13
  assert parameters.get(path="group.c").objects[0].extract() == 1.3
  assert parameters.get(path="group.d").objects[0].extract() == "abc def ghi"
  assert parameters.get(path="group.e").objects[0].extract() == "b"
  try: parameters.get(path="group.e").objects[1].extract()
  except RuntimeError, e:
    assert str(e) \
      == 'Multiple choices where only one is possible (input line 13).'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="group.e").objects[2].extract()
  except RuntimeError, e:
    assert str(e) == 'Unspecified choice (input line 15).'
  else: raise RuntimeError("Exception expected.")
  assert parameters.get(path="group.f").objects[0].extract() == ["b"]
  assert parameters.get(path="group.f").objects[1].extract() == ["b", "c"]
  assert parameters.get(path="group.f").objects[2].extract() == []
  assert parameters.get(path="group.g").objects[0].extract() == "/var/tmp/junk"
  assert parameters.get(path="group.h").objects[0].extract() == "var.tmp.junk"
  assert str(parameters.get(path="group.u").objects[0].extract()) \
      == "(10, 12, 13, 80, 90, 100)"
  assert str(parameters.get(path="group.s").objects[0].extract()) \
      == "P 21 21 21"
  definition = parameters.get(path="group.a").objects[0]
  definition.type = "junk"
  try: definition.extract()
  except RuntimeError, e:
    str(e) == 'No converter for parameter definition type "junk"' \
            + ' required for converting values of "a" (input line 3).'
  else: raise RuntimeError("Exception expected.")
  parameters = iotbx.parameters.parse(input_string="""\
group {
  a yes
    .type bool
  b 13
    .type int
  s P 21 21 21
    .type space_group
}
""")
  group = parameters.get(path="group").objects[0].extract()
  assert group.a is True
  assert group.b == 13
  assert group.s.type().number() == 19

def exercise():
  exercise_parse_and_show()
  exercise_syntax_errors()
  exercise_get()
  exercise_nested()
  exercise_get_with_variable_substitution()
  exercise_include()
  exercise_extract()
  print "OK"

if (__name__ == "__main__"):
  exercise()
