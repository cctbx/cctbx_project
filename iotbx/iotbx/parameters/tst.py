import iotbx.parameters
from cStringIO import StringIO
import copy
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
    input_string="   name\t=value\n\n",
    expected_out="name = value\n")
  recycle(
    input_string='name="value\\\\"',
    expected_out='name = "value\\\\"\n')
  recycle(
    input_string='name="value\\\\\\""',
    expected_out='name = "value\\\\\\""\n')
  r = recycle(
    input_string="   name\t=\tvalue\n\n",
    attributes_level=3,
    expected_out="""\
name = value
  .help = None
  .caption = None
  .short_caption = None
  .required = None
  .type = None
  .input_size = None
  .expert_level = None
""")
  input_string = """
name=value
.help=help message with detailed information
.required=True
.type=path
"""
  recycle(input_string=input_string, attributes_level=3, expected_out="""\
name = value
  .help = "help message with detailed information"
  .caption = None
  .short_caption = None
  .required = %s
  .type = "path"
  .input_size = None
  .expert_level = None
""" % str(True))
  recycle(input_string=input_string, attributes_level=2, expected_out="""\
name = value
  .help = "help message with detailed information"
  .required = %s
  .type = "path"
""" % str(True))
  recycle(input_string=input_string, attributes_level=1, expected_out="""\
name = value
  .help = "help message with detailed information"
""")
  recycle(input_string=input_string, attributes_level=0, expected_out="""\
name = value
""")
  recycle(input_string=input_string, attributes_level=1, print_width=25,
    expected_out="""\
name = value
  .help = "help message"
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
    input_string="name\n.help=message\n.expert_level=3\n{ }",
    attributes_level=3,
    expected_out="""\
name
  .style = None
  .help = "message"
  .caption = None
  .short_caption = None
  .required = None
  .sequential_format = None
  .disable_add = None
  .disable_delete = None
  .expert_level = 3
{
}
""")
  recycle(
    input_string="name{a=b\nc=d\n}", expected_out="""\
name
{
  a = b
  c = d
}
""")
  recycle(
    input_string=\
      "a=b\nc=d\n e {} name { var1=None\nvar2=None\n} f=g",
    prefix=" prefix ",
    expected_out="""\
 prefix a = b
 prefix c = d
 prefix
 prefix e
 prefix {
 prefix }
 prefix
 prefix name
 prefix {
 prefix   var1 = None
 prefix   var2 = None
 prefix }
 prefix
 prefix f = g
""")
  recycle(
    input_string="""\
#name # 1 2 3
  #.help="a"
        "b"
        "c"
      # 1 2 3
  .expert_level=1 # x
        # 1 2 3
{
  #a=b
    .type="int" # y
  c=d # e
    #.help="d"
          "e"
          "f"
    .expert_level=2
}
d=a b c # 1 2 3 \\
4 5 6
""",
    attributes_level=2,
    expected_out="""\
#name
  .expert_level = 1
{
  #a = b
    .type = "int"
  c = d
    .expert_level = 2
}

d = a b c
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
  test_exception("'a'",
    """Unquoted word expected, found 'a' (input line 1)""")
  test_exception("a",
    'Unexpected end of input.')
  test_exception("a=\nb",
    'Missing value for a (input line 1)')
  test_exception("a b",
    'Syntax error: expected "=", found "b" (input line 1)')
  test_exception("{}",
    'Syntax error: unexpected "{" (input line 1)')
  test_exception("a {",
    'Syntax error: no matching "}" for "{" at input line 1')
  test_exception("a=b\n.foo none",
    'Unexpected definition attribute: .foo (input line 2)')
  test_exception('a=b\nc "abc',
    'Syntax error: missing closing quote (input line 2)')
  test_exception('1 {',
    'Syntax error: improper scope name "1" (input line 1)')
  test_exception('scope\n.foo',
    'Unexpected scope attribute: .foo (input line 2)')
  test_exception('a. 2',
    'Syntax error: improper definition name "a." (input line 1)')

def check_get(parameters, path, expected_out=None, with_substitution=False):
  out = StringIO()
  parameters.get(
    path=path,
    with_substitution=with_substitution).show(out=out)
  out = out.getvalue()
  if (expected_out is None or out != expected_out):
    sys.stdout.write(out)
  if (expected_out is not None and out != expected_out):
    raise RuntimeError("out != expected_out")

def check_get_sub(parameters, path, expected_out=None):
  check_get(parameters, path, expected_out, with_substitution=True)

def exercise_get_without_substitution():
  parameters = iotbx.parameters.parse(input_string="""\
a=b
c=d
e {
  a=1
  b=x
}
e=g""")
  check_get(parameters, path="a", expected_out="""\
a = b
""")
  check_get(parameters, path="e", expected_out="""\
e
{
  a = 1
  b = x
}

e = g
""")
  check_get(parameters, path="e.a", expected_out="""\
a = 1
""")
  check_get(parameters, path="e.b", expected_out="""\
b = x
""")
  check_get(parameters, path="e.c", expected_out="")

def exercise_nested():
  parameters = recycle(
    input_string="""\
d0=0
a0 {
  d1=a b c
  a1 {
    t0 {
      c=yes
      t1 {
        x=0
        y=1.
      }
    }
  }
  d2=e f 0g
}
""",
    expected_out="""\
d0 = 0

a0
{
  d1 = a b c

  a1
  {
    t0
    {
      c = yes

      t1
      {
        x = 0
        y = 1.
      }
    }
  }

  d2 = e f 0g
}
""").parameters
  check_get(parameters, path="a0.d1", expected_out="d1 = a b c\n")
  check_get(parameters, path="a0.a1.t0.c", expected_out="c = yes\n")
  check_get(parameters, path="a0.a1.t0.t1.x", expected_out="x = 0\n")
  check_get(parameters, path="a0.a1.t0.t1.y", expected_out="y = 1.\n")
  assert [item.path for item in parameters.all_definitions()] == [
    "d0", "a0.d1", "a0.a1.t0.c", "a0.a1.t0.t1.x", "a0.a1.t0.t1.y", "a0.d2"]
  parameters.automatic_type_assignment(assignment_if_unknown="UNKNOWN")
  out = StringIO()
  parameters.show(out=out, attributes_level=2)
  assert out.getvalue() == """\
d0 = 0
  .type = "int"

a0
{
  d1 = a b c
    .type = "str"

  a1
  {
    t0
    {
      c = yes
        .type = "bool"

      t1
      {
        x = 0
          .type = "int"
        y = 1.
          .type = "float"
      }
    }
  }

  d2 = e f 0g
    .type = "UNKNOWN"
}
"""

def exercise_get_with_substitution():
  parameters = iotbx.parameters.parse(input_string="""\
a=b
c = d   e   2
""")
  check_get_sub(parameters, path="a", expected_out="a = b\n")
  check_get_sub(parameters, path="c", expected_out="c = d e 2\n")
  parameters = iotbx.parameters.parse(input_string="""\
a=$a
b=$c
c=$b
d=$e
e=$f
f=$d
g=$d
h=$_X_Y_Z_
i=0
i=1
j=$i
k=$l
l=x
l=y
""")
  assert len(parameters.get(path="a", with_substitution=False).objects) == 1
  for path in "abcdefg":
    if (path == "g"):
      diag = "d"
    else:
      diag = path
  try: parameters.get(path=path)
  except RuntimeError, e:
    assert str(e) == "Dependency cycle in variable substitution: $%s" % diag
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="h")
  except RuntimeError, e:
    assert str(e) == "Undefined variable: $_X_Y_Z_ (input line 8)"
  else: raise RuntimeError("Exception expected.")
  os.environ["_X_Y_Z_"] = "xyz"
  check_get_sub(parameters, path="h", expected_out='h = "xyz"\n')
  check_get_sub(parameters, path="j", expected_out='j = 1\n')
  check_get_sub(parameters, path="k", expected_out='k = y\n')
  parameters = iotbx.parameters.parse(input_string="""\
a=x
b=$a
c=$a $b.d
d=$a \$b
answer=yes no
e="$a"
f=${a}
g=abc${answer}bc
h=abc${answer}bc 12$answer${a}56
i=$
j=${abc
k=${1bc}
l=${}
m=$@
n='$a'
""")
  check_get_sub(parameters, path="a", expected_out="a = x\n")
  check_get_sub(parameters, path="b", expected_out="b = x\n")
  check_get_sub(parameters, path="c", expected_out='c = x "x.d"\n')
  check_get_sub(parameters, path="d", expected_out="d = x \\$b\n")
  check_get_sub(parameters, path="answer", expected_out="answer = yes no\n")
  check_get_sub(parameters, path="e", expected_out='e = "x"\n')
  check_get_sub(parameters, path="f", expected_out='f = x\n')
  check_get_sub(parameters, path="g", expected_out='g = "abcyes nobc"\n')
  check_get_sub(parameters, path="h",
    expected_out='h = "abcyes nobc" "12yes nox56"\n')
  try: parameters.get(path="i")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: $ must be followed by an identifier:' \
                   + ' "$" (input line 10)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="j")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: missing "}": "${abc" (input line 11)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="k")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: improper variable name "${1bc}"' \
                   + ' (input line 12)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="l")
  except RuntimeError, e:
    assert str(e)=='Syntax error: improper variable name "${}" (input line 13)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="m")
  except RuntimeError, e:
    assert str(e)=='Syntax error: improper variable name "$@" (input line 14)'
  else: raise RuntimeError("Exception expected.")
  check_get_sub(parameters, path="n", expected_out="n = '$a'\n")

def exercise_include():
  print >> open("tmp1.params", "w"), """\
include tmp1.params
#include none
a=x
"""
  print >> open("tmp2.params", "w"), """\
b=y
"""
  print >> open("tmp3.params", "w"), """\
c=z
include tmp2.params
d=$z
"""
  parameters = iotbx.parameters.parse(
    input_string="""\
tmp2=tmp2.params
include tmp1.params $tmp2
r=0
include tmp3.params
s=1
""",
    process_includes=True)
  out = StringIO()
  parameters.show(out=out)
  assert out.getvalue() == """\
tmp2 = tmp2.params
#include none
a = x
b = y
r = 0
c = z
d = $z
s = 1
"""
  try: parameters.get(path="d")
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $z (file "tmp3.params", line 3)'
  else: raise RuntimeError("Exception expected.")
  try: os.makedirs("tmp")
  except OSError: pass
  print >> open("tmp1.params", "w"), """\
a=0
include tmp/tmp1.params
x=1
"""
  print >> open("tmp/tmp1.params", "w"), """\
b=1
include tmp2.params
y=2
"""
  print >> open("tmp/tmp2.params", "w"), """\
c=2
include tmp1.params
z=3
"""
  parameters = iotbx.parameters.parse(
    file_name="tmp1.params",
    process_includes=True)
  out = StringIO()
  parameters.show(out=out)
  assert out.getvalue() == """\
a = 0
b = 1
c = 2
z = 3
y = 2
x = 1
"""

def exercise_extract():
  parameters = iotbx.parameters.parse(input_string="""\
group {
  a=yes
    .type=bool
  a=yes "n o"
  b=13
    .type=int
  c=1.3
    .type=float
  d=abc def ghi
    .type=str
  e=a *b c
    .type=choice
  e=a *b *c
    .type=choice
  e=a b c
    .type=choice
    .required=yes
  f=a *b c
    .type=multi_choice
  f=a *b *c
    .type=multi_choice
  f=a b c
    .type=multi_choice
  g="/var/tmp/foo"
    .type=path
  h="var.tmp.foo"
    .type=key
  i=ceil(4/3)
    .type=int
  i=1/2
    .type=int
  j=1/(5-3)
    .type=float
  j="'a'"
    .type=float
  j=a
    .type=float
  u=10,12 13 80,90 100
    .type=unit_cell
  s=19
    .type=space_group
}
""")
  assert parameters.get(path="group.a",
    with_substitution=False).objects[0].extract() is True
  assert parameters.get(path="group.a",
    with_substitution=False).objects[1].extract() == ["yes", "n o"]
  assert parameters.get(path="group.b",
    with_substitution=False).objects[0].extract() == 13
  assert parameters.get(path="group.c",
    with_substitution=False).objects[0].extract() == 1.3
  assert parameters.get(path="group.d",
    with_substitution=False).objects[0].extract() == "abc def ghi"
  assert parameters.get(path="group.e",
    with_substitution=False).objects[0].extract() == "b"
  try: parameters.get(path="group.e",
    with_substitution=False).objects[1].extract()
  except RuntimeError, e:
    assert str(e) \
      == 'Multiple choices where only one is possible (input line 13)'
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="group.e",
    with_substitution=False).objects[2].extract()
  except RuntimeError, e:
    assert str(e) == 'Unspecified choice (input line 15)'
  else: raise RuntimeError("Exception expected.")
  assert parameters.get(path="group.f",
    with_substitution=False).objects[0].extract() == ["b"]
  assert parameters.get(path="group.f",
    with_substitution=False).objects[1].extract() == ["b", "c"]
  assert parameters.get(path="group.f",
    with_substitution=False).objects[2].extract() == []
  assert parameters.get(path="group.g",
    with_substitution=False).objects[0].extract() == "/var/tmp/foo"
  assert parameters.get(path="group.h",
    with_substitution=False).objects[0].extract() == "var.tmp.foo"
  assert parameters.get(path="group.i",
    with_substitution=False).objects[0].extract() == 2
  try: parameters.get(path="group.i",
    with_substitution=False).objects[1].extract()
  except RuntimeError, e:
    assert str(e) == 'Integer expression expected, "1/2" found (input line 30)'
  else: raise RuntimeError("Exception expected.")
  assert parameters.get(path="group.j",
    with_substitution=False).objects[0].extract() == 0.5
  try: parameters.get(path="group.j",
    with_substitution=False).objects[1].extract()
  except RuntimeError, e:
    assert str(e) == \
      """Floating-point expression expected, "'a'" found (input line 34)"""
  else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="group.j",
    with_substitution=False).objects[2].extract()
  except RuntimeError, e:
    assert str(e) == 'Error interpreting "a" as a numeric expression:' \
                   + " name 'a' is not defined (input line 36)"
  else: raise RuntimeError("Exception expected.")
  assert str(parameters.get(path="group.u",
    with_substitution=False).objects[0].extract()) \
      == "(10, 12, 13, 80, 90, 100)"
  assert str(parameters.get(path="group.s",
    with_substitution=False).objects[0].extract()) \
      == "P 21 21 21"
  definition = parameters.get(path="group.a",
    with_substitution=False).objects[0]
  definition.type = "foo"
  try: definition.extract()
  except RuntimeError, e:
    assert str(e) == 'No converter for parameter definition type "foo"' \
      + ' required for converting words assigned to "a" (input line 2)'
  else: raise RuntimeError("Exception expected.")
  parameters = iotbx.parameters.parse(input_string="""\
group {
  a=yes
    .type=bool
  b=13
    .type=int
  s=P 21 21 21
    .type=space_group
}
""")
  group = parameters.get(path="group",
    with_substitution=False).objects[0].extract()
  assert group.a is True
  assert group.b == 13
  assert group.s.type().number() == 19

def exercise_format():
  parameters = iotbx.parameters.parse(input_string="""\
group {
  n=ab "c d" 'ef '
    .type=None
  a=yes
    .type=bool
  b=13
    .type=int
  c=1.3
    .type=float
  d=abc def ghi
    .type=str
  e=a *b c
    .type=choice
  f=a *b *c
    .type=multi_choice
  g="/var/tmp/foo"
    .type=path
  h="var.tmp.foo"
    .type=key
  u=10,12 13 80,90 100
    .type=unit_cell
  s=19
    .type=space_group
}
""")
  out = StringIO()
  parameters.merge_and_extract(source=parameters).format().show(out=out)
  assert out.getvalue() == """\
group
{
  n = "ab" "c d" "ef "
  a = True
  b = 13
  c = 1.3
  d = "abc def ghi"
  e = a *b c
  f = a *b *c
  g = "/var/tmp/foo"
  h = "var.tmp.foo"
  u = 10 12 13 80 90 100
  s = "P 21 21 21"
}
"""
  definition = parameters.get_last(path="group.a")
  definition.type = "foo"
  try: definition.format(python_object=0)
  except RuntimeError, e:
    assert str(e) == 'No converter for parameter definition type "foo"' \
      + ' required for converting values for "a" (input line 4)'
  else: raise RuntimeError("Exception expected.")

def exercise_deepcopy():
  parameters = iotbx.parameters.parse(input_string="""\
  a=1
  b {
    a=1
  }
""")
  copy.deepcopy(parameters)

def exercise_merge():
  master = iotbx.parameters.parse(input_string="""\
a=None
  .expert_level=1
b=None
  .expert_level=2
c
{
  a=None
    .expert_level=3
  b=None
    .expert_level=4
  c
    .expert_level=5
  {
    x=1
      .expert_level=6
    y=2
      .type="int"
      .expert_level=7
    t
      .expert_level=8
    {
      r=3
        .help="help"
    }
  }
}
t
  .expert_level=9
{
  a=4
    .expert_level=10
  b=5
    .expert_level=11
}

u=a *b c
  .type=choice
""")
  source = iotbx.parameters.parse(input_string="""\
a=7
c
{
  a=${v.x}
  b=1
  c
  {
    y=${v.y}
    t
    {
      r=x
    }
  }
}
v
{
  x=y
  y=3
}
t
{
  a=1
    .expert_level=-1
  b=2
}
a=9
t
{
  b=3
    .expert_level=-2
}
u=e b f *c a g
""")
  merged = master.merge(source=source)
  out = StringIO()
  merged.show(out=out, attributes_level=2)
  assert out.getvalue() == """\
a = 7
  .expert_level = 1

c
{
  a = ${v.x}
    .expert_level = 3
  b = 1
    .expert_level = 4

  c
    .expert_level = 5
  {
    y = ${v.y}
      .type = "int"
      .expert_level = 7

    t
      .expert_level = 8
    {
      r = x
        .help = "help"
    }

    x = 1
      .expert_level = 6
  }
}

v
{
  x = y
  y = 3
}

t
  .expert_level = 9
{
  a = 1
    .expert_level = 10
  b = 2
    .expert_level = 11
}

a = 9
  .expert_level = 1

t
  .expert_level = 9
{
  b = 3
    .expert_level = 11
  a = 4
    .expert_level = 10
}

u = a b *c
  .type = "choice"
b = None
  .expert_level = 2
"""
  assert merged.get(path="c.c.y").objects[0].extract() == 3
  master = iotbx.parameters.parse(input_string="""\
c=a *b c
  .type=choice
""")
  source = iotbx.parameters.parse(input_string="""\
c=a *d c
""")
  try: merged = master.merge(sources=[source])
  except RuntimeError, e:
    assert str(e) == "Not a possible choice: *d (input line 1)"
  else: raise RuntimeError("Exception expected.")

def exercise():
  exercise_parse_and_show()
  exercise_syntax_errors()
  exercise_get_without_substitution()
  exercise_nested()
  exercise_get_with_substitution()
  exercise_include()
  exercise_extract()
  exercise_format()
  exercise_deepcopy()
  exercise_merge()
  print "OK"

if (__name__ == "__main__"):
  exercise()
