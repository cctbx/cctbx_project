from libtbx import phil
import libtbx.phil.command_line
from libtbx.utils import UserError
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
    self.parameters = phil.parse(input_string=input_string)
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
      out_parameters = phil.parse(input_string=self.out)
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
  .optional = None
  .type = None
  .multiple = None
  .input_size = None
  .expert_level = None
""")
  input_string = """
name=value
.help=help message with detailed information
.optional=True
.type=path
"""
  recycle(input_string=input_string, attributes_level=3, expected_out="""\
name = value
  .help = "help message with detailed information"
  .caption = None
  .short_caption = None
  .optional = True
  .type = path
  .multiple = None
  .input_size = None
  .expert_level = None
""")
  recycle(input_string=input_string, attributes_level=2, expected_out="""\
name = value
  .help = "help message with detailed information"
  .optional = True
  .type = path
""")
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
name {
}
""")
  recycle(
    input_string="name\n.help=message\n.expert_level=3\n{ }",
    attributes_level=3,
    expected_out="""\
name
  .style = None
  .help = message
  .caption = None
  .short_caption = None
  .optional = None
  .multiple = None
  .sequential_format = None
  .disable_add = None
  .disable_delete = None
  .expert_level = 3
{
}
""")
  recycle(
    input_string="name{a=b\nc=d\n}", expected_out="""\
name {
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
 prefix e {
 prefix }
 prefix name {
 prefix   var1 = None
 prefix   var2 = None
 prefix }
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
    .type = int
  c = d
    .expert_level = 2
}
d = a b c
""")
  #
  params = phil.parse(input_string="""\
x=1
y=2
  .expert_level=1
z=3
  .expert_level=2
s {
  a=x
}
t
  .expert_level=1
{
  a=y
}
u
  .expert_level=2
{
  a=z
}
""")
  assert params.as_str(expert_level=0) == """\
x = 1
s {
  a = x
}
"""
  assert params.as_str(expert_level=1) == """\
x = 1
y = 2
s {
  a = x
}
t {
  a = y
}
"""
  for expert_level in [-1,2,3,None]:
    assert params.as_str(expert_level=expert_level) == """\
x = 1
y = 2
z = 3
s {
  a = x
}
t {
  a = y
}
u {
  a = z
}
"""

def test_exception(input_string, exception_string=None):
  try: phil.parse(input_string=input_string)
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
  test_exception('a.include=None',
    'Reserved identifier: "include" (input line 1)')
  test_exception('include {}',
    'Reserved identifier: "include" (input line 1)')
  test_exception('a.include.b.c {}',
    'Reserved identifier: "include" (input line 1)')
  test_exception('a=None\n.type=foo',
    'Unexpected definition type: "foo" (input line 2)')

def exercise_deepcopy():
  parameters = phil.parse(input_string="""\
  a=1
  b {
    a=1
  }
""")
  copy.deepcopy(parameters)

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
  parameters = phil.parse(input_string="""\
a=b
c=d
e {
  a=1
  b=x
}
e=g""")
  check_get(parameters, path="", expected_out="""\
a = b
c = d
e {
  a = 1
  b = x
}
e = g
""")
  check_get(parameters, path="a", expected_out="""\
a = b
""")
  check_get(parameters, path="e", expected_out="""\
e {
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
  include file
  a1 {
    t0 {
      c=yes
      t1 {
        x=0
        y=1.
          .multiple=True
      }
    }
  }
  d2=e f 0g
  #d3=x
}
""",
    expected_out="""\
d0 = 0
a0 {
  d1 = a b c
  include file
  a1 {
    t0 {
      c = yes
      t1 {
        x = 0
        y = 1.
      }
    }
  }
  d2 = e f 0g
  #d3 = x
}
""").parameters
  check_get(parameters, path="a0.d1", expected_out="d1 = a b c\n")
  check_get(parameters, path="a0.a1.t0.c", expected_out="c = yes\n")
  check_get(parameters, path="a0.a1.t0.t1.x", expected_out="x = 0\n")
  check_get(parameters, path="a0.a1.t0.t1.y", expected_out="y = 1.\n")
  assert [item.path for item in parameters.all_definitions()] == [
    "d0", "a0.d1", "a0.a1.t0.c", "a0.a1.t0.t1.x", "a0.a1.t0.t1.y", "a0.d2"]
  assert [item.path for item in parameters.all_definitions(
    suppress_multiple=True)] == [
    "d0", "a0.d1", "a0.a1.t0.c", "a0.a1.t0.t1.x", "a0.d2"]
  parameters = phil.parse(input_string="""\
s {
  a=0
}
t.a=0
""")
  check_get(parameters, path="s", expected_out="""\
s {
  a = 0
}
""")
  check_get(parameters, path="s.a", expected_out="""\
a = 0
""")
  check_get(parameters, path="t", expected_out="""\
t.a = 0
""")
  check_get(parameters, path="t.a", expected_out="""\
a = 0
""")
  parameters = phil.parse(input_string="""\
s {
  a {
    b=0
    c.d=1
  }
}
t.a.b=0
t.a.c.d=1
""")
  out = StringIO()
  parameters.show(out=out)
  assert out.getvalue() == """\
s {
  a {
    b = 0
    c.d = 1
  }
}
t.a.b = 0
t.a.c.d = 1
"""
  check_get(parameters, path="s", expected_out="""\
s {
  a {
    b = 0
    c.d = 1
  }
}
""")
  check_get(parameters, path="s.a", expected_out="""\
a {
  b = 0
  c.d = 1
}
""")
  check_get(parameters, path="s.a.b", expected_out="""\
b = 0
""")
  check_get(parameters, path="s.a.c", expected_out="""\
c.d = 1
""")
  check_get(parameters, path="t", expected_out="""\
t.a.b = 0
t.a.c.d = 1
""")
  check_get(parameters, path="t.a", expected_out="""\
a.b = 0
a.c.d = 1
""")
  check_get(parameters, path="t.a.b", expected_out="""\
b = 0
""")
  check_get(parameters, path="t.a.c", expected_out="""\
c.d = 1
""")
  parameters = phil.parse(input_string="""\
a { b { } c { } }
x.y { }
x.z { }
""")
  check_get(parameters, path="a", expected_out="""\
a {
  b {
  }
  c {
  }
}
""")
  check_get(parameters, path="a.b", expected_out="""\
b {
}
""")
  check_get(parameters, path="a.c", expected_out="""\
c {
}
""")
  check_get(parameters, path="x", expected_out="""\
x.y {
}
x.z {
}
""")
  check_get(parameters, path="x.y", expected_out="""\
y {
}
""")
  check_get(parameters, path="x.z", expected_out="""\
z {
}
""")
  parameters = phil.parse(input_string="""\
x.y { a.b { d.e=0
} }
""")
  check_get(parameters, path="x", expected_out="""\
x.y {
  a.b {
    d.e = 0
  }
}
""")
  check_get(parameters, path="x.y", expected_out="""\
y {
  a.b {
    d.e = 0
  }
}
""")
  check_get(parameters, path="x.y.a", expected_out="""\
a.b {
  d.e = 0
}
""")
  check_get(parameters, path="x.y.a.b", expected_out="""\
b {
  d.e = 0
}
""")
  check_get(parameters, path="x.y.a.b.d", expected_out="""\
d.e = 0
""")
  check_get(parameters, path="x.y.a.b.d.e", expected_out="""\
e = 0
""")

def check_resolve_variables(parameters, path, expected_out=None, n_matches=1):
  matches = parameters.get(path=path, with_substitution=False)
  assert len(matches.objects) == n_matches
  result = matches.objects[0].resolve_variables().as_str()
  if (expected_out is None):
    print '  check_resolve_variables(parameters, "%s", "%s")' % (
      path, result.replace("\n", "\\n"))
  elif (result != expected_out):
    raise AssertionError('"%s" != "%s"' % (result, expected_out))

def exercise_get_with_substitution():
  parameters = phil.parse(input_string="""\
a=b
c = d   e   2
""")
  check_get_sub(parameters, path="a", expected_out="a = b\n")
  check_get_sub(parameters, path="c", expected_out="c = d e 2\n")
  check_resolve_variables(parameters, "a", "a = b\n")
  check_resolve_variables(parameters, "c", "c = d e 2\n")
  parameters = phil.parse(input_string="""\
a=1
b=$a
c=$b
d=2
e=${.d} $c
s {
  i=10
  j=$i
  k=$j
  l=20
  m=$l $k
  n=$a
  t {
    x=$l
    y=$k
    z=$n
    n=$a $j $y
  }
  d=9
  e=${d}
  f=${.d}
  g=${.s.t.x}
}
d=x
f=${s.i}
g=${s.a}
s {
  t {
    x=30
    a=40
  }
}
h=${s.t.x}
i=${s.t.a}
j=${s.t.n}
k=${.s.t.x}
l=$s
""")
  check_resolve_variables(parameters, "a", "a = 1\n")
  check_resolve_variables(parameters, "b", "b = 1\n")
  check_resolve_variables(parameters, "c", "c = 1\n")
  check_resolve_variables(parameters, "d", "d = 2\n", n_matches=2)
  check_resolve_variables(parameters, "e", "e = 2 1\n")
  check_resolve_variables(parameters, "s.i", "i = 10\n")
  check_resolve_variables(parameters, "s.j", "j = 10\n")
  check_resolve_variables(parameters, "s.k", "k = 10\n")
  check_resolve_variables(parameters, "s.m", "m = 20 10\n")
  check_resolve_variables(parameters, "s.n", "n = 1\n")
  check_resolve_variables(parameters, "s.t.x", "x = 20\n", n_matches=2)
  check_resolve_variables(parameters, "s.t.y", "y = 10\n")
  check_resolve_variables(parameters, "s.t.z", "z = 1\n")
  check_resolve_variables(parameters, "s.t.n", "n = 1 10 10\n")
  check_resolve_variables(parameters, "f", "f = 10\n")
  try: check_resolve_variables(parameters, "g")
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $s.a (input line 26)'
  else: raise RuntimeError("Exception expected.")
  check_resolve_variables(parameters, "h", "h = 30\n")
  check_resolve_variables(parameters, "i", "i = 40\n")
  check_resolve_variables(parameters, "j", "j = 1 10 10\n")
  check_resolve_variables(parameters, "k", "k = 30\n")
  check_resolve_variables(parameters, "s.e", "e = 9\n")
  check_resolve_variables(parameters, "s.f", "f = 2\n")
  check_resolve_variables(parameters, "s.g", "g = 20\n")
  try: check_resolve_variables(parameters, "l")
  except RuntimeError, e:
    assert str(e) == 'Not a definition: $s (input line 37)'
  else: raise RuntimeError("Exception expected.")
  check_resolve_variables(parameters, "h", "h = 30\n")
  #
  parameters = phil.parse(input_string="""\
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
k=x
k=y
l=$k
k=z
m=$k
""")
  assert len(parameters.get(path="a", with_substitution=False).objects) == 1
  diags = {
    "a": "Undefined variable: $a (input line 1)",
    "b": "Undefined variable: $c (input line 2)",
    "c": "Undefined variable: $c (input line 2)",
    "d": "Undefined variable: $e (input line 4)",
    "e": "Undefined variable: $f (input line 5)",
    "f": "Undefined variable: $e (input line 4)",
    "g": "Undefined variable: $e (input line 4)"}
  for path in "abcdefg":
    try: parameters.get(path=path)
    except RuntimeError, e:
      if (str(e) != diags[path]):
        raise AssertionError("%s != %s" % (repr(str(e)), repr(diags[path])))
    else: raise RuntimeError("Exception expected.")
  try: parameters.get(path="h")
  except RuntimeError, e:
    assert str(e) == "Undefined variable: $_X_Y_Z_ (input line 8)"
  else: raise RuntimeError("Exception expected.")
  os.environ["_X_Y_Z_"] = "xyz"
  check_get_sub(parameters, path="h", expected_out='h = "xyz"\n')
  assert parameters.get(path="h").objects[0].words[0].where_str() \
      == ' (environment: "_X_Y_Z_")'
  check_get_sub(parameters, path="j", expected_out='j = 1\n')
  check_get_sub(parameters, path="l", expected_out='l = y\n')
  check_get_sub(parameters, path="m", expected_out='m = z\n')
  parameters = phil.parse(input_string="""\
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
  #
  parameters = phil.parse(input_string="""\
v=x
w=y
s {
  a=$v
  s {
    a=$w
    b=${s.a}
    c=${s.s.a}
    d=${.s.a}
  }
}
""")
  check_get_sub(parameters, path="s", expected_out="""\
s {
  a = x
  s {
    a = y
    b = y
    c = y
    d = x
  }
}
""")
  #
  parameters = phil.parse(input_string="""\
a=$b
#b=y
""")
  try: parameters.get(path="a")
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $b (input line 1)'
  else: raise RuntimeError("Exception expected.")

def exercise_include():
  print >> open("tmp1.params", "w"), """\
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
  parameters = phil.parse(
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
b = y
d = $z
s = 1
"""
  out = StringIO()
  parameters.unique().show(out=out)
  assert out.getvalue() == """\
tmp2 = tmp2.params
a = x
r = 0
c = z
b = y
d = $z
s = 1
"""
  try: parameters.get(path="d")
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $z (file "tmp3.params", line 3)'
  else: raise RuntimeError("Exception expected.")
  try: os.makedirs("tmp")
  except OSError: pass
  #
  print >> open("tmp1.params", "w"), """\
include tmp3.params
"""
  print >> open("tmp2.params", "w"), """\
include tmp1.params
"""
  print >> open("tmp3.params", "w"), """\
include tmp2.params
"""
  try: parameters = phil.parse(
    file_name="tmp1.params",
    process_includes=True)
  except RuntimeError, e:
    assert str(e).startswith("Include dependency cycle: ")
    assert len(str(e).split(",")) == 4
  else: raise RuntimeError("Exception expected.")
  #
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
z=3
"""
  parameters = phil.parse(
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
  print >> open("tmp4.params", "w"), """\
a=1
include tmp1.params
s {
  a=2
  include tmp1.params
  z=2
}
s {
  a=3
  include tmp1.params
  z=3
}
z=1
"""
  parameters = phil.parse(
    file_name="tmp4.params",
    process_includes=True)
  out = StringIO()
  parameters.show(out=out)
  assert out.getvalue() == """\
a = 1
a = 0
b = 1
c = 2
z = 3
y = 2
x = 1
s {
  a = 2
  a = 0
  b = 1
  c = 2
  z = 3
  y = 2
  x = 1
  z = 2
}
s {
  a = 3
  a = 0
  b = 1
  c = 2
  z = 3
  y = 2
  x = 1
  z = 3
}
z = 1
"""
  out = StringIO()
  parameters.unique().show(out=out)
  assert out.getvalue() == """\
a = 0
b = 1
c = 2
y = 2
x = 1
s {
  a = 0
  b = 1
  c = 2
  y = 2
  x = 1
  z = 3
}
z = 1
"""

def exercise_fetch():
  master = phil.parse(input_string="""\
a=None
  .expert_level=1
""")
  source = phil.parse(input_string="""\
a=1
a=2
""")
  out = StringIO()
  master.fetch(source).show(out=out, attributes_level=2)
  assert out.getvalue() == """\
a = 2
  .expert_level = 1
"""
  source = phil.parse(input_string="""\
a=1
#a=2
""")
  out = StringIO()
  master.fetch(source).show(out=out, attributes_level=2)
  assert out.getvalue() == """\
a = 1
  .expert_level = 1
"""
  #
  master_plain = phil.parse(input_string="""\
s
  .expert_level=1
{
  a=None
    .expert_level=2
}
""")
  master_optional = phil.parse(input_string="""\
s
  .optional=yes
  .expert_level=1
{
  a=None
    .expert_level=2
}
""")
  master_multiple = phil.parse(input_string="""\
s
  .multiple=1
  .expert_level=1
{
  a=None
    .expert_level=2
}
""")
  master_optional_multiple = phil.parse(input_string="""\
s
  .multiple=1
  .optional=yes
  .expert_level=1
{
  a=None
    .expert_level=2
}
""")
  source = phil.parse(input_string="""\
s {
  a=1
  a=2
}
""")
  out = StringIO()
  master_plain.fetch(source).show(out=out, attributes_level=2)
  assert out.getvalue() == """\
s
  .expert_level = 1
{
  a = 2
    .expert_level = 2
}
"""
  source = phil.parse(input_string="""\
s {
  a=1
  #a=2
}
""")
  out = StringIO()
  master_plain.fetch(source).show(out=out, attributes_level=2)
  assert out.getvalue() == """\
s
  .expert_level = 1
{
  a = 1
    .expert_level = 2
}
"""
  source = phil.parse(input_string="""\
s {
  a=1
}
s {
  a=2
}
""")
  out = StringIO()
  master_plain.fetch(source).show(out=out, attributes_level=2)
  assert out.getvalue() == """\
s
  .expert_level = 1
{
  a = 2
    .expert_level = 2
}
"""
  out = StringIO()
  master_multiple.fetch(source).show(out=out)
  assert out.getvalue() == """\
s {
  a = 1
}
s {
  a = 2
}
"""
  source = phil.parse(input_string="""\
s {
  a=1
}
s {
  #a=2
}
""")
  for master in [master_plain, master_optional]:
    out = StringIO()
    master.fetch(source).show(out=out)
    assert out.getvalue() == """\
s {
  a = 1
}
"""
  for master in [master_multiple, master_optional_multiple]:
    out = StringIO()
    master.fetch(source).show(out=out)
    assert out.getvalue() == """\
s {
  a = 1
}
s {
  a = None
}
"""
  source = phil.parse(input_string="""\
#s {
  a=1
}
s {
  #a=2
}
""")
  for master in [master_plain, master_optional, master_multiple]:
    out = StringIO()
    master.fetch(source).show(out=out)
    assert out.getvalue() == """\
s {
  a = None
}
"""
  out = StringIO()
  master_optional_multiple.fetch(source).show(out=out)
  assert out.getvalue() == """\
#s {
  a = None
}
"""
  source = phil.parse(input_string="")
  out = StringIO()
  master_optional_multiple.fetch(source).show(out=out)
  assert out.getvalue() == """\
#s {
  a = None
}
"""
  #
  master = phil.parse(input_string="""\
a=None
b {}
c { a=None
}
d { a {} }
""", source_info="master")
  source = phil.parse(input_string="a { }", source_info="source")
  try: master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == 'Incompatible parameter objects "a" ' \
      '(master, line 1) and "a" (source, line 1)'
  else: raise RuntimeError("Exception expected.")
  source = phil.parse(input_string="b=None")
  try: master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == 'Incompatible parameter objects "b" ' \
      '(master, line 2) and "b" (input line 1)'
  else: raise RuntimeError("Exception expected.")
  source = phil.parse(input_string="c { a { } }")
  try: master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == 'Incompatible parameter objects "a" ' \
      '(master, line 3) and "a" (input line 1)'
  else: raise RuntimeError("Exception expected.")
  source = phil.parse(input_string="d { a=None\n}")
  try: master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == 'Incompatible parameter objects "a" ' \
      '(master, line 5) and "a" (input line 1)'
  else: raise RuntimeError("Exception expected.")
  #
  master = phil.parse(input_string="""\
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
  source = phil.parse(input_string="""\
a=7
v
{
  x=y
  y=3
}
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
  fetched = master.fetch(source=source)
  out = StringIO()
  fetched.show(out=out, attributes_level=2)
  assert out.getvalue() == """\
a = 9
  .expert_level = 1
b = None
  .expert_level = 2
c {
  a = y
    .expert_level = 3
  b = 1
    .expert_level = 4
  c
    .expert_level = 5
  {
    x = 1
      .expert_level = 6
    y = 3
      .type = int
      .expert_level = 7
    t
      .expert_level = 8
    {
      r = x
        .help = help
    }
  }
}
t
  .expert_level = 9
{
  a = 1
    .expert_level = 10
  b = 3
    .expert_level = 11
}
u = a b *c
  .type = choice
"""
  master = phil.parse(input_string="""\
c=a *b c
  .type=choice
""")
  for choice in ["a", "b", "c"]:
    for stat in ["", "*"]:
      source = phil.parse(input_string="c="+choice)
      out = StringIO()
      master.fetch(source=source).show(out=out)
      if (choice == "a"):
        assert out.getvalue() == "c = *a b c\n"
      elif (choice == "b"):
        assert out.getvalue() == "c = a *b c\n"
      else:
        assert out.getvalue() == "c = a b *c\n"
  source = phil.parse(input_string="""\
c=a *d c
""")
  try: fetched = master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == "Not a possible choice: *d (input line 1)"
  else: raise RuntimeError("Exception expected.")
  #
  master = phil.parse(input_string="""\
v
  .multiple=true
{
  x=None
}
c
  .multiple=true
{
  a=None
}
""")
  source = phil.parse(input_string="""\
v {
  x=y
}
c {
  a=${v.x}
}
v {
  x=z
}
c {
  a=${v.x}
}
""")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert out.getvalue() == """\
v {
  x = y
}
v {
  x = z
}
c {
  a = y
}
c {
  a = z
}
"""
  parameters = phil.parse(input_string="""\
c {
  a=${v.x}
}
v {
  x=y
}
""")
  try: parameters.fetch(source=parameters)
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $v.x (input line 2)'
  else: raise RuntimeError("Exception expected.")
  #
  master = phil.parse(input_string="""\
s {
  a=None
}
""")
  source = phil.parse(input_string="s.a=x")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert out.getvalue() == """\
s {
  a = x
}
"""
  master = phil.parse(input_string="""\
s {
  t {
    a=None
  }
}
""")
  source = phil.parse(input_string="s.t.a=x")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert out.getvalue() == """\
s {
  t {
    a = x
  }
}
"""
  master = phil.parse(input_string="""\
s.t {
  a=None
}
""")
  source = phil.parse(input_string="s.t.a=x")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert out.getvalue() == """\
s.t {
  a = x
}
"""
  master = phil.parse(input_string="""\
s.t {
  a.b=None
}
""")
  source = phil.parse(input_string="s.t.a.b=x")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert out.getvalue() == """\
s.t {
  a.b = x
}
"""
  master = phil.parse(input_string="""\
s.t.u {
  a.b.c=None
}
""")
  source = phil.parse(input_string="""\
v.w {
  p.q=z
}
s.t.u.a.b.c=${v.w.p.q}
""")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert out.getvalue() == """\
s.t.u {
  a.b.c = z
}
"""
  for input_string in ["s{t.u.a.b.c=x\n}",
                       "s.t{u.a.b.c=x\n}",
                       "s.t.u{a.b.c=x\n}",
                       "s.t.u.a{b.c=x\n}",
                       "s.t.u.a.b{c=x\n}",
                       "s{t{u.a.b.c=x\n}}",
                       "s{t.u{a.b.c=x\n}}",
                       "s{t.u.a{b.c=x\n}}"]:
    source = phil.parse(input_string=input_string)
    out = StringIO()
    master.fetch(source=source).show(out=out)
    assert out.getvalue() == """\
s.t.u {
  a.b.c = x
}
"""
  #
  master = phil.parse(input_string="""\
s.t
  .multiple=true
{
  a.b=None
}
""")
  source = phil.parse(input_string="""\
s.t.a.b=x
s.t.a.b=y
s{t.a.b=z
  t.a.b=q
}
s.t{a.b=r
    a.b=w
}
s.t.a{b=e
      b=f
}
s.t.a{b=e
      #b=f
}
""")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert out.getvalue() == """\
s.t {
  a.b = x
}
s.t {
  a.b = y
}
s.t {
  a.b = z
}
s.t {
  a.b = q
}
s.t {
  a.b = w
}
s.t {
  a.b = f
}
s.t {
  a.b = e
}
"""
  master = phil.parse(input_string="""\
s {
  a {
    f=None
  }
  b {
    f=None
    g=None
  }
}
""")
  source_af = phil.parse(input_string="s.a.f=1")
  source_bf = phil.parse(input_string="s.b.f=2")
  source_bg = phil.parse(input_string="s.b.g=3")
  assert master.fetch(sources=[source_af, source_bf]).as_str() == """\
s {
  a {
    f = 1
  }
  b {
    f = 2
    g = None
  }
}
"""
  assert master.fetch(sources=[source_bg, source_af]).as_str() == """\
s {
  a {
    f = 1
  }
  b {
    f = None
    g = 3
  }
}
"""
  assert master.fetch(sources=[source_bg, source_af, source_bf]).as_str()=="""\
s {
  a {
    f = 1
  }
  b {
    f = 2
    g = 3
  }
}
"""

def exercise_extract():
  parameters = phil.parse(input_string="""\
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
    .optional=no
  f=a *b c
    .type=choice(multi = True)
  f=a *b *c
    .type=choice(multi = True)
  f=a b c
    .type=choice(multi =True)
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
  n=none
    .type=None
  m=Plain "%^ &*"
    .type=None
  s0=none
    .type=strings
  s1=" a " b
    .type=strings
  v=none
    .type=words
  w=plain "% ^&*"
    .type=words
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
  assert parameters.get(path="group.e").objects[0].type \
      is parameters.get(path="group.e").objects[1].type
  assert parameters.get(path="group.e").objects[0].type \
      is parameters.get(path="group.e").objects[2].type
  assert parameters.get(path="group.f",
    with_substitution=False).objects[0].extract() == ["b"]
  assert parameters.get(path="group.f",
    with_substitution=False).objects[1].extract() == ["b", "c"]
  assert parameters.get(path="group.f",
    with_substitution=False).objects[2].extract() == []
  assert parameters.get(path="group.f").objects[0].type \
      is parameters.get(path="group.f").objects[1].type
  assert parameters.get(path="group.f").objects[0].type \
  is not parameters.get(path="group.f").objects[2].type
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
  assert parameters.get(path="group.n",
    with_substitution=False).objects[0].extract() is None
  assert parameters.get(path="group.m",
    with_substitution=False).objects[0].extract() == ["Plain", "%^ &*"]
  assert parameters.get(path="group.s0",
    with_substitution=False).objects[0].extract() is None
  assert parameters.get(path="group.s1",
    with_substitution=False).objects[0].extract() == [" a ", "b"]
  assert parameters.get(path="group.v",
    with_substitution=False).objects[0].extract() is None
  assert [word.value for word in parameters.get(path="group.w",
    with_substitution=False).objects[0].extract()] == ["plain", "% ^&*"]
  definition = parameters.get(path="group.a",
    with_substitution=False).objects[0]
  parameters = phil.parse(input_string="""\
group {
  a=yes
    .type=bool
  b=13
    .type=int
}
""")
  group = parameters.get(path="group",
    with_substitution=False).objects[0].extract()
  assert group.a is True
  assert group.b == 13
  #
  parameters = phil.parse(input_string="""\
c.d=5
  .type=int
c.e=6
  .type=int
c {
  f=7
    .type=int
}
""")
  extracted = parameters.extract()
  assert extracted.c.d == 5
  assert extracted.c.e == 6
  assert extracted.c.f == 7
  parameters = phil.parse(input_string="""\
c.d=5
  .type=int
c.e=6
  .type=int
c {
  d=7
    .type=int
  f=8
    .type=int
}
""")
  extracted = parameters.extract()
  assert extracted.c.d == 7
  assert extracted.c.e == 6
  assert extracted.c.f == 8
  #
  parameters = phil.parse(input_string="""\
a=1
  .type=int
b {
  a=2
    .type=int
  b {
    a=3
      .type=int
  }
}
c.a {
  b=4
    .type=int
  c.d=5
    .type=int
  c.e=6
    .type=int
}
""")
  extracted = parameters.extract()
  assert extracted.a == 1
  assert extracted.b.a == 2
  assert extracted.b.b.a == 3
  assert extracted.c.a.b == 4
  assert extracted.c.a.c.d == 5
  assert extracted.c.a.c.e == 6
  #
  parameters = phil.parse(input_string="""\
a=1
 .type=int
a=2
 .type=int
b {
  a=3
    .type=int
  a=4
    .type=int
}
b {
  b=5
    .type=int
  b=6
    .type=int
}
""")
  extracted = parameters.extract()
  assert extracted.a == 2
  assert extracted.b.b == 6

def exercise_format():
  parameters = phil.parse(input_string="""\
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
    .type=choice(multi=True)
  g="/var/tmp/foo"
    .type=path
  h="var.tmp.foo"
    .type=key
  m=plain "% ^&*"
    .type=None
  s0=none
    .type=strings
  s1=' a ' b
    .type=strings
  w=plain "% ^&*"
    .type=words
}
""")
  out = StringIO()
  extracted = parameters.fetch(source=parameters).extract()
  parameters.format(extracted).show(out=out)
  assert out.getvalue() == """\
group {
  n = ab "c d" "ef "
  a = True
  b = 13
  c = 1.3
  d = "abc def ghi"
  e = a *b c
  f = a *b *c
  g = "/var/tmp/foo"
  h = "var.tmp.foo"
  m = plain "% ^&*"
  s0 = None
  s1 = " a " b
  w = plain "% ^&*"
}
"""
  #
  master = phil.parse(input_string="""\
a=None
 .type=int
 .multiple=False
""")
  custom = phil.parse(input_string="""\
a=1
a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == 2
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == "a = 2\n"
  custom = phil.parse(input_string="""\
a=1
#a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == 1
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == "a = 1\n"
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == None
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == "a = None\n"
  #
  master = phil.parse(input_string="""\
a=None
 .type=int
 .multiple=True
""")
  custom = phil.parse(input_string="""\
a=1
a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == [1,2]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
a = 1
a = 2
"""
  custom = phil.parse(input_string="""\
a=1
#a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == [1]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
a = 1
"""
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == [None]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
a = None
"""
  #
  master = phil.parse(input_string="""\
a=None
 .type=int
 .multiple=True
 .optional=True
""")
  custom = phil.parse(input_string="""\
a=1
a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == [1,2]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
a = 1
a = 2
"""
  custom = phil.parse(input_string="""\
a=1
#a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == [1]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
a = 1
"""
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == ""
  custom = phil.parse(input_string="""\
a=None
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == ""
  #
  master = phil.parse(input_string="""\
s
  .multiple=False
{
  a=None
   .type=int
   .multiple=False
}
""")
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s.a is None
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = None
}
"""
  custom = phil.parse(input_string="""\
s {
  a=1
}
s {
  a=2
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s.a == 2
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = 2
}
"""
  custom = phil.parse(input_string="""\
s {
  a=1
}
s {
  a=2
  a=3
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s.a == 3
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = 3
}
"""
  #
  master = phil.parse(input_string="""\
s
  .multiple=True
{
  a=None
   .type=int
   .multiple=False
}
""")
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s[0].a == None
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = None
}
"""
  custom = phil.parse(input_string="""\
s {
  a=None
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s[0].a == None
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = None
}
"""
  custom = phil.parse(input_string="""\
s {
  a=1
}
s {
  a=2
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s[0].a == 1
  assert extracted.s[1].a == 2
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = 1
}
s {
  a = 2
}
"""
  custom = phil.parse(input_string="""\
s {
  a=1
}
s {
  a=2
  a=3
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s[0].a == 1
  assert extracted.s[1].a == 3
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = 1
}
s {
  a = 3
}
"""
  master = phil.parse(input_string="""\
s
  .multiple=True
  .optional=True
{
  a=None
   .type=int
   .multiple=False
}
""")
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
"""
  custom = phil.parse(input_string="""\
s {
  a=None
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
"""
  custom = phil.parse(input_string="""\
s {
  a=1
}
s {
  a=2
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s[0].a == 1
  assert extracted.s[1].a == 2
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = 1
}
s {
  a = 2
}
"""
  custom = phil.parse(input_string="""\
s {
  a=1
}
s {
  a=2
  a=3
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s[0].a == 1
  assert extracted.s[1].a == 3
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
s {
  a = 1
}
s {
  a = 3
}
"""
  #
  master = phil.parse(input_string="""\
a=None
 .type=int
 .multiple=True
b
  .multiple=True
{
  a=None
    .type=int
    .multiple=True
  b=None
    .type=int
    .multiple=True
    .optional=True
  c=None
    .type=int
}
c
  .multiple=True
  .optional=True
{
  a=None
}
d
  .multiple=True
  .optional=False
{
  a=None
    .type=int
}
""")
  custom = phil.parse(input_string="""\
a=1
a=2
b {
  a=3
  a=4
  c=10
  c=20
}
b {
  b=5
  b=6
}
c {
  a=None
}
""")
  fetched = master.fetch(custom)
  extracted = fetched.extract()
  assert extracted.a == [1,2]
  assert extracted.b[0].a == [3,4]
  assert extracted.b[0].b == []
  assert extracted.b[0].c == 20
  assert extracted.b[1].a == [None]
  assert extracted.b[1].b == [5,6]
  assert extracted.b[1].c is None
  assert extracted.c == []
  assert extracted.d[0].a is None
  out = StringIO()
  master.format(extracted).show(out=out)
  assert out.getvalue() == """\
a = 1
a = 2
b {
  a = 3
  a = 4
  c = 20
}
b {
  a = None
  b = 5
  b = 6
  c = None
}
d {
  a = None
}
"""
  #
  params = phil.parse(input_string="""\
a=1
  .type=int
s {
  b=2
    .type=int
}
""")
  orig = params.extract()
  assert orig.a == 1
  assert orig.s.b == 2
  clone = params.clone(python_object=orig)
  assert clone.a == 1
  assert clone.s.b == 2
  clone.a = 10
  clone.s.b = 20
  assert orig.a == 1
  assert orig.s.b == 2

class foo1_converters:

  def __init__(self, bar=None):
    if (bar is None):
      raise RuntimeError("foo1 problem")

class foo2_converters:

  def __str__(self): return "foo2"

  def __init__(self, bar=None):
    if (bar is not None):
      raise RuntimeError("foo2 problem")

foo_converter_registry = dict(phil.default_converter_registry)
foo_converter_registry["foo1"] = foo1_converters
foo_converter_registry["foo2"] = foo2_converters

def exercise_type_constructors():
  params = phil.parse(
    input_string="""\
a=None
  .type=foo2
b=None
  .type=foo2(bar=None)
c=None
  .type=foo1(bar=0)
""",
    converter_registry=foo_converter_registry)
  try: params.get(path="a").objects[0].extract()
  except RuntimeError, e:
    assert str(e) == \
      '.type=foo2 does not have a from_words method (input line 1):' \
      " foo2_converters instance has no attribute 'from_words'"
  else: raise RuntimeError("Exception expected.")
  try: params.get(path="a").objects[0].format(python_object=0)
  except RuntimeError, e:
    assert str(e) == \
      '.type=foo2 does not have an as_words method (input line 1):' \
      " foo2_converters instance has no attribute 'as_words'"
  else: raise RuntimeError("Exception expected.")
  try: phil.parse(
    input_string="""\
a=None
  .type=foo1
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e) == \
      'Error constructing definition type "foo1":' \
      ' foo1 problem (input line 2)'
  else: raise RuntimeError("Exception expected.")
  try: phil.parse(
    input_string="""\
a=None
  .type=foo2(bar=1)
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e) == \
      'Error constructing definition type "foo2(bar=1)":' \
      ' foo2 problem (input line 2)'
  else: raise RuntimeError("Exception expected.")
  try: phil.parse(
    input_string="""\
a=None
  .type=foo2(foo=1)
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e) == \
      'Error constructing definition type "foo2(foo=1)":' \
      " __init__() got an unexpected keyword argument 'foo' (input line 2)"
  else: raise RuntimeError("Exception expected.")
  try: phil.parse(
    input_string="""\
a=None
  .type=foo2(foo=1
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e) == \
      'Error constructing definition type "foo2(foo=1":' \
      " unexpected EOF while parsing (line 1) (input line 2)"
  else: raise RuntimeError("Exception expected.")

def exercise_command_line():
  master_params = phil.parse(input_string="""\
foo {
  min=0
  max=10
  index=3
  limit=6
}
bar {
  max=5
  sub {
    limit=8
  }
}
""")
  itpr_bar = phil.command_line.argument_interpreter(
    master_params=master_params,
    home_scope="bar")
  itpr_neutral = phil.command_line.argument_interpreter(
    master_params=master_params)
  for itpr in [itpr_bar, itpr_neutral]:
    itpr.process(arg="foo.limit=4\nbar.max=2").as_str() == """\
foo.limit = 4
bar.max = 2
"""
  assert itpr_bar.process(arg="max=5").as_str() == "bar.max = 5\n"
  try: assert itpr_neutral.process(arg="max=5")
  except UserError, e:
    assert str(e) == """\
Ambiguous parameter definition: max = 5
Best matches:
  foo.max
  bar.max"""
  else: raise RuntimeError("Exception expected.")
  assert itpr_bar.process(arg="limit=0").as_str() == "bar.sub.limit = 0\n"
  assert itpr_bar.process(arg="imit=0").as_str() == "bar.sub.limit = 0\n"
  for itpr in [itpr_bar, itpr_neutral]:
    assert itpr.process(arg="index=0").as_str() == "foo.index = 0\n"
    assert itpr.process(arg="ndex=0").as_str() == "foo.index = 0\n"
  try: itpr_bar.process(arg="xyz=")
  except UserError, e:
    assert str(e) == """\
Error interpreting command line argument as parameter definition:
  "xyz="
  Missing value for xyz (command line argument, line 1)"""
  else: raise RuntimeError("Exception expected.")
  try: itpr_bar.process(arg="xyz=8")
  except UserError, e:
    assert str(e) == "Unknown command line parameter definition: xyz = 8"
  else: raise RuntimeError("Exception expected.")
  try: itpr_bar.process(arg="  ")
  except UserError, e:
    assert str(e) == 'Command line parameter definition has no effect: "  "'
  else: raise RuntimeError("Exception expected.")
  itpr = phil.command_line.argument_interpreter(
    master_params=master_params,
    argument_description="")
  try: itpr.process(arg="bar {}")
  except UserError, e:
    assert str(e) == 'Parameter definition has no effect: "bar {}"'
  else: raise RuntimeError("Exception expected.")

def exercise():
  exercise_parse_and_show()
  exercise_syntax_errors()
  exercise_deepcopy()
  exercise_get_without_substitution()
  exercise_nested()
  exercise_get_with_substitution()
  exercise_include()
  exercise_fetch()
  exercise_extract()
  exercise_format()
  exercise_type_constructors()
  exercise_command_line()
  print "OK"

if (__name__ == "__main__"):
  exercise()
