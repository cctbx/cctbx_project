from __future__ import division
from libtbx import phil
import libtbx.phil
from libtbx.utils import Sorry
from libtbx.test_utils import \
  Exception_expected, show_diff, anchored_block_show_diff
from libtbx import Auto
import warnings
from cStringIO import StringIO
import copy
import sys, os

def exercise_string_quote_and_tokenize():
  n_ok = [0]
  for quote_token in ["'", '"', "'''", '"""']:
    def check(string):
      input_string = phil.tokenizer.quote_python_str(
        quote_token=quote_token, string=string)
      words = phil.tokenize_value_literal(
        input_string=input_string, source_info="test_string")
      assert len(words) == 1
      word = words[0]
      assert word.quote_token == quote_token
      assert word.value == string
      assert str(word) == input_string
      n_ok[0] += 1
    ch_set = ["", "'", '"', "\n", "x"]
    for c0 in ch_set:
      check(string=c0)
      for c1 in ch_set:
        c01 = c0 + c1
        check(string=c01)
        for c2 in ch_set:
          c012 = c01 + c2
          check(string=c012)
          for c3 in ch_set:
            c0123 = c012 + c3
            check(string=c0123)
  assert n_ok[0] == 3120
  #
  def tvl(input_string):
    return phil.tokenize_value_literal(
      input_string=input_string, source_info=None)
  words = tvl(input_string='"\t\n\r\x7a\x8F"')
  assert len(words) == 1
  assert words[0].value == "\t\n\rz\x8f"
  assert str(words[0]) == '"\t\n\rz\x8f"'
  value = "".join([chr(_) for _ in xrange(256)])
  words = tvl(input_string=phil.tokenizer.quote_python_str(
    quote_token='"', string=value))
  assert len(words) == 1
  assert words[0].value == value
  words = tvl(input_string=str(words[0]))
  assert len(words) == 1
  assert words[0].value == value
  word = phil.tokenizer.word(value=value, quote_token='"')
  words = tvl(input_string=str(word))
  assert len(words) == 1
  assert words[0].value == value

class recycle(object):

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
    if (expected_out is not None):
      assert not show_diff(self.out, expected_out)
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
  for input_string in ["", "\n", "   \n", "   \t \n \t ", "#", "\t#"]:
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
  .style = None
  .expert_level = None
""")
  input_string = """\
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
  .style = None
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
  .call = None
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
!name # 1 2 3
  !.help="a"
        "b"
        "c"
      # 1 2 3
  .expert_level=1 # x
        # 1 2 3
{
  !a=b
    .type="int" # y
  c=d # e
    !.help="d"
          "e"
          "f"
    .expert_level=2
}
d=a b c # 1 {2} 3 \\
4 {5 6}}
""",
    attributes_level=2,
    expected_out="""\
!name
  .expert_level = 1
{
  !a = b
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
  assert not show_diff(params.as_str(expert_level=0), """\
x = 1
s {
  a = x
}
""")
  assert not show_diff(params.as_str(expert_level=1), """\
x = 1
y = 2
s {
  a = x
}
t {
  a = y
}
""")
  for expert_level in [-1,2,3,None]:
    assert not show_diff(params.as_str(expert_level=expert_level), """\
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
""")
  #
  recycle(
    input_string="s{a=1}t{b=2}", expected_out="""\
s {
  a = 1
}
t {
  b = 2
}
""")
  recycle(
    input_string="s{a=1;b=2  b ;c=3}t{x=4;y=5;}", expected_out="""\
s {
  a = 1
  b = 2 b
  c = 3
}
t {
  x = 4
  y = 5
}
""")
  recycle(
    input_string=
      "a=None;.type=int;s .multiple=True;.optional=False{c=None}"
      "t .style=column;{d=None;.type=float;.optional=True;e=x;.type=str;}",
    attributes_level=2,
    expected_out="""\
a = None
  .type = int
s
  .optional = False
  .multiple = True
{
  c = None
}
t
  .style = column
{
  d = None
    .optional = True
    .type = float
  e = x
    .type = str
}
""")

improper_phil_converters = None

class int_phil_converters(object):

  def __init__(self, factor=1):
    assert int(factor) == factor
    self.factor = factor

  def __str__(self):
    if (self.factor == 1): return "libtbx.phil.tst.int"
    return "libtbx.phil.tst.int(factor=%d)" % self.factor

  def from_words(self, words, master):
    value = phil.int_from_words(words=words, path=master.full_path())
    if (value is None): return value
    return value * self.factor

  def as_words(self, python_object, master):
    if (python_object is None):
      return [phil.tokenizer.word(value="None")]
    return [phil.tokenizer.word(value=str(python_object/self.factor))]

class converter_implementation(int_phil_converters):

  def __str__(self):
    if (self.factor == 1): return "libtbx.phil.tst.converter_factory"
    return "libtbx.phil.tst.converter_factory(factor=%d)" % self.factor

def converter_factory_phil_converters(**args):
  return converter_implementation(**args)

def exercise_import_converters():
  input_string = """\
x1=None
  .type=libtbx.phil.tst.int
y1=3
  .type=libtbx.phil.tst.int
x2=None
  .type=libtbx.phil.tst.int(factor=2)
y2=3
  .type=libtbx.phil.tst.int( factor = 2 )
z3=4
  .type=libtbx.phil.tst.converter_factory( factor= 3 )
"""
  r = recycle(input_string=input_string, expected_out="""\
x1 = None
  .type = libtbx.phil.tst.int
y1 = 3
  .type = libtbx.phil.tst.int
x2 = None
  .type = libtbx.phil.tst.int(factor=2)
y2 = 3
  .type = libtbx.phil.tst.int(factor=2)
z3 = 4
  .type = libtbx.phil.tst.converter_factory(factor=3)
""",
    attributes_level=2)
  params = r.parameters.extract()
  assert params.x1 is None
  assert params.y1 == 3
  assert params.x2 is None
  assert params.y2 == 6
  assert params.z3 == 12
  #
  try: phil.parse("""\
x=None
  .type=foo.int
""")
  except RuntimeError, e:
    assert str(e) == 'Unexpected definition type: "foo.int" (input line 2)'
  else: raise Exception_expected
  #
  try: phil.parse("""\
x=None
  .type=libtbx.phil.tst.none
""")
  except AttributeError, e:
    assert str(e) == '.type=libtbx.phil.tst.none: object' \
      ' "none_phil_converters" not found in module "libtbx.phil.tst"' \
      ' (input line 2)'
  else: raise Exception_expected
  #
  try: phil.parse("""\
x=None
  .type=libtbx.phil.tst.improper
""")
  except TypeError, e:
    assert str(e) == '"libtbx.phil.tst.improper_phil_converters" is not' \
      ' a callable Python object (input line 2)'
  else: raise Exception_expected
  #
  try: phil.parse("""\
x=None
  .type=libtbx.phil.tst.int(factor=1.5)
""")
  except RuntimeError, e:
    e = str(e)
    assert e.startswith(
      'Error constructing definition type "libtbx.phil.tst.int(factor=1.5)":'
      ' AssertionError: ')
    assert e.endswith(' (input line 2)')
  else: raise Exception_expected
  #
  try: phil.parse("""\
x=None
  .type=libtbx.phil.tst.int(a=1=2)
""")
  except RuntimeError, e:
    assert str(e) == 'Error evaluating definition type "libtbx.phil.tst.int' \
      '(a=1=2)": SyntaxError: invalid syntax (line 1) (input line 2)'
  else: raise Exception_expected

def test_exception(input_string, exception_string=None):
  try: phil.parse(input_string=input_string)
  except KeyboardInterrupt: raise
  except Exception, e:
    if (exception_string is None or str(e) != exception_string):
      print str(e)
      if (exception_string is not None):
        print exception_string
    if (exception_string is not None):
      assert str(e) == exception_string
  else: raise Exception_expected

def exercise_syntax_errors():
  test_exception("'a'",
    """Unquoted word expected, found 'a' (input line 1)""")
  test_exception("a",
    'Unexpected end of input.')
  test_exception("a=\nb",
    'Missing value for a (input line 1)')
  test_exception("x=;",
    'Missing value for x (input line 1)')
  test_exception("s{z=}",
    'Missing value for z (input line 1)')
  test_exception("s{y={",
    'Missing value for y (input line 1)')
  test_exception("s{y=1{",
    'Syntax error: unexpected "{" (input line 1)')
  test_exception("s{y=#",
    'Missing value for y (input line 1)')
  test_exception("a b",
    'Syntax error: expected "=", found "b" (input line 1)')
  test_exception("{}",
    'Syntax error: unexpected "{" (input line 1)')
  test_exception("a {",
    'Syntax error: no matching "}" for "{" at input line 1')
  test_exception(";{}",
    'Syntax error: unexpected ";" (input line 1)')
  test_exception("s{;}",
    'Syntax error: unexpected ";" (input line 1)')
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
  test_exception('__foo__=None',
    'Reserved identifier: "__foo__" (input line 1)')
  test_exception('__foo__ {}',
    'Reserved identifier: "__foo__" (input line 1)')
  test_exception('a=None\n.type=foo',
    'Unexpected definition type: "foo" (input line 2)')
  test_exception('s .multiple=True .optional=False {}',
    'One True or False value expected, .multiple="True .optional=False" found'
    ' (input line 1)')

def exercise_phil_on_off_end():
  assert phil.parse(input_string="#phil __ON__\na=1").as_str() == "a = 1\n"
  assert phil.parse(input_string="#phil __OFF__\na=1").as_str() == ""
  assert phil.parse(input_string="#phil __END__\na=1").as_str() == ""
  assert phil.parse(input_string="#phil __OFF__\n#phil __ON__ \na=1").as_str()\
    == "a = 1\n"
  assert phil.parse(input_string="#phil __END__\n#phil __ON__ \na=1").as_str()\
    == ""
  params = phil.parse(input_string="""\

#phil __OFF__
a b c
#phil __ON__
a=1
#phil __OFF__
c a b
 #phil __ON__
b=2

#phil __ON__
c=3
# phil __OFF__
d=4
#phil __OFF__ b
e=5
#phil __OFF__
f=6
""")
  assert not show_diff(params.as_str(), """\
a = 1
c = 3
d = 4
""")
  assert params.get(path="a").objects[0].words[0].line_number == 5
  assert params.get(path="c").objects[0].words[0].line_number == 12
  assert params.get(path="d").objects[0].words[0].line_number == 14
  try: phil.parse(input_string="#phil __Off__\n")
  except RuntimeError, e:
    assert str(e) == "Unknown: #phil __Off__ (input line 1)"
  else: raise Exception_expected

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
  include file name
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
  !d3=x
}
""",
    expected_out="""\
d0 = 0
a0 {
  d1 = a b c
  include file name
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
  !d3 = x
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
  assert not show_diff(out.getvalue(), """\
s {
  a {
    b = 0
    c.d = 1
  }
}
t.a.b = 0
t.a.c.d = 1
""")
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
e=$(.d) $c
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
  e=$(d)
  f=$(.d)
  g=$(.s.t.x)
}
d=x
f=$(s.i)
g=$(s.a)
s {
  t {
    x=30
    a=40
  }
}
h=$(s.t.x)
i=$(s.t.a)
j=$(s.t.n)
k=$(.s.t.x)
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
  else: raise Exception_expected
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
  else: raise Exception_expected
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
o=$_X_Y_Z_/abc.ext
p="$_X_Y_Z_/abc.ext"
q=$o$_X_Y_Z_/abc.ext$p
r="$o$_X_Y_Z_ abc.ext$p"
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
    else: raise Exception_expected
  try: parameters.get(path="h")
  except RuntimeError, e:
    assert str(e) == "Undefined variable: $_X_Y_Z_ (input line 8)"
  else: raise Exception_expected
  os.environ["_X_Y_Z_"] = "xyz"
  check_get_sub(parameters, path="h", expected_out='h = "xyz"\n')
  assert parameters.get(path="h").objects[0].words[0].where_str() \
      == ' (environment: "_X_Y_Z_")'
  check_get_sub(parameters, path="o", expected_out='o = "xyz/abc.ext"\n')
  check_get_sub(parameters, path="p", expected_out='p = "xyz/abc.ext"\n')
  check_get_sub(parameters, path="q",
    expected_out='q = "xyz/abc.extxyz/abc.extxyz/abc.ext"\n')
  check_get_sub(parameters, path="r",
    expected_out='r = "xyz/abc.extxyz abc.extxyz/abc.ext"\n')
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
f=$(a)
g=abc$(answer)bc
h=abc$(answer)bc 12$answer$(a)56
i=$
j=$(abc
k=$(1bc)
l=$()
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
  else: raise Exception_expected
  try: parameters.get(path="j")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: missing ")": "$(abc" (input line 11)'
  else: raise Exception_expected
  try: parameters.get(path="k")
  except RuntimeError, e:
    assert str(e) == 'Syntax error: improper variable name "$(1bc)"' \
                   + ' (input line 12)'
  else: raise Exception_expected
  try: parameters.get(path="l")
  except RuntimeError, e:
    assert str(e)=='Syntax error: improper variable name "$()" (input line 13)'
  else: raise Exception_expected
  try: parameters.get(path="m")
  except RuntimeError, e:
    assert str(e)=='Syntax error: improper variable name "$@" (input line 14)'
  else: raise Exception_expected
  check_get_sub(parameters, path="n", expected_out="n = '$a'\n")
  #
  parameters = phil.parse(input_string="""\
v=x
w=y
s {
  a=$v
  s {
    a=$w
    b=$(s.a)
    c=$(s.s.a)
    d=$(.s.a)
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
!b=y
""")
  try: parameters.get(path="a")
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $b (input line 1)'
  else: raise Exception_expected

def exercise_full_path():
  params = phil.parse(input_string="""\
d0=0
a0 {
  d1=a
  a1 {
    t0 {
      c=yes
      t1 {
        x=0
      }
    }
  }
  d2=e
}
""")
  for path in ["d0", "a0", "a0.d1", "a0.a1", "a0.a1.t0", "a0.a1.t0.c",
               "a0.a1.t0.t1", "a0.a1.t0.t1.x", "a0.d2"]:
    assert params.get(path).objects[0].full_path() == path

def exercise_include():
  print >> open("tmp1.params", "w"), """\
!include none
a=x
"""
  print >> open("tmp2.params", "w"), """\
b=y
"""
  print >> open("tmp3.params", "w"), """\
c=z
include file tmp2.params
d=$z
"""
  parameters = phil.parse(
    input_string="""\
tmp2=tmp2.params
include file tmp1.params
include file $tmp2
r=0
include file tmp3.params
s=1
""",
    process_includes=True)
  out = StringIO()
  parameters.show(out=out)
  assert not show_diff(out.getvalue(), """\
tmp2 = tmp2.params
!include none
a = x
b = y
r = 0
c = z
b = y
d = $z
s = 1
""")
  out = StringIO()
  parameters.unique().show(out=out)
  assert not show_diff(out.getvalue(), """\
tmp2 = tmp2.params
a = x
r = 0
c = z
b = y
d = $z
s = 1
""")
  try: parameters.get(path="d")
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $z (file "tmp3.params", line 3)'
  else: raise Exception_expected
  try: os.makedirs("tmp")
  except OSError: pass
  #
  print >> open("tmp1.params", "w"), """\
include file tmp3.params
"""
  print >> open("tmp2.params", "w"), """\
include file tmp1.params
"""
  print >> open("tmp3.params", "w"), """\
include file tmp2.params
"""
  try: parameters = phil.parse(
    file_name="tmp1.params",
    process_includes=True)
  except RuntimeError, e:
    assert str(e).startswith("Include dependency cycle: ")
    assert len(str(e).split(",")) == 4
  else: raise Exception_expected
  #
  print >> open("tmp1.params", "w"), """\
a=0
include file tmp/tmp1.params
x=1
"""
  print >> open("tmp/tmp1.params", "w"), """\
b=1
include file tmp2.params
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
  assert not show_diff(out.getvalue(), """\
a = 0
b = 1
c = 2
z = 3
y = 2
x = 1
""")
  print >> open("tmp4.params", "w"), """\
a=1
include file tmp1.params
s {
  a=2
  include file tmp1.params
  z=2
}
s {
  a=3
  include file tmp1.params
  z=3
}
z=1
"""
  parameters = phil.parse(
    file_name="tmp4.params",
    process_includes=True)
  out = StringIO()
  parameters.show(out=out)
  assert not show_diff(out.getvalue(), """\
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
""")
  out = StringIO()
  parameters.unique().show(out=out)
  assert not show_diff(out.getvalue(), """\
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
""")
  #
  try: phil.parse(input_string="""\
include foo
""", process_includes=True)
  except RuntimeError, e:
    assert str(e) \
      == '"include" must be followed by at least two arguments (input line 1)'
  else: raise Exception_expected
  try: phil.parse(input_string="""\
include foo bar
""", process_includes=True)
  except RuntimeError, e:
    assert str(e) == "Unknown include type: foo (input line 1)"
  else: raise Exception_expected
  try: phil.parse(input_string="""\
include file foo foo
""", process_includes=True)
  except RuntimeError, e:
    assert str(e) \
      == '"include file" must be followed exactly one argument (input line 1)'
  else: raise Exception_expected
  #
  params = phil.parse(input_string="""\
include scope libtbx.phil.tst.include_scope_target_1
""", process_includes=True)
  assert not show_diff(params.as_str(attributes_level=1), """\
x = 1
  .help = u
s
  .help = v
{
  y = 2
    .help = w
}
""")
  params = phil.parse(input_string="""\
o {
  include scope libtbx.phil.tst.include_scope_target_1
}
""", process_includes=True)
  assert not show_diff(params.as_str(attributes_level=1), """\
o {
  x = 1
    .help = u
  s
    .help = v
  {
    y = 2
      .help = w
  }
}
""")
  assert not show_diff(params.get('o.s.y').objects[0].full_path(), 'o.s.y')
  params = phil.parse(input_string="""\
include scope libtbx.phil.tst.include_scope_target_1 x
""", process_includes=True)
  assert not show_diff(params.as_str(attributes_level=1), """\
x = 1
  .help = u
""")
  params = phil.parse(input_string="""\
include scope libtbx.phil.tst.include_scope_target_1 s
""", process_includes=True)
  assert not show_diff(params.as_str(attributes_level=1), """\
s
  .help = v
{
  y = 2
    .help = w
}
""")
  params = phil.parse(input_string="""\
a = None
c {
  include scope libtbx.phil.tst.include_scope_target_1 s.y
}
b = None
""", process_includes=True)
  assert not show_diff(params.as_str(attributes_level=1), """\
a = None
c {
  y = 2
    .help = w
}
b = None
""")
  for incl_postfix in ["", "s"]:
    params = phil.parse(input_string="""\
a = None
c {
  include scope libtbx.phil.tst.include_scope_target_2%s
}
b = None
""" % incl_postfix, process_includes=True)
    assert not show_diff(params.as_str(attributes_level=1), """\
a = None
c {
  p = 1
  x = 1
    .help = u
  s
    .help = v
  {
    y = 2
      .help = w
  }
  q = 2
  r {
    y = 2
      .help = w
  }
  x = 3
}
b = None
""")
    extracted = params.extract()
    assert extracted.c.s.__phil_path__() == "c.s"
    try: phil.parse(input_string="""\
include scope foo foo foo
""", process_includes=True)
    except RuntimeError, e:
      assert str(e) \
        == '"include scope" must be followed one or two arguments,' \
           ' i.e. an import path and optionally a phil path (input line 1)'
    else: raise Exception_expected
    try: phil.parse(input_string="""\
include scope libtbx
""", process_includes=True)
    except ValueError, e:
      assert str(e) == 'include scope: import path "libtbx" is too short;' \
        ' target must be a phil scope object or phil string (input line 1)'
    else: raise Exception_expected
    try: phil.parse(input_string="""\
include scope libtbx.phil.t_s_t.include_scope_target_1
""", process_includes=True)
    except ImportError, e:
      assert str(e) \
        == "include scope: no module libtbx.phil.t_s_t (input line 1)"
    else: raise Exception_expected
    try: phil.parse(input_string="""\
include scope libtbx.phil.tst.include_scope_target_none
""", process_includes=True)
    except AttributeError, e:
      assert str(e) == 'include scope: object' \
        ' "include_scope_target_none" not found in module "libtbx.phil.tst"' \
        ' (input line 1)'
    else: raise Exception_expected
    try: phil.parse(input_string="""\
include scope libtbx.phil.tst.include_scope_target_0n
""", process_includes=True)
    except RuntimeError, e:
      assert str(e) == 'include scope: python object' \
        ' "include_scope_target_0n" in module "libtbx.phil.tst"' \
        ' is not a libtbx.phil.scope instance (input line 1)'
    else: raise Exception_expected
    try: phil.parse(input_string="""\
include scope libtbx.phil.tst.include_scope_target_0f
""", process_includes=True)
    except RuntimeError, e:
      assert str(e) == 'include scope: python object' \
        ' "include_scope_target_0f" in module "libtbx.phil.tst"' \
        ' is not a libtbx.phil.scope instance (input line 1)'
    else: raise Exception_expected
    try: phil.parse(input_string="""\
include scope libtbx.phil.tst.include_scope_target_1 t
""", process_includes=True)
    except RuntimeError, e:
      assert str(e) == \
        'include scope: path "t" not found in phil scope object' \
        ' "include_scope_target_1" in module "libtbx.phil.tst" (input line 1)'
    else: raise Exception_expected

include_scope_target_0n = None

include_scope_target_0f = 1.0

include_scope_target_1 = phil.parse("""\
x=1
  .help=u
s
  .help=v
{
  y=2
    .help=w
}
""")

include_scope_target_2s = """\
p=1
include scope libtbx.phil.tst.include_scope_target_1
q=2
r {
  include scope libtbx.phil.tst.include_scope_target_1 s.y
}
x=3
"""

include_scope_target_2 = phil.parse(include_scope_target_2s)

def exercise_fetch():
  master = phil.parse(input_string="""\
a=None
  .expert_level=1
""")
  source = phil.parse(input_string="")
  for f in [master.fetch(source), master.fetch(sources=[]), master.fetch()]:
    assert not show_diff(f.as_str(), """\
a = None
""")
  source = phil.parse(input_string="""\
a=1
a=2
""")
  out = StringIO()
  master.fetch(source).show(out=out, attributes_level=2)
  assert not show_diff(out.getvalue(), """\
a = 2
  .expert_level = 1
""")
  source = phil.parse(input_string="""\
a=1
!a=2
""")
  out = StringIO()
  master.fetch(source).show(out=out, attributes_level=2)
  assert not show_diff(out.getvalue(), """\
a = 1
  .expert_level = 1
""")
  #
  master = phil.parse(input_string="""\
s {
  t {
    v=1
      .type=int
  }
}
s {
  u {
    a=3
      .type=int
  }
}
""")
  source = phil.parse("")
  assert not show_diff(master.fetch(source=source).as_str(), """\
s {
  t {
    v = 1
  }
}
s {
  u {
    a = 3
  }
}
""")
  source = phil.parse("s.t.v=2\ns.u.a=4")
  assert not show_diff(master.fetch(source=source).as_str(), """\
s {
  t {
    v = 2
  }
}
s {
  u {
    a = 4
  }
}
""")
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
  assert not show_diff(out.getvalue(), """\
s
  .expert_level = 1
{
  a = 2
    .expert_level = 2
}
""")
  source = phil.parse(input_string="""\
s {
  a=1
  !a=2
}
""")
  out = StringIO()
  master_plain.fetch(source).show(out=out, attributes_level=2)
  assert not show_diff(out.getvalue(), """\
s
  .expert_level = 1
{
  a = 1
    .expert_level = 2
}
""")
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
  assert not show_diff(out.getvalue(), """\
s
  .expert_level = 1
{
  a = 2
    .expert_level = 2
}
""")
  out = StringIO()
  master_multiple.fetch(source).show(out=out)
  assert not show_diff(out.getvalue(), """\
s {
  a = 1
}
s {
  a = 2
}
""")
  source = phil.parse(input_string="""\
s {
  a=1
}
s {
  !a=2
}
""")
  for master in [master_plain, master_optional]:
    out = StringIO()
    master.fetch(source).show(out=out)
    assert not show_diff(out.getvalue(), """\
s {
  a = 1
}
""")
  for master in [master_multiple, master_optional_multiple]:
    out = StringIO()
    master.fetch(source).show(out=out)
    assert not show_diff(out.getvalue(), """\
s {
  a = 1
}
""")
  source = phil.parse(input_string="""\
!s {
  a=1
}
s {
  !a=2
}
""")
  for master in [master_plain, master_optional,
                 master_multiple, master_optional_multiple]:
    out = StringIO()
    master.fetch(source).show(out=out)
    assert not show_diff(out.getvalue(), """\
s {
  a = None
}
""")
  source = phil.parse(input_string="")
  out = StringIO()
  master_optional_multiple.fetch(source).show(out=out)
  assert not show_diff(out.getvalue(), """\
s {
  a = None
}
""")
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
    assert str(e) == 'Incompatible parameter objects: definition "a" ' \
      '(master, line 1) vs. scope "a" (source, line 1)'
  else: raise Exception_expected
  master.fetch(source=source, skip_incompatible_objects=True)
  source = phil.parse(input_string="b=None")
  try: master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == 'Incompatible parameter objects: scope "b" ' \
      '(master, line 2) vs. definition "b" (input line 1)'
  else: raise Exception_expected
  master.fetch(source=source, skip_incompatible_objects=True)
  source = phil.parse(input_string="c { a { } }")
  try: master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == 'Incompatible parameter objects: definition "a" ' \
      '(master, line 3) vs. scope "a" (input line 1)'
  else: raise Exception_expected
  master.fetch(source=source, skip_incompatible_objects=True)
  source = phil.parse(input_string="d { a=None\n}")
  try: master.fetch(source=source)
  except RuntimeError, e:
    assert str(e) == 'Incompatible parameter objects: scope "a" ' \
      '(master, line 5) vs. definition "a" (input line 1)'
  else: raise Exception_expected
  master.fetch(source=source, skip_incompatible_objects=True)
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
  a=$(v.x)
  b=1
  c
  {
    y=$(v.y)
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
  assert not show_diff(out.getvalue(), """\
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
""")
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
  except Sorry, e:
    assert not show_diff(str(e), """\
Not a possible choice for c: d (input line 1)
  Possible choices are:
    a
    *b
    c""")
  else: raise Exception_expected
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
  a=$(v.x)
}
v {
  x=z
}
c {
  a=$(v.x)
}
""")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert not show_diff(out.getvalue(), """\
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
""")
  parameters = phil.parse(input_string="""\
c {
  a=$(v.x)
}
v {
  x=y
}
""")
  try: parameters.fetch(source=parameters)
  except RuntimeError, e:
    assert str(e) == 'Undefined variable: $v.x (input line 2)'
  else: raise Exception_expected
  #
  master = phil.parse(input_string="""\
s {
  a=None
}
""")
  source = phil.parse(input_string="s.a=x")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert not show_diff(out.getvalue(), """\
s {
  a = x
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  t {
    a = x
  }
}
""")
  master = phil.parse(input_string="""\
s.t {
  a=None
}
""")
  source = phil.parse(input_string="s.t.a=x")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert not show_diff(out.getvalue(), """\
s.t {
  a = x
}
""")
  master = phil.parse(input_string="""\
s.t {
  a.b=None
}
""")
  source = phil.parse(input_string="s.t.a.b=x")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert not show_diff(out.getvalue(), """\
s.t {
  a.b = x
}
""")
  master = phil.parse(input_string="""\
s.t.u {
  a.b.c=None
}
""")
  source = phil.parse(input_string="""\
v.w {
  p.q=z
}
s.t.u.a.b.c=$(v.w.p.q)
""")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert not show_diff(out.getvalue(), """\
s.t.u {
  a.b.c = z
}
""")
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
    assert not show_diff(out.getvalue(), """\
s.t.u {
  a.b.c = x
}
""")
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
      !b=f
}
""")
  out = StringIO()
  master.fetch(source=source).show(out=out)
  assert not show_diff(out.getvalue(), """\
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
""")
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
  assert not show_diff(master.fetch(sources=[source_af, source_bf]).as_str(),
    """\
s {
  a {
    f = 1
  }
  b {
    f = 2
    g = None
  }
}
""")
  assert not show_diff(master.fetch(sources=[source_bg, source_af]).as_str(),
    """\
s {
  a {
    f = 1
  }
  b {
    f = None
    g = 3
  }
}
""")
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
  #
  master = phil.parse("""\
s
  .multiple=True
{
  a=None
    .multiple=True
  b=None
    .multiple=False
}
""")
  custom = phil.parse("""\
s {
  a = 1
  a = 2
  b = 3
  b = 4
}
s {
  a = 2
  a = 3
  b = 4
  b = 5
}
s {
  a = 1
  a = 2
  b = 4
}
""")
  assert not show_diff(master.fetch(sources=[custom, custom]).as_str(), """\
s {
  a = 2
  a = 3
  b = 5
}
s {
  a = 1
  a = 2
  b = 4
}
""")
  #
  master = phil.parse("""\
a=None
  .multiple=True
b=None
  .multiple=False
""")
  custom = phil.parse("""\
a=1
a=2
b=3
b=4
""")
  fetched = master.fetch(sources=[custom, custom])
  assert not show_diff(fetched.as_str(), """\
a = 1
a = 2
b = 4
""")
  master = phil.parse("""\
s
{
  l
    .multiple=False
  {
    a=None
      .multiple=True
    b=None
      .multiple=False
  }
}
""")
  custom = phil.parse("""\
s {
  l {
    a = 1
    a = 2
    b = 3
    b = 4
  }
}
""")
  fetched = master.fetch(sources=[custom, custom])
  assert not show_diff(fetched.as_str(), """\
s {
  l {
    a = 1
    a = 2
    b = 4
  }
}
""")
  #
  master = phil.parse("""\
a=None
  .type=str
  .multiple=True
""")
  custom = phil.parse("""\
a=x
a="x"
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
a = "x"
""")
  custom = phil.parse("""\
a="x"
a=x
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
a = x
""")
  custom = phil.parse("""\
a="x "
a=x
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
a = "x "
a = x
""")
  #
  master = phil.parse("""\
a=None
  .type=str
  .multiple=True
a=x
a=y
""")
  custom = phil.parse("""\
a=z
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
a = x
a = y
a = z
""")
  custom = phil.parse("""\
a=y
a=y
a=z
a=x
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
a = y
a = z
a = x
""")
  master = phil.parse("""\
a=None
  .type=str
  .multiple=True
a=d
a=x
a=y
""")
  custom = phil.parse("""\
a=z
""")
  assert not show_diff(master.fetch(source=custom).as_str(), """\
a = d
a = x
a = y
a = z
""")
  assert not show_diff(master.fetch(sources=[custom, master]).as_str(), """\
a = z
a = d
a = x
a = y
""")
  master = phil.parse("""\
a
  .multiple=True
{
  x=None
  y=None
}
a {
  x=1
}
a {
  y=2
}
""")
  custom = phil.parse("")
  assert not show_diff(master.fetch(sources=[custom]*3).as_str(), """\
a {
  x = 1
  y = None
}
a {
  x = None
  y = 2
}
""")
  custom = phil.parse("a.x=3\na.y=4")
  assert not show_diff(master.fetch(sources=[custom]*3).as_str(), """\
a {
  x = 1
  y = None
}
a {
  x = None
  y = 2
}
a {
  x = 3
  y = None
}
a {
  x = None
  y = 4
}
""")
  custom = phil.parse("a {\nx=3\ny=4\n}")
  assert not show_diff(master.fetch(sources=[custom]*3).as_str(), """\
a {
  x = 1
  y = None
}
a {
  x = None
  y = 2
}
a {
  x = 3
  y = 4
}
""")
  assert not show_diff(master.fetch(
    sources=[custom, master, master, custom]).as_str(), """\
a {
  x = 1
  y = None
}
a {
  x = None
  y = 2
}
a {
  x = 3
  y = 4
}
""")
  assert not show_diff(master.fetch(
    sources=[custom, master, master, custom, master]).as_str(), """\
a {
  x = 3
  y = 4
}
a {
  x = 1
  y = None
}
a {
  x = None
  y = 2
}
""")
  #
  master = phil.parse(input_string="""\
s {
  t
    .multiple=true
  {
    a=None
    b=None
  }
}
""")
  source = phil.parse(input_string="""\
s.t.a=x
s.t.b=y
""")
  assert not show_diff(master.fetch(source=source).as_str(), """\
s {
  t {
    a = x
    b = None
  }
  t {
    a = None
    b = y
  }
}
""")
  master = phil.parse(input_string="""\
s.t
  .multiple=true
{
  a=None
  b=None
}
""")
  source = phil.parse(input_string="""\
s.t.a=x
s.t.b=y
""")
  assert not show_diff(master.fetch(source=source).as_str(), """\
s.t {
  a = x
  b = None
}
s.t {
  a = None
  b = y
}
""")
  master = phil.parse(input_string="""\
s
  .multiple=True
{
  t.u
    .multiple=True
  {
    a = None
  }
}
""")
  for input_string in ["", "s {}", "s { t { } }", "s.t { u { a = None } }"]:
    source = phil.parse(input_string=input_string)
    assert not show_diff(master.fetch(source=source).as_str(), """\
s {
  t.u {
    a = None
  }
}
""")
  for input_string in ["""\
s.t { u { a = 3 } }
""", """\
s.t { u { a = 3 } }
s.t { u { a = None } }
""", """\
s.t { u { a = None } }
s.t { u { a = 3 } }
s.t { u { a = None } }
"""]:
    source = phil.parse(input_string=input_string)
    assert not show_diff(master.fetch(source=source).as_str(), """\
s {
  t.u {
    a = 3
  }
}
""")
  master = phil.parse(input_string="""\
s {
  p = None
  q = None
  t
    .multiple=True
  {
    a = None
    b = None
    c.x
      .multiple=True
    {
      f = None
    }
  }
}
q
  .multiple=True
{
  r
    .multiple=True
  {
    s = None
      .multiple=True
    p = None
  }
}
""")
  source = phil.parse(input_string="""\
s.q = 1
s.p = 2
s.t.c { x { f = 3 } }
q.r.s = 4
s.t.a = 5
q.r.s = 6
s.t.b = 7
s.p = 8
s.t {
  a = 9
  b = 10
  c.x.f = 11
  c.x.f = 12
}
q { r { s=13; s=14; p=15; s=16; p=17 } r.s = 18; r.s = 19 }
s.t.a = 20
q.r.s = 21
""")
  assert not show_diff(master.fetch(source=source).as_str(), """\
s {
  p = 8
  q = 1
  t {
    a = None
    b = None
    c.x {
      f = 3
    }
  }
  t {
    a = 5
    b = None
    c.x {
      f = None
    }
  }
  t {
    a = None
    b = 7
    c.x {
      f = None
    }
  }
  t {
    a = 9
    b = 10
    c.x {
      f = 11
    }
    c.x {
      f = 12
    }
  }
  t {
    a = 20
    b = None
    c.x {
      f = None
    }
  }
}
q {
  r {
    s = 4
    p = None
  }
}
q {
  r {
    s = 6
    p = None
  }
}
q {
  r {
    s = 13
    s = 14
    s = 16
    p = 17
  }
  r {
    s = 18
    p = None
  }
  r {
    s = 19
    p = None
  }
}
q {
  r {
    s = 21
    p = None
  }
}
""")
  #
  master = phil.parse("""\
s
{
  a=None
    .multiple=True
    .optional=True
}
""")
  prev = phil.parse(master.as_str())
  custom = phil.parse("""\
s {
  a=x
}
""")
  eff_phil = master.fetch(sources=[prev, custom])
  assert not show_diff(eff_phil.as_str(), """\
s {
  a = x
}
""")
  assert master.extract().s.a == []
  assert custom.extract().s.a == ["x"]
  eff = eff_phil.extract()
  assert eff.s.a == [["x"]]
  assert not show_diff(master.format(eff).as_str(), """\
s {
  a = x
}
""")
  dft = master.clone(eff)
  assert not show_diff(master.format(dft).as_str(), """\
s {
  a = x
}
""")
  #
  master = phil.parse(input_string="""\
a = None
""")
  source = phil.parse(input_string="""\
a = 1
a = 2
""")
  working, unused = master.fetch(source=source, track_unused_definitions=True)
  assert not show_diff(working.as_str(), """\
a = 2
""")
  assert len(unused) == 0
  master = phil.parse(input_string="""\
a = 1
s {
  b = None
    .multiple = True
  b = 2
  t
    .multiple = True
  {
    c = None
    x = None
  }
}
""")
  custom = phil.parse(input_string="""\
a = 10
s.b = 20
c = 30
d = 40
s {
  t {
    c = 300
    x = 0
  }
  t.c = $c
  e = 5
}
s.t.y = 1
""")
  working, unused = master.fetch(source=custom, track_unused_definitions=True)
  assert not show_diff(working.as_str(), """\
a = 10
s {
  b = 2
  b = 20
  t {
    c = 300
    x = 0
  }
  t {
    c = 30
    x = None
  }
}
""")
  assert not show_diff(
    "\n".join([str(obj_loc) for obj_loc in unused]), """\
d (input line 4)
s.e (input line 11)
s.t.y (input line 13)""")
  #
  master = phil.parse(input_string="""\
a = *x y
  .optional = True
  .type = choice
b = *p *q
  .optional = True
  .type = choice(multi=True)
""")
  custom = phil.parse(input_string="""\
a = none
b = NoNe
""")
  working = master.fetch(source=custom)
  assert not show_diff(working.as_str(), """\
a = x y
b = p q
""")
  #
  master = phil.parse(input_string="""\
a = None
a = None
""")
  try: master.fetch()
  except RuntimeError, e:
    assert not show_diff(str(e), """\
Duplicate definitions in master (first not marked with .multiple=True):
  a (input line 1)
  a (input line 2)""")
  else: raise Exception_expected

def exercise_fetch_diff():
  source = phil.parse(input_string="""\
a = "a" "b"
""")
  master = phil.parse(input_string="""\
a = "a b"
  .type = str
""")
  assert master.objects[0].fetch_diff(source=source.objects[0]) is None
  master = phil.parse(input_string="""\
a = None
  .type = str
""")
  assert not show_diff(
    master.objects[0].fetch(source=source.objects[0], diff=True).as_str(),
    """\
a = "a" "b"
""")
  #
  master = phil.parse(input_string="""\
v = 33/111
  .type = float
""")
  f = master.fetch()
  assert not show_diff(f.objects[0].as_str(), """\
v = 33/111
""")
  assert master.objects[0].fetch_diff(source=f.objects[0]) is None
  f = master.format(master.extract(master.fetch()))
  assert not show_diff(f.objects[0].as_str(), """\
v = 0.2972972973
""")
  assert master.objects[0].fetch_diff(source=f.objects[0]) is None
  source = phil.parse(input_string="""\
v = 0.2972972972
""")
  assert not show_diff(
    master.objects[0].fetch_diff(source=source.objects[0]).as_str(),
    """\
v = 0.2972972972
""")
  #
  master = phil.parse(input_string="""\
s {
  a = None
}
""")
  source = phil.parse(input_string="s.a = x")
  d = master.fetch(source=source, diff=True)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
""")
  for input_string in ["", "s.a = None"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s {
  a = None
  b = y
}
""")
  source = phil.parse(input_string="s.a = x")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
""")
  source = phil.parse(input_string="s.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  b = z
}
""")
  source = phil.parse(input_string="s.a = x; s.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
  b = z
}
""")
  for input_string in ["", "s.a = None", "s.b = y"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s {
  a = None
  t {
    b = y
  }
}
""")
  source = phil.parse(input_string="s.a = x")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
""")
  source = phil.parse(input_string="s.t.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  t {
    b = z
  }
}
""")
  source = phil.parse(input_string="s.a = x; s.t.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
  t {
    b = z
  }
}
""")
  for input_string in ["", "s.a = None", "s.t.b = y"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s {
  a = None
  t.u {
    b = y
  }
}
""")
  source = phil.parse(input_string="s.a = x")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
""")
  source = phil.parse(input_string="s.t.u.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  t.u {
    b = z
  }
}
""")
  source = phil.parse(input_string="s.a = x; s.t.u.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
  t.u {
    b = z
  }
}
""")
  for input_string in ["", "s.a = None", "s.t.u.b = y"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
a = None
  .type = str
  .multiple = True
""")
  source = phil.parse(input_string="a = x")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
a = x
""")
  source = phil.parse(input_string="a = x; a = None; a = y")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
a = x
a = y
""")
  for input_string in ["", "a = None"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s
  .multiple = True
{
  a = None
}
""")
  source = phil.parse(input_string="s.a = x")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
""")
  for input_string in ["", "s.a = None"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s
  .multiple = True
{
  a = None
  b = y
}
""")
  source = phil.parse(input_string="s.a = x")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
""")
  source = phil.parse(input_string="s.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  b = z
}
""")
  source = phil.parse(input_string="s.a = x; s.b = z")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
s {
  b = z
}
""")
  for input_string in ["", "s.a = None", "s.b = y"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s
  .multiple = True
{
  a = None
    .type = str
    .multiple = True
  b = y
}
""")
  source = phil.parse(input_string="s.a = x")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
""")
  source = phil.parse(input_string="s.a = x; s { a = None; b = z } s.a = p")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  a = x
}
s {
  b = z
}
s {
  a = p
}
""")
  for input_string in ["", "s.a = None", "s.b = y", "s { a = None; b = y }"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s
  .multiple = True
{
  t
    .multiple = True
  {
    a = x
  }
  b = None
}
""")
  source = phil.parse(input_string="s.t.a = None")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  t {
    a = None
  }
}
""")
  source = phil.parse(input_string="""\
s.t.a = y
s { t.a = x; b = z }
s.t.a = p
s { t.a = None; b = None }
s { t.a = r; b = q }
""")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  t {
    a = y
  }
}
s {
  b = z
}
s {
  t {
    a = p
  }
}
s {
  t {
    a = None
  }
}
s {
  t {
    a = r
  }
  b = q
}
""")
  for input_string in [
        "", "s.t.a = x", "s.b = None", "s { t.a = x; b = None }"]:
    source = phil.parse(input_string=input_string)
    d = master.fetch_diff(source=source)
    assert d.is_empty()
  #
  master = phil.parse(input_string="""\
s
  .multiple = True
{
  d = None
  t.u
    .multiple = True
  {
    a = x
      .multiple = True
    b = None
    c = v *w
      .type=choice
  }
  e = None
    .multiple = True
  r {
    i = 0
      .multiple = True
    j = 1
  }
}
""")
  source = phil.parse(input_string="")
  d = master.fetch_diff(source=source)
  assert d.is_empty()
  source = phil.parse(input_string="""\
s.t.u.c=v
s {
  t.u {
    a = z
    b = y
    a = x
    c = w
    a = p
  }
  r.i = 0
  t.u {
    a = x
    b = None
    c = w
  }
  r {
    i = 1
    j = 1
  }
}
s {
  d = 2
  e = 1
  t.u {
    c = v
  }
  d = 3
  e = None
  r {
    i = 1
    j = 2
  }
  d = 4
  e = 5
}
""")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  t.u {
    c = *v w
  }
}
s {
  t.u {
    a = z
    a = p
    b = y
  }
  r {
    i = 1
  }
}
s {
  d = 4
  t.u {
    c = *v w
  }
  e = 1
  e = 5
  r {
    i = 1
    j = 2
  }
}
""")
  #
  master = phil.parse(input_string="""\
a = None
  .type = str
""")
  os.environ["_X_Y_Z_"] = "xyz"
  source = phil.parse(input_string="a = $_X_Y_Z_")
  f = master.fetch(source=source)
  assert not show_diff(f.as_str(), """\
a = "xyz"
""")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
a = "$_X_Y_Z_"
""")
  #
  master = phil.parse(input_string="""\
s {
  p = None
    .multiple = True
  p = q1
  p = q2
  t {
    u
      .multiple = True
    {
      a = None
      b = None
        .multiple = True
      b = c1
      b = c2
    }
    u {
      a = x1
      b = y1
    }
    u {
      a = x2
      b = y2
    }
  }
}
""")
  d = master.fetch_diff()
  assert d.is_empty()
  source = phil.parse(input_string="s.t.u.b=c2")
  d = master.fetch_diff(source=source)
  assert d.is_empty()
  source = phil.parse(input_string="s.t.u { b=c3; b=c4; b=c3 }")
  d = master.fetch_diff(source=source)
  assert not show_diff(d.as_str(), """\
s {
  t {
    u {
      b = c4
      b = c3
    }
  }
}
""")

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
  f=a b c
    .type=choice(multi =True)
    .optional=0
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
  int_true = True
    .type=int
  int_false = False
    .type=int
  float_true = True
    .type=float
  float_false = False
    .type=float
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
    assert str(e) == 'Multiple choices for group.e;' \
      ' only one choice can be selected (input line 13)'
  else: raise Exception_expected
  try: parameters.get(path="group.e",
    with_substitution=False).objects[2].extract()
  except RuntimeError, e:
    assert str(e) == 'Unspecified choice for group.e:' \
      ' exactly one choice must be selected (input line 15)'
  else: raise Exception_expected
  try: parameters.get(path="group.f",
    with_substitution=False).objects[3].extract()
  except RuntimeError, e:
    assert str(e) == 'Unspecified choice for group.f:' \
      ' at least one choice must be selected (input line 24)'
  else: raise Exception_expected
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
      is parameters.get(path="group.f").objects[2].type
  assert parameters.get(path="group.g",
    with_substitution=False).objects[0].extract() == "/var/tmp/foo"
  assert parameters.get(path="group.h",
    with_substitution=False).objects[0].extract() == "var.tmp.foo"
  assert parameters.get(path="group.i",
    with_substitution=False).objects[0].extract() == 2
  try: parameters.get(path="group.i",
    with_substitution=False).objects[1].extract()
  except RuntimeError, e:
    assert str(e) == 'Error interpreting group.i="1/2" as' \
      ' an integer expression (input line 33)'
  else: raise Exception_expected
  assert parameters.get(path="group.j",
    with_substitution=False).objects[0].extract() == 0.5
  try: parameters.get(path="group.j",
    with_substitution=False).objects[1].extract()
  except RuntimeError, e:
    assert str(e) == """Error interpreting group.j="'a'" as""" \
      " a floating-point expression (input line 37)"
  else: raise Exception_expected
  try: parameters.get(path="group.j",
    with_substitution=False).objects[2].extract()
  except RuntimeError, e:
    assert str(e) == 'Error interpreting group.j="a" as a numeric expression:'\
                   + " NameError: name 'a' is not defined (input line 39)"
  else: raise Exception_expected
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
  try: parameters.get(path="group.int_true",
    with_substitution=False).objects[0].extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      'Error interpreting group.int_true="True" as a numeric expression'
      ' (input line 53)')
  else: raise Exception_expected
  try: parameters.get(path="group.int_false",
    with_substitution=False).objects[0].extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      'Error interpreting group.int_false="False" as a numeric expression'
      ' (input line 55)')
  else: raise Exception_expected
  try: parameters.get(path="group.float_true",
    with_substitution=False).objects[0].extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      'Error interpreting group.float_true="True" as a numeric expression'
      ' (input line 57)')
  else: raise Exception_expected
  try: parameters.get(path="group.float_false",
    with_substitution=False).objects[0].extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      'Error interpreting group.float_false="False" as a numeric expression'
      ' (input line 59)')
  else: raise Exception_expected
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
  #
  parameters = phil.parse(input_string="""\
s {
  a=None
   .optional=True
   .type=int
   .multiple=True
}
""")
  extracted = parameters.extract()
  assert extracted.s.a == []
  #
  parameters = phil.parse(input_string="""\
s {
  t {
    u {
    }
  }
}
t {
  b.c=None
  c.d.e=f
    .type = str
}
""")
  extracted = parameters.extract()
  assert extracted.__phil_path__() == ""
  assert extracted.__phil_path__(object_name="o1") == "o1"
  assert extracted.s.__phil_path__() == "s"
  assert extracted.s.__phil_path__(object_name="o2") == "s.o2"
  assert extracted.s.t.__phil_path__() == "s.t"
  assert extracted.s.t.__phil_path__(object_name="o3") == "s.t.o3"
  assert extracted.s.t.u.__phil_path__() == "s.t.u"
  assert extracted.t.__phil_path__() == "t"
  assert extracted.t.__phil_path__(object_name="o4") == "t.o4"
  assert extracted.t.b.__phil_path__() == "t.b"
  assert extracted.t.c.__phil_path__() == "t.c"
  assert extracted.t.c.d.__phil_path__() == "t.c.d"
  assert extracted.t.c.d.__phil_path_and_value__(object_name="e") == (
    "t.c.d.e", "f")
  #
  try: extracted.z = 12
  except AttributeError, e:
    assert not show_diff(str(e), """Assignment to non-existing attribute "z"
  Please correct the attribute name, or to create
  a new attribute use: obj.__inject__(name, value)""")
  else: raise Exception_expected
  try: extracted.t.c.z = 13
  except AttributeError, e:
    assert not show_diff(str(e),"""Assignment to non-existing attribute "t.c.z"
  Please correct the attribute name, or to create
  a new attribute use: obj.__inject__(name, value)""")
  else: raise Exception_expected
  extracted.__inject__("z", 14)
  assert extracted.z == 14
  del extracted.z
  extracted.__inject__("z", 15)
  assert extracted.z == 15
  try: extracted.__inject__("z", 16)
  except AttributeError, e:
    assert str(e) == 'Attribute "z" exists already.'
  else: raise Exception_expected
  extracted.t.c.__inject__("z", 17)
  assert extracted.t.c.z == 17
  del extracted.t.c.z
  extracted.t.c.__inject__("z", 18)
  assert extracted.t.c.z == 18
  try: extracted.t.c.__inject__("z", 19)
  except AttributeError, e:
    assert str(e) == 'Attribute "t.c.z" exists already.'
  else: raise Exception_expected
  #
  parameters = phil.parse(input_string="""\
s.t.x {
  a=13
    .type=int
}
s.t.y {
  b=45
    .type=int
}
""")
  extracted = parameters.extract()
  assert extracted.s.t.x.a == 13
  assert extracted.s.t.y.b == 45
  parameters = phil.parse(input_string="""\
s.t.x {
  a=13
    .type=int
    .multiple=True
}
s.t.x {
  a=None
    .type=int
    .multiple=True
}
s.t.x {
  a=45
    .type=int
    .multiple=True
}
""")
  extracted = parameters.extract()
  assert extracted.s.t.x.a == [13, 45]
  parameters = phil.parse(input_string="""\
s.t.x {
  a=None
    .type=int
    .multiple=True
}
s.t.x {
  a=45
    .type=int
    .multiple=True
}
s.t.x {
  a=13
    .type=int
    .multiple=True
}
""")
  extracted = parameters.extract()
  assert extracted.s.t.x.a == [45, 13]
  #
  parameters = phil.parse(input_string="""\
s.t {
}
s.u
  .multiple=True
  .optional=True
{
}
""")
  for object in parameters.get(path="s").objects:
    assert object.optional is None
    assert object.multiple is None
  for object in parameters.get(path="s.u").objects:
    assert object.optional
    assert object.multiple
  extracted = parameters.extract()
  assert len(extracted.s.u) == 1
  #
  parameters = phil.parse(input_string="""\
a=008 # eval fails since this is an invalid octal literal
  .type=int
b=008.34 # just to be sure (eval works for this)
  .type=float
""")
  extracted = parameters.extract()
  assert extracted.a == 8
  assert extracted.b == 8.34
  #
  master = phil.parse(input_string="""\
a=None
  .type=qstr
b=Auto
  .type=qstr
c= 1  2
  .type=qstr
ds="1"'2'
  .type=str
dq="1"'2'
  .type=qstr
es="1\\""'2\\"\\''
  .type=str
eq="1\\""'2\\"\\''
  .type=qstr
""")
  work = master
  for i_cycle in xrange(2):
    ex = work.extract()
    assert ex.a is None
    assert ex.b is Auto
    assert ex.c == "1 2"
    assert ex.ds == "1 2"
    assert ex.dq == """ "1" '2' """[1:-1]
    assert ex.es == """ 1" 2\\"' """[1:-1]
    assert ex.eq == """ "1\\"" '2\\\\"\\'' """[1:-1]
    work = master.format(ex)
    assert not show_diff(work.as_str(), """\
a = None
b = Auto
c = 1 2
ds = "1 2"
dq = "1" '2'
es = "1\\" 2\\\\\\"'"
eq = "1\\"" '2\\\\"\\''
""")

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
  assert not show_diff(out.getvalue(), """\
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
""")
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
!a=2
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
  assert not show_diff(out.getvalue(), """\
a = 1
a = 2
""")
  custom = phil.parse(input_string="""\
a=1
!a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == [1]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
a = 1
""")
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
a = None
""")
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
  assert not show_diff(out.getvalue(), """\
a = 1
a = 2
""")
  custom = phil.parse(input_string="""\
a=1
!a=2
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == [1]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
a = 1
""")
  custom = phil.parse(input_string="""\
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
a = None
""")
  custom = phil.parse(input_string="""\
a=None
""")
  extracted = master.fetch(custom).extract()
  assert extracted.a == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
a = None
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = None
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = 2
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = 3
}
""")
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
  assert len(extracted.s) == 0
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
s {
  a = None
}
""")
  custom = phil.parse(input_string="""\
s {
  a=None
}
""")
  extracted = master.fetch(custom).extract()
  assert len(extracted.s) == 0
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
s {
  a = None
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = 1
}
s {
  a = 2
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = 1
}
s {
  a = 3
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = None
}
""")
  custom = phil.parse(input_string="""\
s {
  a=None
}
""")
  extracted = master.fetch(custom).extract()
  assert extracted.s == []
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
s {
  a = None
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = 1
}
s {
  a = 2
}
""")
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
  assert not show_diff(out.getvalue(), """\
s {
  a = 1
}
s {
  a = 3
}
""")
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
  assert extracted.b[1].a == []
  assert extracted.b[1].b == [5,6]
  assert extracted.b[1].c is None
  assert extracted.c == []
  assert extracted.d[0].a is None
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
a = 1
a = 2
b {
  a = 3
  a = 4
  b = None
  c = 20
}
b {
  a = None
  b = 5
  b = 6
  c = None
}
c {
  a = None
}
d {
  a = None
}
""")
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
  #
  master = phil.parse(input_string="""\
a = x y
  .optional = True
  .type = choice
b = *p q
  .optional = True
  .type = choice
c = *r s
  .optional = True
  .type = choice
d = g h
  .optional = True
  .type = choice(multi=True)
e = *j *k
  .optional = True
  .type = choice(multi=True)
f = *l *m
  .optional = True
  .type = choice(multi=True)
""")
  extracted = master.extract()
  assert extracted.a is None
  assert extracted.b == "p"
  assert extracted.c == "r"
  assert extracted.d == []
  assert extracted.e == ["j", "k"]
  assert extracted.f == ["l", "m"]
  extracted.c = None
  del extracted.f[:]
  out = StringIO()
  master.format(extracted).show(out=out)
  assert not show_diff(out.getvalue(), """\
a = x y
b = *p q
c = r s
d = g h
e = *j *k
f = l m
""")
  #
  for input_string in ["""
a = None
  .type = str
  .multiple = True
""", """\
a = None
  .optional = True
  .type = str
  .multiple = True
""", """\
a = None
  .optional = False
  .type = str
  .multiple = True
"""]:
    master = phil.parse(input_string=input_string)
    source = phil.parse(input_string="")
    w = master.fetch(source=source)
    e = w.extract()
    opt = master.objects[0].optional
    if (opt is None or opt):
      assert len(e.a) == 0
    else:
      assert len(e.a) == 1
      e.a[0] is None
    source = phil.parse(input_string="a=None")
    w = master.fetch(source=source)
    e = w.extract()
    if (opt is None or opt):
      assert len(e.a) == 0
    else:
      assert len(e.a) == 1
      e.a[0] is None
    source = phil.parse(input_string="a=x")
    w = master.fetch(source=source)
    e = w.extract()
    if (opt is None or opt):
      assert e.a == ["x"]
    else:
      assert len(e.a) == 2
      e.a[0] is None
      e.a[1] == "x"
  #
  for input_string in ["""\
s
  .multiple = True
{
  a = None
    .type = str
}
""", """
s
  .optional = True
  .multiple = True
{
  a = None
    .type = str
}
""", """
s
  .optional = False
  .multiple = True
{
  a = None
    .type = str
}
"""]:
    master = phil.parse(input_string=input_string)
    source = phil.parse(input_string="")
    w = master.fetch(source=source)
    e = w.extract()
    opt = master.objects[0].optional
    if (opt is None or opt):
      assert len(e.s) == 0
    else:
      assert len(e.s) == 1
      assert e.s[0].a is None
    source = phil.parse(input_string="s {}")
    w = master.fetch(source=source)
    e = w.extract()
    if (opt is None or opt):
      assert len(e.s) == 0
    else:
      assert len(e.s) == 1
      assert e.s[0].a is None
    source = phil.parse(input_string="s { a=x }")
    w = master.fetch(source=source)
    e = w.extract()
    if (opt is None or opt):
      assert len(e.s) == 1
      assert e.s[0].a == "x"
    else:
      assert len(e.s) == 2
      assert e.s[0].a is None
      assert e.s[1].a == "x"
  #
  master = phil.parse(input_string="""\
s
  .optional = True
  .multiple = True
{
  a = None
    .type = str
  b = 0
    .type = float
  c = *x *y
    .optional = True
    .type = choice(multi=True)
}
""")
  for input_string in ["", "s {}"]:
    source = phil.parse(input_string=input_string)
    w = master.fetch(source=source)
    e = w.extract()
    assert len(e.s) == 0
    f = master.format(python_object=e)
    assert f.objects[0].is_template == 1
    assert not show_diff(f.as_str(), """\
s {
  a = None
  b = 0
  c = *x *y
}
""")
    assert not show_diff(f.as_str(attributes_level=2), """\
s
  .optional = True
  .multiple = True
{
  a = None
    .type = str
  b = 0
    .type = float
  c = *x *y
    .optional = True
    .type = choice(multi=True)
}
""")
  source = phil.parse(input_string="""\
s.b=1
""")
  w = master.fetch(source=source)
  e = w.extract()
  assert len(e.s) == 1
  assert e.s[0].b == 1
  f = master.format(python_object=e)
  assert f.objects[0].is_template == -1
  assert f.objects[1].is_template == 0
  assert not show_diff(f.as_str(), """\
s {
  a = None
  b = 1
  c = *x *y
}
""")
  assert not show_diff(f.as_str(attributes_level=2), """\
s
  .optional = True
  .multiple = True
{
  a = None
    .type = str
  b = 0
    .type = float
  c = *x *y
    .optional = True
    .type = choice(multi=True)
}
s
  .optional = True
  .multiple = True
{
  a = None
    .type = str
  b = 1
    .type = float
  c = *x *y
    .optional = True
    .type = choice(multi=True)
}
""")
  source = phil.parse(input_string="""\
s.c=*x
s.c=*y
""")
  w = master.fetch(source=source)
  e = w.extract()
  assert len(e.s) == 2
  assert e.s[0].c == ["x"]
  assert e.s[1].c == ["y"]
  f = master.format(python_object=e)
  assert f.objects[0].is_template == -1
  assert f.objects[1].is_template == 0
  assert f.objects[2].is_template == 0
  assert not show_diff(f.as_str(), """\
s {
  a = None
  b = 0
  c = *x y
}
s {
  a = None
  b = 0
  c = x *y
}
""")
  #
  master = phil.parse(input_string="""\
s
  .multiple = True
{
  a = x
    .type = str
  b = 2
    .type = int
  c = x *y
    .optional = True
    .type = choice(multi=True)
  d = None
    .type = float
}
s {
}
""")
  source = phil.parse(input_string="""\
""")
  w = master.fetch(source=source)
  e = w.extract()
  assert len(e.s) == 0
  f = master.format(python_object=e)
  assert not show_diff(f.as_str(), """\
s {
  a = x
  b = 2
  c = x *y
  d = None
}
""")
  assert not show_diff(f.as_str(attributes_level=2), """\
s
  .multiple = True
{
  a = x
    .type = str
  b = 2
    .type = int
  c = x *y
    .optional = True
    .type = choice(multi=True)
  d = None
    .type = float
}
""")
  #
  master = phil.parse(input_string="""\
a = None
  .optional = True
  .multiple = True
""")
  source = phil.parse(input_string="""\
""")
  w = master.fetch(source=source)
  e = w.extract()
  assert len(e.a) == 0
  f = master.format(python_object=e)
  assert not show_diff(f.as_str(), """\
a = None
""")
  assert not show_diff(f.as_str(attributes_level=2), """\
a = None
  .optional = True
  .multiple = True
""")
  source = phil.parse(input_string="""\
a=x
""")
  w = master.fetch(source=source)
  e = w.extract()
  assert e.a == [["x"]]
  f = master.format(python_object=e)
  assert not show_diff(f.as_str(), """\
a = x
""")
  assert not show_diff(f.as_str(attributes_level=2), """\
a = x
  .optional = True
  .multiple = True
""")
  source = phil.parse(input_string="""\
a=x
a=y
""")
  w = master.fetch(source=source)
  e = w.extract()
  assert len(e.a) == 2
  assert e.a == [["x"], ["y"]]
  f = master.format(python_object=e)
  assert not show_diff(f.as_str(), """\
a = x
a = y
""")
  assert not show_diff(f.as_str(attributes_level=2), """\
a = x
  .optional = True
  .multiple = True
a = y
  .optional = True
  .multiple = True
""")

class foo1_converters(object):

  phil_type = "foo1"

  def __str__(self): return "foo1"

  def __init__(self, bar=None):
    if (bar is None):
      raise RuntimeError("foo1 problem")

class foo2_converters(object):

  def __str__(self): return "foo2"

  def __init__(self, bar=None):
    if (bar is not None):
      raise RuntimeError("foo2 problem")

foo_converter_registry = phil.extended_converter_registry(
  additional_converters=[foo1_converters, foo2_converters])

def exercise_choice():
  master = phil.parse(input_string="""\
x=*a b
  .type=choice
y=*a b
  .type=choice
  .optional=True
""")
  py_params = master.extract()
  py_params.x = "c"
  try: master.format(py_params)
  except RuntimeError, e:
    assert str(e) == "Invalid choice: x=c"
  else: raise Exception_expected
  py_params.x = "b"
  py_params.y = "c"
  try: master.format(py_params)
  except RuntimeError, e:
    assert str(e) == "Invalid choice: y=c"
  else: raise Exception_expected
  master = phil.parse(input_string="""\
x=*a a
  .type=choice
""")
  py_params = master.extract()
  try: master.format(py_params)
  except RuntimeError, e:
    assert str(e) \
        == "Improper master choice definition: x = *a a (input line 1)"
  else: raise Exception_expected
  #
  master = phil.parse(input_string="""\
x=*a b
  .type=choice(multi=True)
  .optional=False
y=a b
  .type=choice(multi=True)
  .optional=True
""")
  py_params = master.extract()
  py_params.x = ["a", "c", "d"]
  try: master.format(py_params)
  except RuntimeError, e:
    assert str(e) == "Invalid choice(multi=True): x=['c', 'd']"
  else: raise Exception_expected
  py_params.x = []
  try: master.format(py_params)
  except RuntimeError, e:
    assert str(e) == "Empty list for mandatory choice(multi=True): x"
  else: raise Exception_expected
  master = phil.parse(input_string="""\
x=*a b
  .type=choice(multi=True)
""")
  py_params = master.extract()
  py_params.x = []
  assert not show_diff(master.format(py_params).as_str(), """\
x = a b
""")
  #
  master = phil.parse(input_string="""\
x=*a b
  .type=choice
""")
  py_params = master.extract()
  py_params.x = None
  assert not show_diff(master.format(py_params).as_str(), """\
x = a b
""")
  master = phil.parse(input_string="""\
x=*a b
  .type=choice
  .optional=False
""")
  py_params = master.extract()
  py_params.x = None
  try: master.format(py_params)
  except RuntimeError, e:
    assert str(e) == "Invalid choice: x=None"
  else: raise Exception_expected
  master_phil = phil.parse("""\
multiwordchoice = 'a b' 'b c'
  .type = choice
""")
  other = phil.parse("""
multiwordchoice = 'a b'
""")
  s = StringIO()
  master_phil.fetch_diff(source=other).show(out=s)
  diff_phil = phil.parse(input_string=s.getvalue())
  working_phil = master_phil.fetch(source=diff_phil)
  assert not show_diff(working_phil.as_str(),"multiwordchoice = '*a b' 'b c'\n")
  master_phil = phil.parse("""\
x = a b
  .type = choice
y = *a B
  .type = choice
""")
  other = phil.parse("""\
x = A
y = b
""")
  params = master_phil.extract()
  assert params.x is None
  assert params.y == 'a'
  params = master_phil.fetch(source=other).extract()
  assert params.x == 'a'
  assert params.y == 'B'

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
      " AttributeError: 'foo2_converters' object has no attribute 'from_words'"
  else: raise Exception_expected
  try: params.get(path="a").objects[0].format(python_object=0)
  except RuntimeError, e:
    assert str(e) == \
      '.type=foo2 does not have an as_words method (input line 1):' \
      " AttributeError: 'foo2_converters' object has no attribute 'as_words'"
  else: raise Exception_expected
  try: phil.parse(
    input_string="""\
a=None
  .type=foo1
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e) == \
      'Error constructing definition type "foo1": RuntimeError:' \
      ' foo1 problem (input line 2)'
  else: raise Exception_expected
  try: phil.parse(
    input_string="""\
a=None
  .type=foo2(bar=1)
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    e = str(e)
    assert e.startswith(
      'Error constructing definition type "foo2(bar=1)": RuntimeError: ')
    assert e.endswith(' foo2 problem (input line 2)')
  else: raise Exception_expected
  try: phil.parse(
    input_string="""\
a=None
  .type=foo2(foo=1)
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e) == \
      'Error constructing definition type "foo2(foo=1)": TypeError:' \
      " __init__() got an unexpected keyword argument 'foo' (input line 2)"
  else: raise Exception_expected
  try: phil.parse("""\
a=None
  .type=foo2(foo=1
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e).startswith('Error evaluating definition type "foo2(foo=1": ')
    assert str(e).endswith(' (input line 2)')
  else: raise Exception_expected
  try: phil.parse("""\
a=None
  .type=foo2(foo=1=2)
""",
    converter_registry=foo_converter_registry)
  except RuntimeError, e:
    assert str(e) == 'Error constructing definition type "foo2(foo=1=2)":' \
      ' SyntaxError: invalid syntax (line 1) (input line 2)'
  else: raise Exception_expected

def exercise_auto():
  for nao in [None, Auto]:
    na = str(nao)
    master = phil.parse(input_string="""\
abool=%(na)s
  .type=bool
aint=%(na)s
  .type=int
afloat=%(na)s
  .type=float
astr=%(na)s
  .type=str
achoice=a b
  .type=choice
achoicemand=*a b
  .type=choice
  .optional=False
achoicemult=a b
  .type=choice(multi=True)
achoicemultmand=*a b
  .type=choice(multi=True)
  .optional=False
apath=%(na)s
  .type=path
akey=%(na)s
  .type=key
anone=%(na)s
  .type=None
astrings=%(na)s
  .type=strings
awords=%(na)s
  .type=words
""" % vars())
    for i in xrange(3):
      params = master.extract()
      assert params.abool is nao
      assert params.aint is nao
      assert params.afloat is nao
      assert params.astr is nao
      assert params.achoice is None
      assert params.achoicemand == "a"
      assert params.achoicemult == []
      assert params.achoicemultmand == ["a"]
      assert params.apath is nao
      assert params.akey is nao
      assert params.anone is nao
      assert params.astrings is nao
      assert params.awords is nao
      master = master.fetch()
    assert not show_diff(master.format(params).as_str(), """\
abool = %(na)s
aint = %(na)s
afloat = %(na)s
astr = %(na)s
achoice = a b
achoicemand = *a b
achoicemult = a b
achoicemultmand = *a b
apath = %(na)s
akey = %(na)s
anone = %(na)s
astrings = %(na)s
awords = %(na)s
""" % vars())
    custom = phil.parse(input_string="""\
abool=Auto
aint=Auto
afloat=Auto
astr=Auto
achoice=Auto
achoicemand=Auto
achoicemult=Auto
achoicemultmand=Auto
apath=Auto
akey=Auto
anone=Auto
astrings=Auto
awords=Auto
""")
    params = master.fetch(source=custom).extract()
    assert params.abool is Auto
    assert params.aint is Auto
    assert params.afloat is Auto
    assert params.astr is Auto
    assert params.achoice is Auto
    assert params.achoicemand is Auto
    assert params.achoicemult is Auto
    assert params.achoicemultmand is Auto
    assert params.apath is Auto
    assert params.akey is Auto
    assert params.anone is Auto
    assert params.astrings is Auto
    assert params.awords is Auto
    assert not show_diff(master.format(params).as_str(), """\
abool = Auto
aint = Auto
afloat = Auto
astr = Auto
achoice = Auto
achoicemand = Auto
achoicemult = Auto
achoicemultmand = Auto
apath = Auto
akey = Auto
anone = Auto
astrings = Auto
awords = Auto
""")

def exercise_int_and_float():
  master_phil = phil.parse(input_string="""\
a=None
  .type=float
b=1.0
  .type=int
c=1.5
  .type=float(value_min=1.0)
d=4.5
  .type=float(value_max=5.0)
e=6.2
  .type=float(value_min=3.2, value_max=7.6)
f = 5
  .type = int(value_min=0)
g = 3
  .type = int(value_max=4)
h = 10
  .type = int(value_min=4, value_max=12)
""")
  work_params = master_phil.extract()
  work_phil = master_phil.format(python_object=work_params)
  assert not show_diff(work_phil.as_str(attributes_level=2),"""\
a = None
  .type = float
b = 1
  .type = int
c = 1.5
  .type = float(value_min=1)
d = 4.5
  .type = float(value_max=5)
e = 6.2
  .type = float(value_min=3.2, value_max=7.6)
f = 5
  .type = int(value_min=0)
g = 3
  .type = int(value_max=4)
h = 10
  .type = int(value_min=4, value_max=12)
""")
  #
  master_phil = phil.parse(input_string="""\
a=-1.5
  .type=float(value_min=1.0)
""")
  try: master_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "a element is less than the minimum allowed value:"
      " -1.5 < 1 (input line 1)")
  else: raise Exception_expected
  #
  master_phil = phil.parse(input_string="""\
b=5.1
  .type=float(value_min=0.0, value_max=5.0)
""")
  try: master_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "b element is greater than the maximum allowed value:"
      " 5.1 > 5 (input line 1)")
  else: raise Exception_expected
  #
  master_phil = phil.parse(input_string="""\
c=7
  .type=int(value_max=6)
""")
  try: master_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "c element is greater than the maximum allowed value:"
      " 7 > 6 (input line 1)")
  else: raise Exception_expected
  #
  master_phil = phil.parse(input_string="""\
d=-6
  .type=int(value_min=3, value_max=25)
""")
  try: master_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "d element is less than the minimum allowed value:"
      " -6 < 3 (input line 1)")
  else: raise Exception_expected

def exercise_ints_and_floats():
  master_phil = phil.parse(input_string="""\
a=None
  .type=floats
b=Auto
  .type=ints
c=1
  .type=floats
d=1 1/2
  .type=floats
e=2,3/4,4/5
  .type=floats
f="1;1/2"
  .type=floats
g=None
  .type=floats
""")
  work_params = master_phil.extract()
  work_params.g = [1/4]
  work_phil = master_phil.format(python_object=work_params)
  assert not show_diff(work_phil.as_str(attributes_level=2), """\
a = None
  .type = floats
b = Auto
  .type = ints
c = 1
  .type = floats
d = 1 0.5
  .type = floats
e = 2 0.75 0.8
  .type = floats
f = 1 0.5
  .type = floats
g = 0.25
  .type = floats
""")
  #
  master_phil = phil.parse(input_string="""\
a=1 v
  .type=floats
""")
  try: master_phil.extract()
  except RuntimeError, e:
    assert str(e) == """\
Error interpreting a="v" as a numeric expression:\
 NameError: name 'v' is not defined (input line 1)"""
  else: raise Exception_expected
  #
  master_phil = phil.parse(input_string="""\
a=1 2
  .type=ints(2)
b=3
  .type=floats(size=1)
c=4 5
  .type=ints(size_min=2)
d=6
  .type=floats(size_max=2)
e=7 8 9
  .type=floats(size_max=4, size_min=2)
f=10 20
  .type=floats(value_max=20, value_min=10)
g=nOnE aUtO
  .type=floats(allow_none_elements=True, allow_auto_elements=True)
""")
  work_params = master_phil.extract()
  work_phil = master_phil.format(python_object=work_params)
  assert not show_diff(work_phil.as_str(attributes_level=2), """\
a = 1 2
  .type = ints(size=2)
b = 3
  .type = floats(size=1)
c = 4 5
  .type = ints(size_min=2)
d = 6
  .type = floats(size_max=2)
e = 7 8 9
  .type = floats(size_min=2, size_max=4)
f = 10 20
  .type = floats(value_min=10, value_max=20)
g = None Auto
  .type = floats(allow_none_elements=True, allow_auto_elements=True)
""")
  #
  user_phil = phil.parse(input_string="""\
f=9
""")
  work_phil = master_phil.fetch(source=user_phil)
  try: work_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "f element is less than the minimum allowed value:"
      " 9 < 10 (input line 1)")
  else: raise Exception_expected
  user_phil = phil.parse(input_string="""\
f=10
""")
  work_phil = master_phil.fetch(source=user_phil)
  work_params = work_phil.extract()
  work_params.f[0] = 9
  try: master_phil.format(python_object=work_params)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "f element is less than the minimum allowed value: 9 < 10")
  else: raise Exception_expected
  #
  user_phil = phil.parse(input_string="""\
f=21
""")
  work_phil = master_phil.fetch(source=user_phil)
  try: work_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "f element is greater than the maximum allowed value:"
      " 21 > 20 (input line 1)")
  else: raise Exception_expected
  user_phil = phil.parse(input_string="""\
f=20
""")
  work_phil = master_phil.fetch(source=user_phil)
  work_params = work_phil.extract()
  work_params.f[0] = 21
  try: master_phil.format(python_object=work_params)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "f element is greater than the maximum allowed value: 21 > 20")
  else: raise Exception_expected
  #
  user_phil = phil.parse(input_string="""\
e=1
""")
  work_phil = master_phil.fetch(source=user_phil)
  try: work_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Not enough values for e: 1 given, at least 2 required (input line 1)")
  else: raise Exception_expected
  work_params = master_phil.extract()
  work_params.e = [1]
  try: master_phil.format(python_object=work_params)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Not enough values for e: 1 given, at least 2 required")
  else: raise Exception_expected
  #
  user_phil = phil.parse(input_string="""\
a=1
""")
  work_phil = master_phil.fetch(source=user_phil)
  try: work_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Not enough values for a: 1 given, exactly 2 required (input line 1)")
  else: raise Exception_expected
  work_params = master_phil.extract()
  work_params.a = [1]
  try: master_phil.format(python_object=work_params)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Not enough values for a: 1 given, exactly 2 required")
  else: raise Exception_expected
  #
  user_phil = phil.parse(input_string="""\
e=1 2 3 4 5
""")
  work_phil = master_phil.fetch(source=user_phil)
  try: work_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Too many values for e: 5 given, 4 allowed at most (input line 1)")
  else: raise Exception_expected
  work_params = master_phil.extract()
  work_params.e = [1,2,3,4,5]
  try: master_phil.format(python_object=work_params)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Too many values for e: 5 given, 4 allowed at most")
  else: raise Exception_expected
  #
  user_phil = phil.parse(input_string="""\
a=1 2 3
""")
  work_phil = master_phil.fetch(source=user_phil)
  try: work_phil.extract()
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Too many values for a: 3 given, exactly 2 required (input line 1)")
  else: raise Exception_expected
  work_params = master_phil.extract()
  work_params.a = [1,2,3]
  try: master_phil.format(python_object=work_params)
  except RuntimeError, e:
    assert not show_diff(str(e),
      "Too many values for a: 3 given, exactly 2 required")
  else: raise Exception_expected
  #
  for value in ["None", "Auto"]:
    user_phil = phil.parse(input_string="""\
a=1 %s
""" % value.lower())
    work_phil = master_phil.fetch(source=user_phil)
    try: work_phil.extract()
    except RuntimeError, e:
      assert not show_diff(str(e),
        "a element cannot be %s (input line 1)" % value)
    else: raise Exception_expected
    work_params = master_phil.extract()
    work_params.a = [1, eval(value)]
    try: master_phil.format(python_object=work_params)
    except RuntimeError, e:
      assert not show_diff(str(e),
        "a element cannot be %s" % value)
    else: raise Exception_expected
  #
  master_phil = phil.parse(input_string="""\
a = None
  .type = floats
b = None
  .type = ints
""")
  user_phil = phil.parse(input_string="""\
a = ([12.34, 45.67])
b = [(3,4)]
""")
  work_phil = master_phil.fetch(user_phil)
  assert not show_diff(work_phil.extract_format().as_str(), """\
a = 12.34 45.67
b = 3 4
""")

def exercise_definition_validate_etc():
  working_phil = phil.parse(input_string="""\
a=None
  .type=ints
""")
  a = working_phil.objects[0]
  proxy = a.try_tokenize(input_string=" 3 ")
  assert proxy.error_message is None
  assert not show_diff(proxy.tokenized.as_str(), "a = 3\n")
  proxy = a.try_tokenize(input_string="'", source_info="src")
  assert not show_diff(
    proxy.error_message, 'Syntax error: missing closing quote (src, line 1)')
  assert proxy.tokenized is None
  proxy = a.try_tokenize(input_string=" 3,  2")
  assert proxy.error_message is None
  assert not show_diff(proxy.tokenized.as_str(), "a = 3, 2\n")
  proxy2 = proxy.tokenized.try_extract()
  assert proxy2.error_message is None
  assert proxy2.extracted == [3, 2]
  proxy2 = proxy.tokenized.try_extract_format()
  assert proxy2.error_message is None
  assert not show_diff(proxy2.formatted.as_str(), "a = 3 2\n")
  proxy = a.try_tokenize(input_string="'x y'", source_info="Src")
  assert proxy.error_message is None
  assert not show_diff(proxy.tokenized.as_str(), "a = 'x y'\n")
  proxy2 = proxy.tokenized.try_extract()
  assert not show_diff(proxy2.error_message, """\
Error interpreting a="x" as a numeric expression:\
 NameError: name 'x' is not defined (Src, line 1)""")
  assert proxy2.extracted == None
  proxy2f = proxy.tokenized.try_extract_format()
  assert not show_diff(proxy2f.error_message, proxy2.error_message)
  assert proxy2f.formatted is None
  proxy = a.validate(input_string="1;2,3")
  assert proxy.error_message is None
  assert proxy.extracted == [1,2,3]
  proxy = a.validate_and_format(input_string="1;2,3\n4")
  assert proxy.error_message is None
  assert not show_diff(proxy.formatted.as_str(), "a = 1 2 3 4\n")
  for v,m in [(a.validate, "extracted"), (a.validate_and_format, "formatted")]:
    proxy = v(input_string="*", source_info="si")
    assert not show_diff(proxy.error_message, """\
Error interpreting a="*" as a numeric expression:\
 SyntaxError: unexpected EOF while parsing (line 1) (si, line 1)""")
    assert getattr(proxy, m) is None
  proxy = a.validate(input_string="1;2,3")
  #
  proxy = a.validate_and_format(input_string=" 5*2 3 7 ")
  a.words = proxy.formatted.words
  assert not show_diff(working_phil.as_str(), "a = 10 3 7\n")
  #
  working_phil = phil.parse(input_string="""\
a=None
  .type=None
""")
  a = working_phil.objects[0]
  proxy = a.try_tokenize(input_string="")
  assert proxy.error_message is None
  proxy2 = proxy.tokenized.try_extract()
  assert proxy2.error_message is None
  assert proxy2.extracted is None

def exercise_command_line():
  master_string = """\
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
  flag=None
}
"""
  master_phil = phil.parse(input_string=master_string)
  itpr_bar = master_phil.command_line_argument_interpreter(
    home_scope="bar")
  itpr_neutral = master_phil.command_line_argument_interpreter()
  for itpr in [itpr_bar, itpr_neutral]:
    itpr.process(arg="foo.limit=4\nbar.max=2").as_str() == """\
foo.limit = 4
bar.max = 2
"""
  assert itpr_bar.process(arg="max=6").as_str() == "bar.max = 6\n"
  assert itpr_bar.process(arg="ax=7").as_str() == "bar.max = 7\n"
  try: assert itpr_neutral.process(arg="max=5")
  except Sorry, e:
    assert not show_diff(str(e), """\
Ambiguous parameter definition: max = 5
Best matches:
  foo.max
  bar.max""")
  else: raise Exception_expected
  assert itpr_bar.process(arg="limit=0").as_str() == "bar.sub.limit = 0\n"
  assert itpr_bar.process(arg="imit=0").as_str() == "bar.sub.limit = 0\n"
  for itpr in [itpr_bar, itpr_neutral]:
    assert itpr.process(arg="index=0").as_str() == "foo.index = 0\n"
    assert itpr.process(arg="ndex=0").as_str() == "foo.index = 0\n"
  try: itpr_bar.process(arg="xyz=")
  except Sorry, e:
    assert not show_diff(str(e), """\
Error interpreting command line argument as parameter definition:
  "xyz="
  RuntimeError: Missing value for xyz (command line argument, line 1)""")
  else: raise Exception_expected
  try: itpr_bar.process(arg="xyz=8")
  except Sorry, e:
    assert str(e) == "Unknown command line parameter definition: xyz = 8"
  else: raise Exception_expected
  try: itpr_bar.process(arg="  ")
  except Sorry, e:
    assert str(e) == 'Command line parameter definition has no effect: "  "'
  else: raise Exception_expected
  itpr = master_phil.command_line_argument_interpreter(
    argument_description="")
  try: itpr.process(arg="bar {}")
  except Sorry, e:
    assert str(e) == 'Parameter definition has no effect: "bar {}"'
  else: raise Exception_expected
  #
  print >> open("tmp0d5f6e10.phil", "w"), "foo.limit=-3"
  user_phils = itpr_bar.process(args=[
    "",
    "--flag",
    "--flag=no",
    "tmp0d5f6e10.phil",
    "max=8",
    "limit=9"])
  assert len(user_phils) == 5
  for i,expected in enumerate([
        "bar.flag = True\n",
        "bar.flag = no\n",
        "foo.limit = -3\n",
        "bar.max = 8\n",
        "bar.sub.limit = 9\n"]):
    assert not show_diff(user_phils[i].as_str(), expected)
  os.remove("tmp0d5f6e10.phil")
  assert not os.path.exists("tmp0d5f6e10.phil")
  for arg in ["tmp0d5f6e10.phil", "lmit=3"]:
    try: itpr_bar.process(args=[arg])
    except Sorry, e:
      assert not show_diff(str(e),
        'Uninterpretable command line argument: "%s"' % arg)
    else: raise Exception_expected
  print >> open("tmp0d5f6e10.phil", "w"), "foo$limit=0"
  try: itpr_bar.process(args=["tmp0d5f6e10.phil"])
  except RuntimeError, e:
    assert not show_diff(str(e),
      'Syntax error: improper definition name "foo$limit"'
      ' (file "tmp0d5f6e10.phil", line 1)')
  else: raise Exception_expected
  intercepted = []
  def custom_processor(arg):
    intercepted.append(arg)
    return True
  args = ["tmp0d5f6e10.phil", "lmit=3"]
  user_phils = itpr_bar.process(args=args, custom_processor=custom_processor)
  assert len(user_phils) == 0
  assert intercepted == args
  intercepted = []
  def custom_processor(arg):
    if (not os.path.isfile(arg)):
      return False
    intercepted.append(arg)
    return True
  try: itpr_bar.process(args=args, custom_processor=custom_processor)
  except Sorry, e:
    assert not show_diff(str(e),
      'Uninterpretable command line argument: "lmit=3"')
  else: raise Exception_expected
  assert intercepted == args[:1]
  user_phil = itpr_bar.process_and_fetch(args=["limit=12"])
  assert not anchored_block_show_diff(user_phil.as_str(), 9, """\
    limit = 12
""")
  #
  pcl = phil.process_command_line(
    args=["bar.max=4943"], master_string=master_string)
  assert pcl.parse is phil.parse
  assert not show_diff(pcl.master.as_str(), master_phil.as_str())
  s = StringIO()
  pcl.show(out=s)
  assert not anchored_block_show_diff(s.getvalue(), 7, """\
  max = 4943
""")
  assert pcl.remaining_args == []
  pcl = phil.process_command_line(
    args=["892c8632"], master_string=master_string)
  assert pcl.remaining_args == ["892c8632"]
  master_phil = phil.parse("""
sites = None
  .type = int
""")
  def custom_processor (arg) :
    try :
      val = int(arg)
    except ValueError :
      return None
    else :
      return phil.parse("sites=%d" % val)
  pcl = master_phil.command_line_argument_interpreter()
  args = ["2"]
  user_phil = pcl.process_args(args=["2"],
    custom_processor=custom_processor)
  params = master_phil.fetch(sources=user_phil).extract()
  assert (params.sites == 2)

def exercise_choice_multi_plus_support():
  master_phil = libtbx.phil.parse("""\
  u = a b c
    .type = choice(multi=True)
    .optional = False
  """)
  # argument_interpreter used only for convenience
  # (i.e. it is not exercised here)
  cai = master_phil.command_line_argument_interpreter()
  for arg,expected_result in [
        ("u=a", "u = *a b c"),
        ("u=b", "u = a *b c"),
        ("u=c", "u = a b *c"),
        ("u=a+b", "u = *a *b c"),
        ("u=b +c", "u = a *b *c"),
        ("u=a+ c", "u = *a b *c"),
        ("u=a + b +c", "u = *a *b *c"),
        ("u=+a + b +c", "u = *a *b *c"),
        ("u=+a ++ b +c", "u = a b c"),
        ("u=a + b + *c", "u = a b *c"),
        ("u=a + b + *C", "u = a b *c"),
        ("u=a + b + 'c'", "u = a b c")]:
    work_params = master_phil.fetch(source=cai.process(arg=arg))
    assert not show_diff(work_params.as_str(), expected_result+"\n")
  for arg,err in [("u=a+d", "d"), ("u=e + b", "e")]:
    try: master_phil.fetch(source=cai.process(arg=arg))
    except Sorry, e:
      assert not show_diff(str(e), """\
Not a possible choice for u: %s (command line argument, line 1)
  Possible choices are:
    a
    b
    c""" % err)
    else: raise Exception_expected
  for val in ["++a", "a++", "a++b", "a+b+"]:
    try: master_phil.fetch(source=cai.process(arg="u="+val))
    except Sorry, e:
      assert not show_diff(str(e), """\
Not a possible choice for u: %s (command line argument, line 1)
  Possible choices are:
    a
    b
    c""" % val)
    else: raise Exception_expected

def exercise_scope_call():
  try: phil.parse("""\
s
  .call=libtbx.phil.tst.scope_call_none
{
}
""")
  except AttributeError, e:
    assert str(e) == 'scope "s" .call: object "scope_call_none" not found in' \
      ' module "libtbx.phil.tst" (input line 2)'
  else: raise Exception_expected
  #
  try: phil.parse("""\
s
  .call=libtbx.phil.tst.scope_call_not_callable
{
}
""")
  except TypeError, e:
    assert str(e) \
      == 'scope "s" .call: "libtbx.phil.tst.scope_call_not_callable"' \
        ' is not a callable Python object (input line 2)'
  else: raise Exception_expected
  #
  try: phil.parse("""\
s
  .call=libtbx.phil.tst.scope_call_func(
{
}
""")
  except RuntimeError, e:
    assert str(e).startswith(
      'scope "s" .call=libtbx.phil.tst.scope_call_func(: ')
    assert str(e).endswith(' (input line 2)')
  else: raise Exception_expected
  #
  try: phil.parse("""\
s
  .call=libtbx.phil.tst.scope_call_func(a=b=c)
{
}
""")
  except RuntimeError, e:
    assert str(e) == 'scope "s" .call=libtbx.phil.tst.scope_call_func(a=b=c):'\
      ' SyntaxError: invalid syntax (line 1) (input line 2)'
  else: raise Exception_expected
  #
  master = phil.parse("""\
s
  .call=libtbx.phil.tst.scope_call_func
{
}
t
  .call=libtbx.phil.tst.scope_call_func(a =1)
{
}
u
  .call=libtbx.phil.tst.scope_call_func (a= 1 , b=2)
{
}
v
{
}
""")
  assert not show_diff(master.as_str(attributes_level=2), """\
s
  .call = libtbx.phil.tst.scope_call_func
{
}
t
  .call = libtbx.phil.tst.scope_call_func(a=1)
{
}
u
  .call = libtbx.phil.tst.scope_call_func(a=1, b=2)
{
}
v {
}
""")
  params = master.extract()
  c = params.s()
  assert c[0] is params.s
  assert c[1] == {}
  c = params.t()
  assert c[0] is params.t
  assert c[1] == {"a": 1}
  c = params.t(x=3)
  assert c[0] is params.t
  assert c[1] == {"a": 1, "x": 3}
  c = params.t(a=3)
  assert c[0] is params.t
  assert c[1] == {"a": 3}
  c = params.u()
  assert c[0] is params.u
  assert c[1] == {"a": 1, "b": 2}
  c = params.u(x=3)
  assert c[0] is params.u
  assert c[1] == {"a": 1, "b": 2, "x": 3}
  c = params.u(a=3)
  assert c[0] is params.u
  assert c[1] == {"a": 3, "b": 2}
  try: params.u(action="raise")
  except RuntimeError, e:
    assert str(e) \
      == 'scope "u" .call=libtbx.phil.tst.scope_call_func(a=1, b=2)' \
        ' execution: ValueError: action==raise (input line 10)'
  else: raise Exception_expected
  try: params.v()
  except RuntimeError, e:
    assert str(e) == 'scope "v" is not callable.'
  else: raise Exception_expected
  #
  master = phil.parse("""\
s
  .call=libtbx.phil.tst.scope_call_class
{
}
u
  .call = libtbx . phil .tst. scope_call_class ( a = 1 , b = 2 )
{
}
""")
  assert not show_diff(master.as_str(attributes_level=2), """\
s
  .call = libtbx.phil.tst.scope_call_class
{
}
u
  .call = libtbx.phil.tst.scope_call_class(a=1, b=2)
{
}
""")
  params = master.extract()
  c = params.s()
  assert c.scope_extract is params.s
  assert c.keyword_args == {}
  c = params.u(b=3, x=4)
  assert c.scope_extract is params.u
  assert c.keyword_args == {"a": 1, "b": 3, "x": 4}
  #
  master = phil.parse("""\
s
  .call=libtbx.phil.tst.scope_call_class_object
{
}
u
  .call=libtbx.phil.tst.scope_call_class_object(a=1, b=2)
{
}
""")
  assert not show_diff(master.as_str(attributes_level=2), """\
s
  .call = libtbx.phil.tst.scope_call_class_object
{
}
u
  .call = libtbx.phil.tst.scope_call_class_object(a=1, b=2)
{
}
""")
  params = master.extract()
  c = params.s()
  assert c.scope_extract is params.s
  assert c.keyword_args == {}
  c = params.u(a=3, b=4)
  assert c.scope_extract is params.u
  assert c.keyword_args == {"a": 3, "b": 4}

def exercise_deprecation () :
  master = phil.parse("""
foo {
  bar = None
    .type = str
}
fubar = None
  .type = str
  .deprecated = True
strategy = *xyz *adp tls
  .type = choice(multi=True)
  .deprecated = True
""")
  class _showwarning (object) :
    def __init__ (self) :
      self.n = 0
      self.message = None

    def __call__ (self, message, category, *args, **kwds) :
      if (category is phil.PhilDeprecationWarning) :
        self.n += 1
        self.message = str(message)
  warn = _showwarning()
  warnings.showwarning = warn.__call__
  w0 = master.fetch()
  out = StringIO()
  w0.show(out=out, attributes_level=2)
  assert (not "fubar" in out.getvalue())
  user1 = phil.parse("""fubar=None""")
  w1 = master.fetch(source=user1)
  out = StringIO()
  w1.show(out=out, attributes_level=3)
  assert (not "fubar" in out.getvalue())
  user2 = phil.parse("""fubar=abcedf""")
  w2 = master.fetch(source=user2)
  assert (warn.n == 1)
  assert (warn.message == 'fubar is deprecated - not recommended for use.')
  out = StringIO()
  w2.show(out=out, attributes_level=3)
  assert ("fubar" in out.getvalue())
  user3 = phil.parse("""strategy=xyz+tls""")
  w2 = master.fetch(source=user3)
  assert (warn.n == 2)
  assert (warn.message == 'strategy is deprecated - not recommended for use.')

scope_call_not_callable = None

def scope_call_func(scope_extract, **keyword_args):
  if (keyword_args.get("action") == "raise"):
    raise ValueError("action==raise")
  return scope_extract, keyword_args

class scope_call_class:
  def __init__(self, scope_extract, **keyword_args):
    self.scope_extract = scope_extract
    self.keyword_args = keyword_args

class scope_call_class_object(object):
  def __init__(self, scope_extract, **keyword_args):
    self.scope_extract = scope_extract
    self.keyword_args = keyword_args

def exercise():
  exercise_string_quote_and_tokenize()
  exercise_parse_and_show()
  exercise_import_converters()
  exercise_syntax_errors()
  exercise_phil_on_off_end()
  exercise_deepcopy()
  exercise_get_without_substitution()
  exercise_nested()
  exercise_get_with_substitution()
  exercise_include()
  exercise_full_path()
  exercise_fetch()
  exercise_fetch_diff()
  exercise_extract()
  exercise_format()
  exercise_type_constructors()
  exercise_choice()
  exercise_scope_call()
  exercise_auto()
  exercise_int_and_float()
  exercise_ints_and_floats()
  exercise_definition_validate_etc()
  exercise_command_line()
  exercise_choice_multi_plus_support()
  exercise_deprecation()
  print "OK"

if (__name__ == "__main__"):
  exercise()
