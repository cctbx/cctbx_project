from __future__ import absolute_import, division, print_function


def exercise_unused_imports():
  from libtbx.python_code_parsing import unused_imports, imported_name
  unused = unused_imports(unused_imports_test_case_1_1)
  assert unused.names == set((
    'foo', 'far', 'close', 'baz', 'barbar', 'bozboz' ))
  unused = unused_imports(unused_imports_test_case_1_2)
  assert unused.names == set((
    'far', 'close', 'baz', 'barbar', 'bozboz' ))
  unused = unused_imports(unused_imports_test_case_1_3)
  assert unused.names == set((
    'foo', 'far', 'close', 'barbar', 'bozboz' ))
  unused = unused_imports(unused_imports_test_case_1_4)
  assert unused.names == set((
    'far', 'close', 'barbar', 'bozboz' ))
  unused = unused_imports(unused_imports_test_case_1_5)
  assert unused.names == set((
    'far', 'close', 'baz', 'barbar', 'bozboz' ))
  unused = unused_imports(unused_imports_test_case_2)
  assert not unused
  unused = unused_imports(unused_imports_test_case_3)
  assert set(unused) == set(( imported_name('foo', lineno=1), ))
  unused = unused_imports(unused_imports_test_case_4)
  assert set(unused) == set(( imported_name('foo', lineno=4), ))
  unused = unused_imports(unused_imports_test_case_5)
  assert set(unused) == set(( imported_name('foo', lineno=4), ))
  unused = unused_imports(unused_imports_test_case_6)
  assert set(unused) == set(( imported_name('bar', lineno=1), ))
  unused = unused_imports(unused_imports_test_case_7)
  assert set(unused) == set(( imported_name('bar', lineno=1),
                              imported_name('foo', lineno=6), ))
  unused = unused_imports(unused_imports_test_case_7_bis)
  assert set(unused) == set(( imported_name('bar', lineno=1), ))
  unused = unused_imports(unused_imports_test_case_8)
  assert set(unused) == set(( imported_name('foo', lineno=1), ))
  unused = unused_imports(unused_imports_test_case_9)
  assert not unused
  unused = unused_imports(unused_imports_test_case_10)
  assert set(unused) == set(( imported_name('foo', lineno=10), ))
  unused = unused_imports(unused_imports_test_case_11)
  assert unused
  unused = unused_imports(unused_imports_test_case_11,
                          ignored_imports=('libtbx.load_env',))
  assert not unused
  unused = unused_imports(unused_imports_test_case_12)
  assert not unused
  unused = unused_imports(unused_imports_test_case_13)
  assert set(unused) == set(( imported_name('bar', lineno=1), ))
  unused = unused_imports(unused_imports_test_case_14)
  assert set(unused) == set(( imported_name('bar', lineno=1), ))
  unused = unused_imports(unused_imports_test_case_15)
  assert set(unused) == set(( imported_name('bar', lineno=1), ))
  unused = unused_imports(unused_imports_test_case_16)
  assert not unused
  unused = unused_imports(unused_imports_test_case_17)
  assert not unused
  unused = unused_imports(unused_imports_test_case_18)
  assert not unused
  unused = unused_imports(unused_imports_test_case_19)
  assert not unused
  unused = unused_imports(unused_imports_test_case_20)
  assert not unused
  unused = unused_imports(unused_imports_test_case_21)
  assert not unused
  unused = unused_imports(unused_imports_test_case_22)
  assert not unused
  unused = unused_imports(unused_imports_test_case_23)
  assert unused
  unused = unused_imports(
    unused_imports_test_case_23,
    ignore_imports_flagged_by_comments=('# import dependency',))
  assert not unused
  unused = unused_imports(unused_imports_test_case_24)
  assert not unused
  unused = unused_imports(unused_imports_test_case_25)
  assert not unused
  unused = unused_imports(unused_imports_test_case_26)
  assert not unused


unused_imports_test_case_1_header = """\
import foo
import far, near as close
from bar import baz
from foobar import barbar, bazbaz as bozboz
from buz import *
"""

unused_imports_test_case_1_1 = """\
%s
""" % unused_imports_test_case_1_header

unused_imports_test_case_1_2 = """\
%s
type(foo.z)
""" % unused_imports_test_case_1_header

unused_imports_test_case_1_3 = """\
%s
type(baz.x)
""" % unused_imports_test_case_1_header

unused_imports_test_case_1_4 = """\
%s
type(foo.z)
type(baz.x)
""" % unused_imports_test_case_1_header

unused_imports_test_case_1_5 = """\
%s
type(foo.baz)
""" % unused_imports_test_case_1_header


unused_imports_test_case_2 = """\
import foo

def f():
  return type(foo.z)
"""

unused_imports_test_case_3 = """\
import foo

def f():
  return type(bar.foo.z)
"""

unused_imports_test_case_4 = """\
import bar

def f():
  import foo
  return type(bar.foo.z)
"""

unused_imports_test_case_5 = """\
import bar

def f():
  import foo
  def g():
    return type(bar.foo.z)
  return g
"""

unused_imports_test_case_6 = """\
import bar

def f():
  import foo
  def g():
    return type(foo.z)
  return g
"""

unused_imports_test_case_7 = """\
import bar

def f():
  def g():
    return type(foo.z)
  import foo
  return g
"""

unused_imports_test_case_7_bis = """\
import bar

def f():
  def g():
    return type(foo.z)
  import foo
  def h():
    return type(foo.y)
  return g, h
"""

unused_imports_test_case_8 = """\
import foo, bar

class klass(object):

  def f(self):
    self.x = bar.z
"""

unused_imports_test_case_9 = """\
import foo, bar

class klass(object):

  attr = foo.x

  def f(self):
    self.x = bar.z
"""

unused_imports_test_case_10 = """\
import bar

class klass(object):

  attr = foo.x

  def f(self):
    self.x = bar.z

import foo
"""

unused_imports_test_case_11 = """\
import libtbx.load_env

def run():
  libtbx.env.build_options.report()
"""

unused_imports_test_case_12 = """\
import foo.bar
type(foo.bar.x)
"""

unused_imports_test_case_13 = """\
import foo, bar
type(foo.moo.maz)
"""

unused_imports_test_case_14 = """\
import foo, bar
type(foo.moo().moz)
"""

unused_imports_test_case_15 = """\
import foo, bar
type(foo.moo[5].moz)
"""

unused_imports_test_case_16 = """\
from foo import bobar
type(bobar)
"""

unused_imports_test_case_17 = """\
import foobar
text = '.'.join(foobar.z)
"""

unused_imports_test_case_18 = """\
import foz
def faaz(arg=foz.x.yyy):
  pass
"""

unused_imports_test_case_19 = """\
import foo.barr
hasattr(foo.barr, 'bok')
"""

unused_imports_test_case_20 = """\
import foo.barr
class bokka(kool, foo.barr.klass):
  pass
"""

unused_imports_test_case_21 = """\
from foo.bar import bakar
import foo.bar.bakar.baaz

type(bakar.baaz)
"""

unused_imports_test_case_22 = """\
from foo.bar import boz
def f():
  return g( boz(q) ).k
"""

unused_imports_test_case_23 = """\
from scitbx.array_family import flex # import dependency
"""

unused_imports_test_case_24 = """\
from foo.bar import boz
import baz.buz
def f():
  boz.x = 1
  baz.buz.y = 2
"""

unused_imports_test_case_25 = """\
from foo.bar import boz

def f():
  return g(a=boz.h()).x
"""

unused_imports_test_case_26 = """\
from foo.bar import boz
import fooo.baz

def f():
  return g[boz.h()].x

def h():
  return g[fooo.baz.x:1000]
"""

def exercise_old_style_classes():
  from libtbx.python_code_parsing import find_old_style_classes, old_style_class
  old_style = find_old_style_classes(old_style_class_test_case_1)
  assert old_style.names == set(('foo',))
  old_style = find_old_style_classes(old_style_class_test_case_2)
  assert old_style.names == set()
  old_style = find_old_style_classes(old_style_class_test_case_3)
  assert old_style.names == set(('bar',))
  old_style = find_old_style_classes(old_style_class_test_case_4)
  assert set(old_style) == set((old_style_class(name='foo', lineno=4),
                                old_style_class(name='bar', lineno=1)))
  old_style = find_old_style_classes(old_style_class_test_case_5)
  assert old_style.names == set(('foo',))

old_style_class_test_case_1 = """\
class foo:
  pass
"""

old_style_class_test_case_2 = """\
class foo(object):
  pass
"""

old_style_class_test_case_3 = """\
class bar:
  pass

class foo(bar):
  pass
"""

old_style_class_test_case_4 = """\
class bar:
  pass

class foo():
  pass
"""

old_style_class_test_case_5 = """\
def bar():
  class foo:
    pass
"""

def run():
  exercise_unused_imports()
  exercise_old_style_classes()
  print('OK')

if __name__ == '__main__':
  run()
