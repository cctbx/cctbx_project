

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
  assert unused


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
print foo.z
""" % unused_imports_test_case_1_header

unused_imports_test_case_1_3 = """\
%s
print baz.x
""" % unused_imports_test_case_1_header

unused_imports_test_case_1_4 = """\
%s
print foo.z
print baz.x
""" % unused_imports_test_case_1_header

unused_imports_test_case_1_5 = """\
%s
print foo.baz
""" % unused_imports_test_case_1_header


unused_imports_test_case_2 = """\
import foo

def f():
  print foo.z
"""

unused_imports_test_case_3 = """\
import foo

def f():
  print bar.foo.z
"""

unused_imports_test_case_4 = """\
import bar

def f():
  import foo
  print bar.foo.z
"""

unused_imports_test_case_5 = """\
import bar

def f():
  import foo
  def g():
    print bar.foo.z
  return g
"""

unused_imports_test_case_6 = """\
import bar

def f():
  import foo
  def g():
    print foo.z
  return g
"""

unused_imports_test_case_7 = """\
import bar

def f():
  def g():
    print foo.z
  import foo
  return g
"""

unused_imports_test_case_7_bis = """\
import bar

def f():
  def g():
    print foo.z
  import foo
  def h():
    print foo.y
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
print foo.bar.x
"""

unused_imports_test_case_13 = """\
import foo, bar
print foo.moo.maz
"""

unused_imports_test_case_14 = """\
import foo, bar
print foo.moo().moz
"""

unused_imports_test_case_15 = """\
import foo, bar
print foo.moo[5].moz
"""

unused_imports_test_case_16 = """\
from foo import bobar
print bobar
"""

unused_imports_test_case_17 = """\
import foobar
print '.'.join(foobar.z)
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

print bakar.baaz
"""

unused_imports_test_case_22 = """\
from foo.bar import boz
def f():
  print g( boz(q) ).k
"""

unused_imports_test_case_23 = """\
from scitbx.array_family import flex # import dependency
"""

def run():
  exercise_unused_imports()
  print 'OK'

if __name__ == '__main__':
  run()
