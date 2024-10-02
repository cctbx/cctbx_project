""" A library of object-oriented patterns """
from __future__ import absolute_import, division, print_function

from builtins import object

import sys
import weakref
from types import MethodType
import six


class meta_injector(type):

  def __init__(cls, name, bases, namespace, **kwds):
    for target_class in bases[1:]:
      for name, attribute in namespace.items():
        if name in ('__module__', '__doc__'):
          continue
        if sys.hexversion >= 0x03000000:
          if name in ('__qualname__',):
            continue
        if sys.version_info.major == 3 and sys.version_info.minor >= 13:
          if name in ('__firstlineno__', '__static_attributes__'):
            continue
        assert not hasattr(target_class, name), (
          "class %s already has attribute '%s'"
          % (target_class.__name__, name))
        setattr(target_class, name, attribute)


class injector(six.with_metaclass(meta_injector)):
  """ Injection of new methods into an existing class

  * synopsis *

  class some_descriptive_text(injector, some_existing_class,
                                        another_existing_class, ...):

    def new_method(self, ...): ...

    def another_new_method(self, ...): ...

  obj = some_existing_class(...)
  obj.new_method(...)
  obj = another_existing_class(...)
  obj.new_method(...)
  # class 'some_descriptive_text' is an empty shell with no use by itself.

  * motivation *

  The traditional way to add methods to an existing class involves typing
  three times the same information:

  def new_method(self, ...): ...
  existing_class.new_method = new_method

  or to defer the naming to after the definition:

  def foo(self, ...): ...
  existing_class.new_method = foo

  A bit of metaclass trickery results in a cleaner syntax.
  """


class memoize(object):
  """ Memoize the result returned by a function """

  def __init__(self, func):
    self.cached = {}
    self.func = func
    self.__doc__ = self.func.__doc__

  def __call__(self, *args):
    try:
      return self.cached[args]
    except KeyError:
      self.cached[args] = result = self.func(*args)
      return result
    except TypeError:
      return self.func(*args)


class memoize_method(object):
  """ Memoize the result returned by a bound method.
      This is to be used with immutable objects only.
  """

  def __init__(self, meth):
    self.cache = '_memoized_%s' % meth.__name__
    self.meth = meth

  def __get__(self, obj, type=None):
    if obj is None:
      return self
    try:
      return getattr(obj, self.cache)
    except AttributeError:
      # We use weakref.proxy to break the following cycle
      # obj._memoized_xxx.func.im_self is obj
      # It's always better to enable reference counting to collect
      # unreachable object as soon as they become so instead of relying
      # on a later gc collection.
      try:
        # Python 2: 3 arguments
        func = MethodType(self.meth, weakref.proxy(obj), type)
      except TypeError:
        # Python 3: 2 arguments
        func = MethodType(self.meth, weakref.proxy(obj))
      memoized = memoize(func)
      setattr(obj, self.cache, memoized)
      return memoized


class null(object):

  def __init__(self, *args, **kwds): pass

  def __getattr__(self, a): return self
  def __setattr__(self, a, v): return self
  def __delattr__(self, a): return self
  def __getnewargs__(self):  return () # allow for pickling during multiprocessing
  def __call__(self, *args, **kwds): return self

  def __getitem__(self, i): return self
  def __setitem__(self, i, v): return self
  def __delitem__(self, i): return self

  def __repr__(self): return 'null()'

  def __bool__(self): return False


class proxy(object):

  def __init__(self, subject):
    self.subject = subject

  def __getattr__(self, attr):
    return getattr(self.subject, attr)


class journal_mixin(object):
  """ An easy way to store the history of an attribute as it changes
      through the course of a routine.
  """

  __journal__ = []
  __journal_suffix__ = "_history"

  def __getattr__(self, name):
    if name in self.__journal__:
      key = name+self.__journal_suffix__
    else:
      key = name
    if key in self.__dict__:
      return self.__dict__[key][-1]
    else: raise AttributeError(name)

  def __setattr__(self, name, value):
    if name in self.__journal__:
      key = name+self.__journal_suffix__
      if key not in self.__dict__:
        self.__dict__[key] = [value]
      else:
        self.__dict__[key].append(value)
    else:
      self.__dict__[name] = value

  def __delattr__(self, name):
    if name in self.__journal__:
      key = name+self.__journal_suffix__
    else:
      key = name
    del self.__dict__[key]


def calculate_state(holder):

  holder.data = holder.calculation()
  holder.state = retrieve_state
  return holder.data


def retrieve_state(holder):

  return holder.data


class lazy_initialization(object):
  """
  Defers initialization until the value is accessed (state pattern)
  """

  def __init__(self, func, *args, **kwargs):

    import functools
    self.calculation = functools.partial( func, *args, **kwargs )
    self.data = None
    self.state = calculate_state


  def __call__(self):

    return self.state( holder = self )
