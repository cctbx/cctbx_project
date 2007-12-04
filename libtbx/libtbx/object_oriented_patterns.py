""" A library of object-oriented patterns """

import new

class injector:
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
  class __metaclass__(type):
    def __init__(cls, classname, bases, classdict):
      for target_class in bases[1:]:
        for name, attribute in classdict.items():
          if name in ('__module__', '__doc__'): continue
          assert not hasattr(target_class, name), (
            "class %s has already attribute '%s'"
            % (target_class.__name__, name))
          setattr(target_class, name, attribute)


class memoize(object):
  """ Memoize the result returned by a function """

  def __init__(self, func):
    self.cached = {}
    self.func = func

  def __call__(self, *args, **kwds):
    if kwds:
      return self.func(*args, **kwds)
    try:
      return self.cached[args]
    except KeyError:
      self.cached[args] = result = self.func(*args)
      return result
    except TypeError:
      return self.func(*args)


class memoize_method(object):
  """ Memoize the result returned by a bound method """

  def __init__(self, meth):
    self.cache = '_memoized_%s' % meth.__name__
    self.meth = meth

  def __get__(self, obj, type=None):
    if obj is None:
      return self
    try:
      return getattr(obj, self.cache)
    except AttributeError:
      func = new.instancemethod(self.meth, obj, type)
      memoized = memoize(func)
      setattr(obj, self.cache, memoized)
      return memoized
