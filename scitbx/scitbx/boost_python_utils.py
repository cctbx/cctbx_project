from scitbx.array_family import flex

# see boost/libs/python/doc/tutorial/doc/quickstart.txt

boost_python_metaclass = flex.double.__class__

class injector(object):

  class __metaclass__(boost_python_metaclass):

    def __init__(self, name, bases, dict):
      for b in bases:
        if type(b) not in (self, type):
          for k,v in dict.items():
            if (k in ("__init__", "__del__", "__module__", "__file__")):
              continue
            setattr(b,k,v)
      return type.__init__(self, name, bases, dict)
