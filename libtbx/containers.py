import sys

if (sys.version_info[:2] == (2,3)):
  MutableSet = None
else:
  import collections
  MutableSet = getattr(collections, "MutableSet", None)
if (MutableSet is None):
  class MutableSet(set):

    def add(self, value):
        """Add an element."""
        raise NotImplementedError

    def discard(self, value):
        """Remove an element.  Do not raise an exception if absent."""
        raise NotImplementedError

    def remove(self, value):
        """Remove an element. If not a member, raise a KeyError."""
        if value not in self:
            raise KeyError(value)
        self.discard(value)

    def pop(self):
        """Return the popped value.  Raise KeyError if empty."""
        it = iter(self)
        try:
            value = next(it)
        except StopIteration:
            raise KeyError
        self.discard(value)
        return value

    def clear(self):
        """This is slow (creates N new iterators!) but effective."""
        try:
            while True:
                self.pop()
        except KeyError:
            pass

    def __ior__(self, it):
        for value in it:
            self.add(value)
        return self

    def __iand__(self, it):
        for value in (self - it):
            self.discard(value)
        return self

    def __ixor__(self, it):
        if not isinstance(it, Set):
            it = self._from_iterable(it)
        for value in it:
            if value in self:
                self.discard(value)
            else:
                self.add(value)
        return self

    def __isub__(self, it):
        for value in it:
            self.discard(value)
        return self

class OrderedSet(MutableSet):
  """
  http://code.activestate.com/recipes/528878-ordered-set/
  Minor changes for Python 2.3 compatibility.
  """

  def __init__(self, iterable=None):
    self.end = end = []
    end += [None, end, end] # sentinel node for doubly linked list
    self.map = {}           # key --> [key, prev, next]
    if iterable is not None:
      self |= iterable

  def __len__(self):
    return len(self.map)

  def __contains__(self, key):
    return key in self.map

  def add(self, key):
    if key not in self.map:
      end = self.end
      PREV, NEXT = 1,2
      curr = end[PREV]
      curr[NEXT] = end[PREV] = self.map[key] = [key, curr, end]

  def discard(self, key):
    if key in self.map:
      key, prev, next = self.map.pop(key)
      PREV, NEXT = 1,2
      prev[NEXT] = next
      next[PREV] = prev

  def __iter__(self):
    end = self.end
    KEY, NEXT = 0,2
    curr = end[NEXT]
    while curr is not end:
      yield curr[KEY]
      curr = curr[NEXT]

  def __reversed__(self):
    end = self.end
    KEY, PREV = 0,1
    curr = end[PREV]
    while curr is not end:
      yield curr[KEY]
      curr = curr[PREV]

  def pop(self, last=True):
    if not self:
      raise KeyError('set is empty')
    if (last):
      for key in self.__reversed__(): break
    else:
      for key in self.__iter__(): break
    self.discard(key)
    return key

  def __repr__(self):
    if not self:
      return '%s()' % (self.__class__.__name__,)
    return '%s(%r)' % (self.__class__.__name__, list(self))

  def __eq__(self, other):
    if isinstance(other, OrderedSet):
      return len(self) == len(other) and list(self) == list(other)
    return not self.isdisjoint(other)

  def __ne__(self, other):
    return not self.__eq__(other)

  def __del__(self):
    self.clear() # remove circular references

  def __copy__(self):
    from copy import copy
    result = OrderedSet()
    for elt in self:
      result.add(elt)
    return result

  copy = __copy__

  def __deepcopy__(self, memo):
    from copy import deepcopy
    result = OrderedSet()
    for elt in self:
      result.add(deepcopy(elt, memo))
    return result


from UserDict import DictMixin
class OrderedDict(dict, DictMixin):
  """
  An equivalent recipe for the OrderedDict class introduced in Python 2.7
  from ActiveState:
  http://code.activestate.com/recipes/576693/
  With minor modifications for Python 2.3 compatibility.
  """

  def __init__(self, *args, **kwds):
    if len(args) > 1:
      raise TypeError('expected at most 1 arguments, got %d' % len(args))
    try:
      self.__end
    except AttributeError:
      self.clear()
    self.update(*args, **kwds)

  def clear(self):
    self.__end = end = []
    end += [None, end, end]         # sentinel node for doubly linked list
    self.__map = {}                 # key --> [key, prev, next]
    dict.clear(self)

  def __setitem__(self, key, value):
    if key not in self:
      end = self.__end
      curr = end[1]
      curr[2] = end[1] = self.__map[key] = [key, curr, end]
    dict.__setitem__(self, key, value)

  def __delitem__(self, key):
    dict.__delitem__(self, key)
    key, prev, next = self.__map.pop(key)
    prev[2] = next
    next[1] = prev

  def __iter__(self):
    end = self.__end
    curr = end[2]
    while curr is not end:
      yield curr[0]
      curr = curr[2]

  def __reversed__(self):
    end = self.__end
    curr = end[1]
    while curr is not end:
      yield curr[0]
      curr = curr[1]

  def popitem(self, last=True):
    if not self:
      raise KeyError('dictionary is empty')
    if last:
      key = reversed(self).next()
    else:
      key = iter(self).next()
    value = self.pop(key)
    return key, value

  def __reduce__(self):
    items = [[k, self[k]] for k in self]
    tmp = self.__map, self.__end
    del self.__map, self.__end
    inst_dict = vars(self).copy()
    self.__map, self.__end = tmp
    if inst_dict:
      return (self.__class__, (items,), inst_dict)
    return self.__class__, (items,)

  def keys(self):
    return list(self)

  setdefault = DictMixin.setdefault
  #update = DictMixin.update
  pop = DictMixin.pop
  values = DictMixin.values
  items = DictMixin.items
  iterkeys = DictMixin.iterkeys
  itervalues = DictMixin.itervalues
  iteritems = DictMixin.iteritems

  # this method copied from Python26 DictMixin source
  # needed for compatibility with Python 2.3
  def update(self, other=None, **kwargs):
    # Make progressively weaker assumptions about "other"
    if other is None:
      pass
    elif hasattr(other, 'iteritems'):  # iteritems saves memory and lookups
      for k, v in other.iteritems():
        self[k] = v
    elif hasattr(other, 'keys'):
      for k in other.keys():
        self[k] = other[k]
    else:
      for k, v in other:
        self[k] = v
    if kwargs:
      self.update(kwargs)

  def __repr__(self):
    if not self:
      return '%s()' % (self.__class__.__name__,)
    return '%s(%r)' % (self.__class__.__name__, self.items())

  def copy(self):
      return self.__class__(self)

  def fromkeys(cls, iterable, value=None):
    d = cls()
    for key in iterable:
      d[key] = value
    return d
  fromkeys = classmethod(fromkeys)

  def __eq__(self, other):
    if isinstance(other, OrderedDict):
      return len(self)==len(other) and self.items() == other.items()
    return dict.__eq__(self, other)

  def __ne__(self, other):
    return not self == other

class deque_template(object):

  __slots__ = ('_list_proxy', '_set_proxy')

  def __init__(self):
    self._list_proxy = self.__class__.list_proxy_type()
    if self.__class__.set_proxy_type is not None:
      self._set_proxy = self.__class__.set_proxy_type()
    else:
      self._set_proxy = self._list_proxy

  def push(self, item):
    self._list_proxy.append(item)
    if self._set_proxy is not self._list_proxy:
      self._set_proxy.add(item)
    return self

  def __nonzero__(self):
    return bool(self._list_proxy)

  def __contains__(self, item):
    return item in self._set_proxy


class stack(deque_template):
  list_proxy_type = list
  set_proxy_type  = None

  def pull(self):
    result = self._l.pop()
    if self._set_proxy is not self._list_proxy:
      self._set_proxy.discard(result)
    return result


class hashed_stack(stack):
  set_proxy_type = set


if (sys.version_info[:2] > (2,3)):
  import collections

  class queue(deque_template):
    list_proxy_type = collections.deque
    set_proxy_type  = None

    def pull(self):
      result = self._list_proxy.popleft()
      if self._set_proxy is not self._list_proxy:
        self._set_proxy.discard(result)
      return result

else:

  class queue(deque_template):
    list_proxy_type = list
    set_proxy_type  = None

    def pull(self):
      result = self._list_proxy.pop(0)
      if self._set_proxy is not self._list_proxy:
        try:
          self._set_proxy.remove(result)
        except ValueError:
          pass
      return result


class hashed_queue(queue):
  set_proxy_type = set
