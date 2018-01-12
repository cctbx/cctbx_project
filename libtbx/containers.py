from __future__ import absolute_import, division

import collections
MutableSet = collections.MutableSet
OrderedDict = collections.OrderedDict

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

class queue(deque_template):
  list_proxy_type = collections.deque
  set_proxy_type  = None

  def pull(self):
    result = self._list_proxy.popleft()
    if self._set_proxy is not self._list_proxy:
      self._set_proxy.discard(result)
    return result

class hashed_queue(queue):
  set_proxy_type = set
