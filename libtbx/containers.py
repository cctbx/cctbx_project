from __future__ import absolute_import, division, print_function

import collections
from collections import OrderedDict  # special import
try:
  from collections.abc import MutableSet
except ImportError:
  from collections import MutableSet

class OrderedSet(MutableSet):
  """
  http://code.activestate.com/recipes/576694/ rev9
  recommended replacement: https://pypi.python.org/pypi/orderedset
  """

  def __init__(self, iterable=None):
    self.end = end = []
    end += [None, end, end]         # sentinel node for doubly linked list
    self.map = {}                   # key --> [key, prev, next]
    if iterable is not None:
      self |= iterable

  def __len__(self):
    return len(self.map)

  def __contains__(self, key):
    return key in self.map

  def add(self, key):
    if key not in self.map:
      end = self.end
      curr = end[1]
      curr[2] = end[1] = self.map[key] = [key, curr, end]

  def discard(self, key):
    if key in self.map:
      key, prev, next = self.map.pop(key)
      prev[2] = next
      next[1] = prev

  def __iter__(self):
    end = self.end
    curr = end[2]
    while curr is not end:
      yield curr[0]
      curr = curr[2]

  def __reversed__(self):
    end = self.end
    curr = end[1]
    while curr is not end:
      yield curr[0]
      curr = curr[1]

  def pop(self, last=True):
    if not self:
      raise KeyError('set is empty')
    key = self.end[1][0] if last else self.end[2][0]
    self.discard(key)
    return key

  def __repr__(self):
    if not self:
      return '%s()' % (self.__class__.__name__,)
    return '%s(%r)' % (self.__class__.__name__, list(self))

  def __eq__(self, other):
    if isinstance(other, OrderedSet):
      return len(self) == len(other) and list(self) == list(other)
    return set(self) == set(other)

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

  def __bool__(self):
    return bool(self._list_proxy)

  __nonzero__ = __bool__

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
