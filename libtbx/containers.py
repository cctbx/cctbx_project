""" Minimalistic (and efficient) wrappers around Python standard containers
"""

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


try:
  import collections

  class queue(deque_template):
    list_proxy_type = collections.deque
    set_proxy_type  = None

    def pull(self):
      result = self._list_proxy.popleft()
      if self._set_proxy is not self._list_proxy:
        self._set_proxy.discard(result)
      return result

except ImportError:

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
