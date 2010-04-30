import sys
vers_info = sys.version_info[:2]
if (vers_info == (2,3)):
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
