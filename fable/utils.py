class unique_list(object):

  __slots__ = ["value_list", "value_set"]

  def __init__(O):
    O.value_list = []
    O.value_set = set()

  def __contains__(O, value):
    return value in O.value_set

  def append(O, value):
    if (value not in O.value_set):
      O.value_list.append(value)
      O.value_set.add(value)

class keyed_lists(object):

  __slots__ = ["indices_by_key", "keys", "lists"]

  def __init__(O):
    O.indices_by_key = {}
    O.keys = []
    O.lists = []

  def get(O, key):
    class undef(object): pass
    i = O.indices_by_key.get(key, undef)
    if (i is undef):
      O.indices_by_key[key] = len(O.keys)
      O.keys.append(key)
      result = []
      O.lists.append(result)
      return result
    return O.lists[i]

  def items(O):
    return zip(O.keys, O.lists)
