class value(dict):

  def __init__(self, default_value):
    self.default_value = default_value

  def __getitem__(self, key):
    try: return dict.__getitem__(self, key)
    except: pass
    val = self.default_value
    dict.__setitem__(self, key, val)
    return val

class factory(dict):

  def __init__(self, default_factory):
    self.default_factory = default_factory

  def __getitem__(self, key):
    try: return dict.__getitem__(self, key)
    except: pass
    val = self.default_factory()
    dict.__setitem__(self, key, val)
    return val
