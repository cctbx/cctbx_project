class DotDict(dict):
  def __getattr__(self, name):
    if name in self:
      value = self[name]
      if isinstance(value, dict):
        return DotDict(value)
      return value
    raise AttributeError(f"No such attribute: {name}")

  def __setattr__(self, name, value):
    self[name] = value

  def __delattr__(self, name):
    if name in self:
      del self[name]
    else:
      raise AttributeError(f"No such attribute: {name}")

  def __repr__(self):
    return f"<DotDict {super().__repr__()}>"
