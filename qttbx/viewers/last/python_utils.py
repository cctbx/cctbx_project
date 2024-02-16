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


def flatten_dict(d, parent_key='', sep='.'):
  items = []
  for k, v in d.items():
    new_key = f"{parent_key}{sep}{k}" if parent_key else k
    if isinstance(v, dict):
      items.extend(flatten_dict(v, new_key, sep=sep).items())
    else:
      items.append((new_key, v))
  return dict(items)

