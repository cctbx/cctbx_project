
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

def find_key_path(nested_dict, target_key, path=None):
  if path is None:
      path = []

  # Check if the current level of the nested dictionary contains the target key
  if target_key in nested_dict:
      return path + [target_key]

  # Recursively search for the key in dictionaries within the current level
  for key, value in nested_dict.items():
      if isinstance(value, dict):
          # If the value is a dictionary, search it for the target key
          new_path = find_key_path(value, target_key, path + [key])
          if new_path is not None:
              return new_path

  # Return None if the key was not found at the current level or in any nested dictionaries
  return None

def get_value_by_path(nested_dict, key_path):
    current_level = nested_dict
    for key in key_path:
        # Navigate deeper into the nested dictionary
        current_level = current_level[key]
    return current_level

