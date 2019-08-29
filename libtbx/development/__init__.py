from __future__ import absolute_import, division, print_function

def show_pickle_sizes(obj, indent="",
    display_if_size_greater_than=10000):
  """
  Function for debugging issues with pickle file size, e.g. detecting which
  objects are contributing to bloat.  This is very approximate, and somewhat
  misleading since child objects may appear to be larger than their parents,
  probably for some reason involving duplicated references.
  """
  from six.moves.cPickle import dumps
  if hasattr(obj, "__slots__"):
    attr_names = obj.__slots__
  else :
    attr_names = list(obj.__dict__.keys())
  for attr in attr_names :
    child_obj = getattr(obj, attr, None)
    if (child_obj is not None):
      pkl = dumps(child_obj)
      n_bytes = len(pkl)
      if (n_bytes > display_if_size_greater_than):
        print(indent+attr, n_bytes)
        if isinstance(child_obj, object) and hasattr(child_obj, "__dict__"):
          show_pickle_sizes(child_obj, indent+"  ",
            display_if_size_greater_than=display_if_size_greater_than)
        elif isinstance(child_obj, list):
          for item in child_obj :
            n_bytes_item = len(dumps(item))
            if (n_bytes_item > display_if_size_greater_than):
              print(indent+"  "+type(item).__name__, n_bytes_item)
              if isinstance(item, object) and hasattr(item, "__dict__"):
                show_pickle_sizes(item, indent+"    ",
                  display_if_size_greater_than=display_if_size_greater_than)
