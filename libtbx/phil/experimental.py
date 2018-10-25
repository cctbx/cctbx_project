from __future__ import absolute_import, division, print_function

from libtbx import phil


def _merge_param_into_list(dst_list, src, key):
  """The _merge_param_into_list() helper function returns @c True if
  the phil extract @p src was merged into the extract list @p dst_list
  and @c False otherwise.
  """

  for dst in dst_list:
    if (merge_params_by_key(dst, src, key)):
      return True
  return False


def merge_params_by_key(dst, src, key):
  """The merge_params_by_key() function recursively merges the phil
  extract @p src into @p dst.  Two scopes with attributes named @p key
  are merged only if the values of @p key are equal.  The function
  returns @c True if @p src was merged into @p dst, and @c False
  otherwise.
  """

  # Base case: if dst and src both have keys, they must match.
  if (hasattr(dst, key) and
      hasattr(src, key) and
      getattr(dst, key) != getattr(src, key)):
    return False

  for attr in dir(src):
    # Ignore "private" attributes.
    if (len(attr) <= 0 or attr[0] == "_"):
      continue

    src_val = getattr(src, attr)
    if (src_val is None):
      # Do not merge undefined values.
      continue

    elif (isinstance(src_val, phil.scope_extract)):
      # If dst has a scope with the same name, merge it recursively.
      if (hasattr(dst, attr)):
        dst_val = getattr(dst, attr)
        assert isinstance(dst_val, phil.scope_extract)
        merge_params_by_key(dst_val, src_val, key)
      else:
        setattr(dst, attr, src_val)

    elif (isinstance(src_val, phil.scope_extract_list)):
      # If dst has a list of scopes with the same name, merge them
      # recursively, item by item.
      if (hasattr(dst, attr)):
        dst_list = getattr(dst, attr)
        assert isinstance(dst_list, phil.scope_extract_list)
        for src_item in src_val:
          if (not _merge_param_into_list(dst_list, src_item, key)):
            dst_list.append(src_item)
      else:
        setattr(dst, attr, src_val)

    else:
      # Merge non-scope values.
      setattr(dst, attr, src_val)

  return True
