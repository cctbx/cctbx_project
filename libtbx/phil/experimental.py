# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-

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
  otherwise.  XXX Experimental!  XXX param or extract in name?
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

    val = getattr(src, attr)
    if (val is None):
      # Do not merge undefined values.
      continue

    elif (isinstance(val, phil.scope_extract)):
      # If dst has a scope with the same name, merge it recursively.
      # XXX Not sure how to handle the case when dst has a non-scope
      # item with the same name.  That can probably never occur for
      # proper phil objects?  Same applies to the scope_extract_list
      # below.
      if (hasattr(dst, attr) and
          isinstance(getattr(dst, attr), phil.scope_extract)):
        merge_params_by_key(getattr(dst, attr), val, key)
      else:
        setattr(dst, attr, val)

    elif (isinstance(val, phil.scope_extract_list)):
      # If dst has a list of scopes with the same name, merge them
      # recursively, item by item.
      if (hasattr(dst, attr) and
          isinstance(getattr(dst, attr), phil.scope_extract_list)):
        dst_list = getattr(dst, attr)
        for src_item in val:
          if (not _merge_param_into_list(dst_list, src_item, key)):
            dst_list.append(src_item)
      else:
        setattr(dst, attr, val)

    else:
      # Merge non-scope values.
      setattr(dst, attr, val)

  return True
