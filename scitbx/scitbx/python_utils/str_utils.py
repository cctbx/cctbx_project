def overwrite_at(s, offset, replacement):
  return s[:offset] + replacement + s[offset+len(replacement):]
