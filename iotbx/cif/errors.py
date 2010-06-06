# Lists of error messages thrown during CIF validation

error_dicts = {}

def add_error_dict(language, error_dict):
  error_dicts[language.lower()] = error_dict

def get_error_dict(language):
  return error_dicts[language.lower()]

english = {

  # warnings
  1001: "Unknown data item: '%(key)s' not present in dictionary",

  # Type errors
  2001: "Type error for %(key)s: '%(value)s' is not interpretable as type '%(data_type)s'",
  2002: "Type error for %(key)s: su is not allowed",

  # Enumeration errors
  2101: "Invalid enumeration value for %(key)s: '%(value)s' is outside of the permitted range: %(enum)s",
  2102: "Invalid enumeration value for %(key)s: '%(value)s' not in %(enum)s",

  # Related errors
  2201: "Both data item '%(key)s' and exclusive alternate '%(related_item)s' present in data block",

  # Loop errors
  2501: "Invalid loop: data item '%(key)s' cannot be declared in a looped list",
  2502: "Invalid loop: data names from multiple categories present",
  2503: "Invalid loop: value '%(value)s' present in _list_link_child '%(child)s' but not found in _list_link_parent '%(parent)s'",
  2504: "Invalid loop: missing _list_parent for loop containing '%(child)s': '%(parent)s' required but not present in data block",
  2505: "Invalid loop: missing _list_reference for loop containing '%(key)s': '%(reference)s' required but not present",
  2506: "data item '%(key)s' can only be declared in a looped list",

  }

add_error_dict("en", english)
