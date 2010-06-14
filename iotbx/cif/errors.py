# Lists of error messages thrown during CIF validation

error_dicts = {}

def add_error_dict(language, error_dict):
  error_dicts[language.lower()] = error_dict

def get_error_dict(language):
  return error_dicts[language.lower()]

english = {

  # warnings
  1001: "Unknown data item: '%(key)s' not present in dictionary",
  1002: "Warning: case-sensitive match failure for value '%(value)s' for '%(key)s'",
  1003: "Warning: obsolete definition: '%(key)s' replaced by '%(related_item)s'",

  # Type errors
  2001: "Type error for %(key)s: '%(value)s' is not interpretable as type '%(item_type)s'",
  2002: "Type error for %(key)s: su is not allowed",

  # Enumeration errors
  2101: "Invalid enumeration value for %(key)s: '%(value)s' is outside of the permitted range: %(enum)s",
  2102: "Invalid enumeration value for %(key)s: '%(value)s' not in %(enum)s",

  # Related errors
  2201: "Both data item '%(key)s' and exclusive alternate '%(related_item)s' present in data block",
  2202: "Associated value item '%(related_item)s' missing for '%(key)s'",
  2203: "Mandatory category key item '%(key)s' for category '%(category)s' not present in data block",

  # Dependent errors
  2301: "Dependent item '%(dependent)s' not present in data block for '%(key)s'",


  # Loop errors
  2501: "Invalid loop: data item '%(key)s' cannot be declared in a looped list",
  2502: "Invalid loop: data names from multiple categories present",
  2503: "Invalid loop: value '%(value)s' present in child '%(child)s' but not found in parent '%(parent)s'",
  2504: "Invalid loop: missing parent for loop containing '%(child)s': '%(parent)s' required but not present in data block",
  2505: "Invalid loop: missing reference for loop containing '%(key)s': '%(reference)s' required but not present",
  2506: "data item '%(key)s' can only be declared in a looped list",

  }

add_error_dict("en", english)
