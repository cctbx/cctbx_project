from __future__ import absolute_import, division, print_function
try:
  import boost_adaptbx.boost.python as bp
except Exception:
  ext = None
else:
  ext = bp.import_ext("fable_ext", optional=True)
from six.moves import range

# compare with fem/utils/string.hpp
def py_fem_utils_unsigned_integer_scan(code, start=0, stop=-1):
  i = start
  while (i < stop):
    c = code[i]
    if (not c.isdigit()): break
    i += 1
  if (i == start): return -1
  return i

# compare with ext.cpp
def py_ext_get_code_stop(code, stop):
  len_code = len(code)
  if (stop < 0): return len_code
  assert stop <= len_code
  return stop

# compare with ext.cpp
def py_unsigned_integer_scan(code, start=0, stop=-1):
  return py_fem_utils_unsigned_integer_scan(
    code=code, start=start, stop=py_ext_get_code_stop(code, stop))

# compare with ext.cpp
def py_floating_point_scan_after_exponent_char(code, start=0, stop=-1):
  code_stop = py_ext_get_code_stop(code=code, stop=stop)
  i = start
  if (i < code_stop):
    c = code[i]
    if (c == '+' or c == '-'):
      i += 1
    return py_unsigned_integer_scan(code=code, start=i, stop=stop)
  return -1

# compare with ext.cpp
def py_floating_point_scan_after_dot(code, start=0, stop=-1):
  code_stop = py_ext_get_code_stop(code=code, stop=stop)
  i = py_unsigned_integer_scan(code=code, start=start, stop=stop)
  if (i < 0): i = start
  if (i < code_stop):
    c = code[i]
    if (c == 'e' or c == 'd'):
      return py_floating_point_scan_after_exponent_char(
        code=code, start=i+1, stop=stop)
  return i

# compare with ext.cpp
def py_identifier_scan(code, start=0, stop=-1):
  code_stop = py_ext_get_code_stop(code=code, stop=stop)
  i = start
  if (i < code_stop):
    c = code[i]; i += 1
    if ((c < 'a' or c > 'z') and c != '_'): return -1
    while (i < code_stop):
      c = code[i]; i += 1
      if (    (c < 'a' or c > 'z')
          and (c < '0' or c > '9') and c != '_'): return i-1
    return i
  return -1

def py_find_closing_parenthesis(code, start=0, stop=-1):
  code_stop = py_ext_get_code_stop(code=code, stop=stop)
  n_inner = 0
  for i in range(start, code_stop):
    c = code[i]
    if (c == ')'):
      if (n_inner == 0): return i
      n_inner -= 1
    elif (c == '('):
      n_inner += 1
  return -1

if (ext is not None):
  from fable_ext import *
else:
  unsigned_integer_scan = py_unsigned_integer_scan
  floating_point_scan_after_exponent_char = \
    py_floating_point_scan_after_exponent_char
  floating_point_scan_after_dot = py_floating_point_scan_after_dot
  identifier_scan = py_identifier_scan
  find_closing_parenthesis = py_find_closing_parenthesis

class SemanticError(Exception): pass
