from __future__ import absolute_import, division, print_function
from libtbx import easy_run

def demangle(symbol):
  """ Return the demangled C++ function name corresponding to symbol
      on platforms featuring c++filt. On those which do not, return symbol
      untouched
  """
  try:
    result = easy_run.fully_buffered("c++filt %s" % symbol)\
                     .raise_if_errors()\
                     .stdout_lines[0]
  except AssertionError:
    result = symbol
  return result

if __name__ == '__main__':
  def exercise_demangle():
    assert (demangle("__ZN5iotbx2xd17master_file_input8on_error"
                     "EPKcRKNS_10flex_bison8locationEPv")
            == "iotbx::xd::master_file_input::on_error"
               "(char const*, iotbx::flex_bison::location const&, void*)")

  exercise_demangle()
