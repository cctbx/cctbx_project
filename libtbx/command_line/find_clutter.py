from __future__ import absolute_import, division, print_function
import sys
from libtbx.file_clutter import gather

def run(args):
  flag_x = False
  flag_ni = False
  flag_dos_format = True
  flag_indentation = False
  verbose = False
  #
  only_whitespace = False
  only_dos = False
  only_future = False
  flag_absolute_import = False
  flag_print_function = False
  #
  paths = []
  for arg in args:
    if (arg == "-x"):
      flag_x = True
    elif (arg == "-ni"):
      flag_ni = True
    elif (arg == "-ndos"):
      flag_dos_format = False
    elif (arg == "--verbose") or (arg == '-v'):
      verbose = True
    elif (arg == "--indentation"):
      flag_indentation = True
    elif (arg == "--only_whitespace"):
      only_whitespace = True
    elif (arg == "--only_dos"):
      only_dos = True
    elif (arg == "--only_future"):
      only_future = True
    elif (arg == "--absolute_import"):
      flag_absolute_import = True
    elif (arg == "--print_function"):
      flag_print_function = True
    else:
      paths.append(arg)
  if (len(paths) == 0): paths = ["."]
  n_is_cluttered = 0
  n_bare_excepts = 0
  n_has_unused_imports = 0
  message_lines = []
  n_missing_from_future_import_division = 0
  n_too_many_from_future_import_division = 0
  n_missing_from_future_import_absolute_import = 0
  n_too_many_from_future_import_absolute_import = 0
  n_missing_from_future_import_print_function = 0
  n_too_many_from_future_import_print_function = 0
  n_bad_indentation = 0
  for info in gather(paths=paths, find_unused_imports=not flag_ni,
      find_bad_indentation=flag_indentation, flag_absolute_import=flag_absolute_import,
      flag_print_function=flag_print_function):
    if (info.is_cluttered(flag_x=flag_x)):
      n_is_cluttered += 1
    if (info.n_bare_excepts > 0):
      n_bare_excepts += info.n_bare_excepts
    if (info.has_unused_imports()):
      n_has_unused_imports += 1
    if info.n_from_future_import_division == 0:
      n_missing_from_future_import_division += 1
    elif info.n_from_future_import_division > 1:
      n_too_many_from_future_import_division += 1
    if info.n_from_future_import_absolute_import == 0:
      n_missing_from_future_import_absolute_import += 1
    elif info.n_from_future_import_absolute_import > 1:
      n_too_many_from_future_import_absolute_import += 1
    if info.n_from_future_import_print_function == 0:
      n_missing_from_future_import_print_function += 1
    elif info.n_from_future_import_print_function > 1:
      n_too_many_from_future_import_print_function += 1
    if (info.bad_indentation is not None) and (flag_indentation):
      n_bad_indentation += 1
    info.show(
      flag_x=flag_x,
      flag_dos_format=flag_dos_format,
      append=message_lines.append,
      flag_indentation=flag_indentation,
      verbose=verbose)
  please_use = []
  if (n_is_cluttered != 0):
    please_use.append("libtbx.clean_clutter")
  if only_whitespace:
    def _is_whitespace(s):
      if s.find("tabs or trailing")>-1: return True
      return False
    message_lines = list(filter(_is_whitespace, message_lines))
  elif only_dos:
    def _is_dos(s):
      if s.find("dos format")>-1: return True
      return False
    message_lines = list(filter(_is_dos, message_lines))
  elif only_future:
    def _is_future(s):
      if s.find("from __future__")>-1: return True
      return False
    message_lines = list(filter(_is_future, message_lines))
  else:
    if (n_has_unused_imports != 0):
      please_use.append("libtbx.find_unused_imports_crude")
    # if n_missing_from_future_import_division:
      # please_use.append('libtbx.add_from_future_import_division')
    if (len(please_use) != 0):
      message_lines.append("")
      message_lines.append(
        "*** To clean up please use: %s ***" % ", ".join(please_use))
    if (n_bare_excepts > 0):
      message_lines.append("")
      message_lines.extend("""\
  *** Please change bare excepts: ***
        Usually best:
          except Exception:
        Rarely necessary:
          except: # intentional
  """.splitlines())
    if (n_bad_indentation != 0):
      message_lines.append("")
      message_lines.append("*** Please fix indentation in a text editor ***")
  if (len(message_lines) != 0):
    print()
    print("\n".join(message_lines))
    print()
    return (1)
  return (0)

if (__name__ == "__main__"):
  sys.exit(run(sys.argv[1:]))
