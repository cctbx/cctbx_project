import sys
from libtbx.command_line.file_clutter import gather

def run(args):
  flag_x = False
  flag_ni = False
  flag_dos_format = True
  paths = []
  for arg in args:
    if (arg == "-x"):
      flag_x = True
    elif (arg == "-ni"):
      flag_ni = True
    elif (arg == "-ndos"):
      flag_dos_format = False
    else:
      paths.append(arg)
  if (len(paths) == 0): paths = ["."]
  n_is_cluttered = 0
  n_bare_excepts = 0
  n_has_unused_imports = 0
  message_lines = []
  for info in gather(paths=paths, find_unused_imports=not flag_ni):
    if (info.is_cluttered(flag_x=flag_x)):
      n_is_cluttered += 1
    if (info.n_bare_excepts > 0):
      n_bare_excepts += info.n_bare_excepts
    if (info.has_unused_imports()):
      n_has_unused_imports += 1
    info.show(
      flag_x=flag_x,
      flag_dos_format=flag_dos_format,
      append=message_lines.append)
  please_use = []
  if (n_is_cluttered != 0):
    please_use.append("libtbx.clean_clutter")
  if (n_has_unused_imports != 0):
    please_use.append("libtbx.find_unused_imports_crude")
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
  if (len(message_lines) != 0):
    print
    print "\n".join(message_lines)
    print
    return (1)
  return (0)

if (__name__ == "__main__"):
  sys.exit(run(sys.argv[1:]))
