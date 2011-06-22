from libtbx.command_line import find_unused_imports_crude
import sys, os

class file_clutter(object):

  def __init__(self, path, find_unused_imports=False):
    self.path = path
    self.is_executable = os.access(path, os.X_OK)
    self.dos_format = False
    self.n_tabs_or_trailing_whitespace = []
    self.n_trailing_empty_lines = 0
    self.missing_eol = False
    self.n_bare_excepts = 0
    self.unused_imports = None
    bytes = open(path, "rb").read()
    if (len(bytes) > 0):
      if (bytes[-1] != "\n"):
        self.missing_eol = True
      else:
        bytes = bytes[:-1]
      text = bytes.split("\n")
      for i, line in enumerate(text):
        if (line.endswith("\r")):
          line = line[:-1]
          self.dos_format = True
        clean_line = line.expandtabs().rstrip()
        if (clean_line != line): self.n_tabs_or_trailing_whitespace.append(i+1)
        if (len(clean_line) == 0): self.n_trailing_empty_lines += 1
        else: self.n_trailing_empty_lines = 0
      if (path.endswith(".py")):
        py_lines = bytes.splitlines()
        for line in py_lines:
          ls = line.strip()
          if (    ls.startswith("except")
              and ls[6:].strip().startswith(":")
              and not ls.endswith(" # intentional")):
            self.n_bare_excepts += 1
        if (find_unused_imports and path.endswith(".py")):
          self.unused_imports = find_unused_imports_crude.inspect(
            py_lines=py_lines)

  def is_cluttered(self, flag_x):
    return ((self.is_executable and flag_x)
            or self.dos_format
            or len(self.n_tabs_or_trailing_whitespace) > 0
            or self.n_trailing_empty_lines > 1
            or self.missing_eol)

  def has_unused_imports(self):
    return (self.unused_imports is not None and len(self.unused_imports) != 0)

  def status(self, flag_x, flag_dos_format=True):
    status = []
    def sapp(s): status.append(s)
    if (self.is_executable and flag_x
        and self.path.lower().find("command_line") < 0
        and not self.path.endswith(".csh")
        and not self.path.endswith(".sh")):
      sapp("is executable")
    if (flag_dos_format and self.dos_format):
      sapp("dos format")
    if (len(self.n_tabs_or_trailing_whitespace) > 0):
      line = "tabs or trailing whitespace=%d" % \
        len(self.n_tabs_or_trailing_whitespace)
      for cnt, i in enumerate(self.n_tabs_or_trailing_whitespace):
        line += ", #" + str(i)
        if cnt >= 9:
          break
      sapp(line)
    if (self.n_trailing_empty_lines > 1):
      sapp("trailing empty lines=%d" % self.n_trailing_empty_lines)
    if (self.missing_eol):
      sapp("missing end-of-line")
    if (self.n_bare_excepts > 0):
      sapp("bare excepts=%d" % self.n_bare_excepts)
    if (self.has_unused_imports()):
      sapp("unused imports=%d" % len(self.unused_imports))
    return ", ".join(status)

  def show(self, flag_x, flag_dos_format=True, append=None):
    status = self.status(flag_x, flag_dos_format)
    if (len(status) != 0):
      msg = "%s: %s" % (self.path, status)
      if (append is not None):
        append(msg)
      else:
        print msg

def is_text_file(file_name):
  name = file_name.lower()
  for extension in (".c", ".cpp", ".h", ".hpp", ".py", ".java", ".params",
                    ".dox", ".txt", ".html", ".csh", ".sh", ".cif"):
    if (name.endswith(extension)): return True
  return False

def gather(paths, find_unused_imports=False):
  clutter = []
  def capp():
    clutter.append(file_clutter(path, find_unused_imports))
  for path in paths:
    if (not os.path.exists(path)):
      print >> sys.stderr, "No such file or directory:", path
    elif (os.path.isfile(path)):
      capp()
    else:
      for root, dirs, files in os.walk(path):
        for f in files:
          if (is_text_file(f)):
            path = os.path.normpath(os.path.join(root, f))
            capp()
  return clutter
