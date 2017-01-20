from __future__ import division
from libtbx.command_line import find_unused_imports_crude
import sys, os
import re

class file_clutter(object):

  from_future_import_division_pat = re.compile(
    '^ from [ ]+ __future__ [ ]+ import [ \w,]+ division', re.VERBOSE)

  def __init__(self, path, find_unused_imports=False,
      find_bad_indentation=True):
    self.path = path
    self.is_executable = os.access(path, os.X_OK)
    self.dos_format = False
    self.n_tabs_or_trailing_whitespace = []
    self.n_trailing_empty_lines = 0
    self.missing_eol = False
    self.n_bare_excepts = 0
    self.unused_imports = None
    self.n_from_future_import_division = None
    self.bad_indentation = None
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
        self.n_from_future_import_division = 0
        py_lines = bytes.splitlines()
        for line in py_lines:
          if self.from_future_import_division_pat.search(line):
            self.n_from_future_import_division += 1
          ls = line.strip()
          if (    ls.startswith("except")
              and ls[6:].strip().startswith(":")
              and not ls.endswith(" # intentional")):
            self.n_bare_excepts += 1
        if (find_unused_imports and path.endswith(".py")):
          self.unused_imports = find_unused_imports_crude.inspect(
            py_lines=py_lines)
        if (find_bad_indentation) :
          self.bad_indentation = detect_indentation_problems(path)

  def is_cluttered(self, flag_x):
    return ((self.is_executable and flag_x)
            or self.dos_format
            or len(self.n_tabs_or_trailing_whitespace) > 0
            or self.n_trailing_empty_lines > 1
            or self.missing_eol
            or self.bad_indentation)

  def has_unused_imports(self):
    return (self.unused_imports is not None and len(self.unused_imports) != 0)

  def status(self, flag_x, flag_dos_format=True, flag_indentation=False):
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
    if self.n_from_future_import_division == 0:
      sapp("missing 'from __future__ import division'")
    elif self.n_from_future_import_division > 1:
      sapp("more than one appearance of 'from __future__ import division'")
    if (self.bad_indentation is not None) and (flag_indentation) :
      n_tab, n_space = self.bad_indentation
      sapp("non-standard indentation: %d space, %d tab" % (n_space, n_tab))
    return ", ".join(status)

  def show(self, flag_x, flag_dos_format=True, append=None, verbose=False,
      flag_indentation=False):
    status = self.status(flag_x, flag_dos_format, flag_indentation)
    if (len(status) != 0):
      msg = "%s: %s" % (self.path, status)
      if (append is not None):
        append(msg)
      else:
        print msg
      if (verbose) and (self.has_unused_imports()) :
        msg2 = "  unused imports: %s" % ", ".join(self.unused_imports)
        if (append is not None):
          append(msg2)
        else :
          print msg2

def is_text_file(file_name):
  name = file_name.lower()
  for extension in (".c", ".cpp", ".h", ".hpp", ".py", ".java", ".params",
                    ".dox", ".txt", ".html", ".csh", ".sh", ".cif", ".cc"):
    if (name.endswith(extension)): return True
  return False

def gather(paths, find_unused_imports=False, find_bad_indentation=False):
  clutter = []
  def capp():
    clutter.append(file_clutter(path, find_unused_imports,
      find_bad_indentation=find_bad_indentation))
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

def detect_indentation_problems (file_name) :
  try :
    import indent_finder
  except ImportError :
    return None
  fi = indent_finder.IndentFinder()
  fi.clear()
  fi.parse_file(file_name)
  result = fi.results()
  if (result is fi.default_result) :
    return None
  itype, ival = result
  n_tab = n_space = 0
  if (itype != "mixed") :
    if (itype == "space") :
      n_space = ival
    else :
      n_tab = ival
    if (n_tab == 0) and (n_space == 2) : # this is our "standard"
      return None
    return (n_tab, n_space)
  else :
    return ival
  return None
