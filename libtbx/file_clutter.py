from __future__ import absolute_import, division, print_function
from libtbx.command_line import find_unused_imports_crude
from libtbx.utils import to_unicode
import sys, os
import re

class file_clutter(object):

  from_future_pat = re.compile(
    '^ from [ ]+ __future__ ', re.VERBOSE)
  from_future_import_division_pat = re.compile(
    r'^ from [ ]+ __future__ [ ]+ import [ \w,]+ division', re.VERBOSE)
  from_future_import_absolute_import_pat = re.compile(
    r'^ from [ ]+ __future__ [ ]+ import [ \w,]+ absolute_import', re.VERBOSE)
  from_future_import_print_function_pat = re.compile(
    r'^ from [ ]+ __future__ [ ]+ import [ \w,]+ print_function', re.VERBOSE)

  def __init__(self, path, find_unused_imports=False,
      find_bad_indentation=True, flag_absolute_import=False,
      flag_print_function=False):
    self.path = path
    self.is_executable = os.access(path, os.X_OK)
    self.dos_format = False
    self.n_tabs_or_trailing_whitespace = []
    self.n_trailing_empty_lines = 0
    self.missing_eol = False
    self.n_bare_excepts = 0
    self.unused_imports = None
    self.n_from_future_import_division = -1
    self.flag_absolute_import = flag_absolute_import
    self.n_from_future_import_absolute_import = -1
    self.flag_print_function = flag_print_function
    self.n_from_future_import_print_function = -1
    self.bad_indentation = None
    self.file_should_be_empty = False

    if self.ignore_file():
      return

    bytes = to_unicode(open(path, "rb").read())
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
        self.n_from_future_import_absolute_import = 0
        self.n_from_future_import_print_function = 0
        py_lines = bytes.splitlines()
        self.file_should_be_empty = True
        for line in py_lines:
          if self.file_should_be_empty and line.strip() != '' and not self.from_future_pat.search(line):
            self.file_should_be_empty = False
          if self.from_future_import_division_pat.search(line):
            self.n_from_future_import_division += 1
          if self.from_future_import_absolute_import_pat.search(line):
            self.n_from_future_import_absolute_import += 1
          if self.from_future_import_print_function_pat.search(line):
            self.n_from_future_import_print_function += 1
          ls = line.strip()
          if (    ls.startswith("except")
              and ls[6:].strip().startswith(":")
              and not ls.endswith(" # intentional")):
            self.n_bare_excepts += 1
        if (find_unused_imports and path.endswith(".py")):
          self.unused_imports = find_unused_imports_crude.inspect(
            py_lines=py_lines)
        if (find_bad_indentation):
          self.bad_indentation = detect_indentation_problems(path)

  def ignore_file(self):
    root = os.path.dirname(self.path)
    ignore_path = os.path.join(root, ".clutterignore")
    if not os.path.exists(ignore_path):
      return False
    with open(ignore_path) as f:
      for line in f:
        if len(line) == 0: continue
        if line.startswith('#'): continue
        if os.path.basename(self.path) == line.strip():
          return True
    return False

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
    sapp = status.append
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
    if self.file_should_be_empty:
      if self.n_from_future_import_division == 0 and self.n_from_future_import_absolute_import == 0 and \
         self.n_from_future_import_print_function:
        sapp("file is empty, should be 0 byte file")
      else:
        sapp("file contains only 'from __future__ import' and should be empty instead")
    # elif self.n_from_future_import_division == 0:
    #   sapp("missing 'from __future__ import division'")
    elif self.n_from_future_import_division > 1:
      sapp("more than one appearance of 'from __future__ import division'")
    if self.flag_absolute_import and not self.file_should_be_empty:
      if self.n_from_future_import_absolute_import == 0:
        sapp("missing 'from __future__ import absolute_import'")
      elif self.n_from_future_import_absolute_import > 1:
        sapp("more than one appearance of 'from __future__ import absolute_import'")
    if self.flag_print_function and not self.file_should_be_empty:
      if self.n_from_future_import_print_function == 0:
        sapp("missing 'from __future__ import print_function'")
      elif self.n_from_future_import_print_function > 1:
        sapp("more than one appearance of 'from __future__ import print_function'")
    if (self.bad_indentation is not None) and (flag_indentation):
      n_tab, n_space = self.bad_indentation
      sapp("non-standard indentation: %d space, %d tab" % (n_space, n_tab))
    return ", ".join(status)

  def show(self, flag_x, flag_dos_format=True, append=None, verbose=False,
      flag_indentation=False):
    status = self.status(flag_x, flag_dos_format, flag_indentation)
    if status:
      msg = "%s: %s" % (self.path, status)
      if append:
        append(msg)
      else:
        print(msg)
      if (verbose) and (self.has_unused_imports()):
        msg2 = "  unused imports: %s" % ", ".join(self.unused_imports)
        if append:
          append(msg2)
        else :
          print(msg2)

def is_text_file(file_name):
  name = file_name.lower()
  for extension in (".c", ".cpp", ".h", ".hpp", ".py", ".java", ".params",
                    ".dox", ".txt", ".html", ".csh", ".sh", ".cif", ".cc",
                    '.mol2', '.frcmod',
                    ):
    if (name.endswith(extension)): return True
  return False

def gather(paths, find_unused_imports=False, find_bad_indentation=False, flag_absolute_import=False, flag_print_function=False):
  clutter = []
  def capp():
    clutter.append(file_clutter(path, find_unused_imports,
      find_bad_indentation=find_bad_indentation, flag_absolute_import=flag_absolute_import,
      flag_print_function=flag_print_function))
  for path in paths:
    if (not os.path.exists(path)):
      print("No such file or directory:", path, file=sys.stderr)
    elif (os.path.isfile(path)):
      capp()
    else:
      for root, dirs, files in os.walk(path):
        for f in sorted(files):
          if (is_text_file(f)):
            path = os.path.normpath(os.path.join(root, f))
            capp()
  return clutter

def detect_indentation_problems(file_name):
  try :
    import indent_finder
  except ImportError :
    return None
  fi = indent_finder.IndentFinder()
  fi.clear()
  fi.parse_file(file_name)
  result = fi.results()
  if (result is fi.default_result):
    return None
  itype, ival = result
  n_tab = n_space = 0
  if (itype != "mixed"):
    if (itype == "space"):
      n_space = ival
    else :
      n_tab = ival
    if (n_tab == 0) and (n_space == 2) : # this is our "standard"
      return None
    return (n_tab, n_space)
  else :
    return ival
  return None
