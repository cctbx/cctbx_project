import sys, os

class file_clutter(object):

  def __init__(self, path):
    self.path = path
    self.is_executable = os.access(path, os.X_OK)
    self.dos_format = False
    self.n_tabs_or_trailing_whitespace = 0
    self.n_trailing_empty_lines = 0
    self.missing_eol = False
    stream = open(path, "rb").read()
    if (len(stream) > 0):
      if (stream[-1] != "\n"):
        self.missing_eol = True
      else:
        stream = stream[:-1]
      text = stream.split("\n")
      for line in text:
        if (line.endswith("\r")):
          line = line[:-1]
          self.dos_format = True
        clean_line = line.expandtabs().rstrip()
        if (clean_line != line): self.n_tabs_or_trailing_whitespace += 1
        if (len(clean_line) == 0): self.n_trailing_empty_lines += 1
        else: self.n_trailing_empty_lines = 0

  def is_cluttered(self, flag_x):
    return ((self.is_executable and flag_x)
            or self.dos_format
            or self.n_tabs_or_trailing_whitespace > 0
            or self.n_trailing_empty_lines > 1
            or self.missing_eol)

  def status(self, flag_x, flag_dos_format=True):
    status = ""
    if (self.is_executable and flag_x
        and self.path.lower().find("command_line") < 0
        and not self.path.endswith(".csh")
        and not self.path.endswith(".sh")):
      status += "is executable, "
    if (flag_dos_format and self.dos_format):
      status += "dos format, "
    if (self.n_tabs_or_trailing_whitespace > 0):
      status += "tabs or trailing whitespace=%d, " \
             % self.n_tabs_or_trailing_whitespace
    if (self.n_trailing_empty_lines > 1):
      status += "trailing empty lines=%d, " % self.n_trailing_empty_lines
    if (self.missing_eol):
      status += "missing end-of-line, "
    return status

  def show(self, flag_x, flag_dos_format=True):
    status = self.status(flag_x, flag_dos_format)
    if status: print "%s: %s" % (self.path, status)


def is_text_file(file_name):
  name = file_name.lower()
  for extension in (".c", ".cpp", ".h", ".hpp", ".py", ".java", ".params",
                    ".dox", ".txt", ".html", ".csh", ".sh", ".cif"):
    if (name.endswith(extension)): return True
  return False

def gather(paths):
  clutter = []
  for path in paths:
    if (not os.path.exists(path)):
      print >> sys.stderr, "No such file or directory:", path
    elif (os.path.isfile(path)):
      clutter.append(file_clutter(path))
    else:
      for root, dirs, files in os.walk(path):
        for f in files:
          if (is_text_file(f)):
            path = os.path.normpath(os.path.join(root, f))
            clutter.append(file_clutter(path))
  return clutter
