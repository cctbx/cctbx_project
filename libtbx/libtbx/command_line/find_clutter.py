import sys, os

def show_status(path):
  stream = open(path, "rb").read()
  if (len(stream) == 0): return
  dos_format = 00000
  n_tabs_or_trailing_whitespace = 0
  n_trailing_empty_lines = 0
  missing_eol = 00000
  if (stream[-1] != "\n"):
    missing_eol = 0001
  else:
    stream = stream[:-1]
  text = stream.split("\n")
  for line in text:
    if (line.endswith("\r")):
      line = line[:-1]
      dos_format = 0001
    clean_line = line.expandtabs().rstrip()
    if (clean_line != line): n_tabs_or_trailing_whitespace += 1
    if (len(clean_line) == 0): n_trailing_empty_lines += 1
    else: n_trailing_empty_lines = 0
  status = ""
  if (dos_format):
    status += "dos format, "
  if (n_tabs_or_trailing_whitespace > 0):
    status += "tabs or trailing whitespace=%d, " \
           % n_tabs_or_trailing_whitespace
  if (n_trailing_empty_lines > 0):
    status += "trailing empty lines=%d, " % n_trailing_empty_lines
  if (missing_eol):
    status += "missing end-of-line, "
  if (len(status) > 0):
    print "%s:" % path, status[:-2]

def is_text_file(file_name):
  name = file_name.lower()
  if (name.startswith("scons")): return 0001
  for extension in (".c", ".cpp", ".h", ".hpp", ".py",
                    ".dox", ".txt", ".html", ".csh", ".sh"):
    if (name.endswith(extension)): return 0001
  return 00000

def visitor(arg, dirname, names):
  for file_name in names:
    if (is_text_file(file_name)):
      path = os.path.normpath(os.path.join(dirname, file_name))
      show_status(path)

def run(paths):
  if (len(paths) == 0): paths = ["."]
  for path in paths:
    if (os.path.isfile(path)):
      show_status(path)
    else:
      os.path.walk(path, visitor, None)

if (__name__ == "__main__"):
  run(sys.argv[1:])
