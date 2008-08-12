from libtbx.utils import detect_binary_file
from libtbx.path import walk_source_tree
from libtbx.option_parser import option_parser
from libtbx.str_utils import show_string
import sys, os

def show_status(path, text, binary, quote):
  def show():
    if (quote): print show_string(path)
    else: print path
  if (text and binary):
    show()
  else:
    status = detect_binary_file.from_initial_block(file_name=path)
    if (status is None or status is binary):
      show()

def run(args, command_name="libtbx.list_files"):
  if (len(args) == 0): args = ["."]
  command_line = (option_parser(
    usage="%s [options] path ..." % command_name,
    description="Recursively lists all files,"
      " excluding CVS and .svn directories and .pyc files.")
    .option("-t", "--text",
      action="store_true",
      default=False,
      help="list text files only")
    .option("-b", "--binary",
      action="store_true",
      default=False,
      help="list binary files only")
    .option("-q", "--quote",
      action="store_true",
      default=False,
      help="quote file names")
  ).process(args=args)
  paths = command_line.args
  co = command_line.options
  text = co.text
  binary = co.binary
  quote = co.quote
  if (not (text or binary)):
    binary = True
    text = True
  if (len(paths) == 0): paths = ["."]
  for path in paths:
    if (not os.path.exists(path)):
      print >> sys.stderr, "No such file or directory:", path
    elif (os.path.isfile(path)):
      show_status(path=path, text=text, binary=binary, quote=quote)
    else:
      for file_path in walk_source_tree(top=path):
        show_status(path=file_path, text=text, binary=binary, quote=quote)

if (__name__ == "__main__"):
  run(sys.argv[1:])
