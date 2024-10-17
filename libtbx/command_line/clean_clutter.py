from __future__ import absolute_import, division, print_function
import sys
import os
import shutil
from libtbx.option_parser import option_parser
import libtbx.file_clutter

def clean_clutter_in(files, tabsize=8):
  if not files: return
  for fname in files:
    tmpname = fname + '.bak'
    if os.path.isfile(tmpname):
      print("found temporary file {temp}, ignoring {original}.".format(temp=tmpname, original=fname))
      continue
    if os.path.isfile(fname):
      try:
        print(fname)
        with open(fname, 'rb') as ifh, open(tmpname, 'wb') as ofh:
          # explicitly convert Windows linebreaks into Unix linebreaks
          lines = ifh.read().replace(b'\r\n', b'\n').split(b'\n')
          n_empty = 0
          for line in lines:
            clean_line = line.expandtabs(tabsize).rstrip()
            if clean_line:
              ofh.write(b"\n" * n_empty + clean_line + b"\n")
              n_empty = 0
            else:
              n_empty += 1
        shutil.move(tmpname, fname)
      except: # intentional
              # to trap KeyboardInterrupt, too
        os.remove(tmpname)
        raise

def isort(path):
  # Potential ImportErrors are caught upstream
  import mock
  from isort.main import main
  return # Disable isort pending resolution of https://github.com/timothycrosley/isort/issues/606
  with mock.patch.object(sys, 'argv', ['isort', '-y', '-ac', '-vb']):
    oldcwd = os.getcwd()
    try:
      os.chdir(path)
      main()
    finally:
      os.chdir(oldcwd)

def run():
  opt_parser = (option_parser(
    usage="""
clean_clutter [-t n | --tabsize=n] file1 file2 ...
clean_clutter [-t n | --tabsize=n] [directory]""",
    description="""The first form cleans the specified files whereas the second
form cleans all files in the hierarchy rooted in the given directory or
the current directory is none is given.""")
    .option("-t", "--tabsize",
      action="store",
      type="int",
      default=8,
      help="the number of spaces a tab is to be replaced by",
      metavar="INT")
  )
  command_line = opt_parser.process(args=sys.argv[1:])
  co = command_line.options
  files = command_line.args
  run_isort_in_path = False
  if len(files) <= 1:
    if not files: dir = '.'
    else: dir = files[0]
    files = [ c.path for c in libtbx.file_clutter.gather([dir])
              if c.is_cluttered(flag_x=False) ]
    if os.path.exists(os.path.join(dir, '.isort.cfg')):
      run_isort_in_path = dir
  clean_clutter_in(files, tabsize=co.tabsize)
  if run_isort_in_path:
    try:
      isort(run_isort_in_path)
    except Exception as e:
      print("Did not run isort (%s)" % str(e))

if __name__ == "__main__":
  run()
