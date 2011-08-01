import os, fnmatch
from libtbx import python_code_parsing

def run(args):
  if not args: args = [ '.' ]
  work = set()
  arg_filenames = []
  for arg in args:
    if os.path.isdir(arg):
      for dirpath, dirnames, filenames in os.walk(arg):
        work.update( os.path.join(dirpath, f)
                     for f in fnmatch.filter(filenames, '*.py') )
    else:
      arg_filenames.append(arg)
  work.update(fnmatch.filter(arg_filenames, '*.py'))
  for filename in work:
    unused = python_code_parsing.unused_imports(
      python_source_filename=filename,
      ignored_imports=('libtbx.load_env',
                       'libtbx.forward_compatibility',
                       'import libtbx.start_print_trace',
                       'import libtbx.callbacks'),
      ignored_imports_from=('__future__',),
      ignore_imports_flagged_by_comments=('# import dependency',
                                          '# implicit import'))
    if unused:
      print 'In file %s:' % filename
      print unused
      print


if __name__ == '__main__':
  import sys
  run(args=sys.argv[1:])
