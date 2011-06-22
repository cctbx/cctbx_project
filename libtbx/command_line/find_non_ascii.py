import os, sys, string

def visit(arg, dirname, names):
  if not os.path.basename(dirname)=='CVS':
    for name in names:
      if not name.endswith('.py'):
        continue
      filename = os.path.join(dirname,name)
      dot_filename = filename[len(working_dir)+1:]
      try:
        f = file(filename, 'r')
        lines = f.readlines()
        f.close()
      except IOError, e:
        continue

      for j, line in enumerate(lines):
        for i, l in enumerate(line):
          if not l in string.printable:
            print '-'*80
            print "Non-ASCII letter at position %d in line %d of file\n\t%s" \
                  % (i, j, dot_filename)
            print line
            print "%s^" % (' '*i)
            break

def run(dir):
  global working_dir
  working_dir = dir
  os.path.walk(dir, visit, 'show')

if (__name__ == "__main__"):
  try:
    run(sys.argv[1])
  except Exception:
    run(os.getcwd())
