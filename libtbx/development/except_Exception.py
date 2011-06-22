def run(args):
  if (len(args) == 0):
    folders = ["."]
  else:
    folders = args
  from libtbx.path import walk_source_tree
  mod_count_total = 0
  mod_file_count = 0
  for folder in folders:
    for path in walk_source_tree(folder):
      if (not path.endswith(".py")): continue
      txt = open(path).read()
      mod_lines = []
      mod_count = 0
      for line in txt.splitlines():
        ls = line.strip()
        if (    ls.startswith("except")
            and ls[6:].strip().startswith(":")
            and not ls.endswith(" # intentional")):
          line = line.replace("except", "except Exception", 1)
          mod_count += 1
        mod_lines.append(line)
      if (mod_count != 0):
        print >> open(path, "w"), "\n".join(mod_lines)
        mod_count_total += mod_count
        mod_file_count += 1
  print "Number of modifications: %d in %d files" % (
    mod_count_total, mod_file_count)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
