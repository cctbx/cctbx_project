import os
op = os.path

def inspect(py_lines):
  imports_to_ignore = set([
    "import libtbx.forward_compatibility",
    "  import libtbx.start_print_trace",
    "    import libtbx.callbacks"])
  combined_lines = []
  block = []
  for line in py_lines:
    if (line in imports_to_ignore):
      continue
    l = line.strip()
    if (not l.endswith("\\")):
      block.append(l)
      combined_lines.append(" ".join(block))
      block = []
    else:
      block.append(l[:-1])
  if (len(block) != 0):
    combined_lines.append("".join(block))
  imported_names_dict = {}
  non_import_lines = []
  for l in combined_lines:
    def split():
      return l.replace(",", " , ").split()
    if (l.startswith("import ")):
      if (l.endswith(" # import dependency")):
        continue
      if (l.endswith(" # implicit import")):
        continue
      if (l.endswith(" # special import")):
        continue
      flds = split()
    elif (l.startswith("from ")):
      if (l.startswith("from __future__ ")):
        continue
      flds = split()
      for i,fld in enumerate(flds):
        if (fld == "import"):
          flds = flds[i:]
          break
      else:
        continue
    else:
      non_import_lines.append(l)
      continue
    assert flds[0] == "import", flds
    flds = flds[1:]
    flds.append(",")
    #
    def collect(flds):
      if (len(flds) == 1):
        name = flds[0]
      elif (len(flds) == 3):
        name = flds[2]
      else:
        return
      if (name == "libtbx.load_env"):
        name = "env"
      else:
        name = name.split(".")[-1] # XXX very crude
      if (name != "*" and name not in imported_names_dict):
        imported_names_dict[name] = len(imported_names_dict)
    i = 0
    for j,fld in enumerate(flds):
      if (fld == ","):
        if (j > i):
          collect(flds[i:j])
        i = j+1
  idc = "_0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
  imported_names = set(imported_names_dict.keys())
  used_names = set()
  for l in non_import_lines:
    filtered = []
    for c in l:
      if (idc.find(c) < 0):
        c = " "
      filtered.append(c)
    flds = "".join(filtered).split()
    used_names.update(imported_names.intersection(flds))
  unused_names = list(imported_names - used_names)
  def cmp_input_order(a, b):
    return cmp(*(imported_names_dict[a], imported_names_dict[b]))
  unused_names.sort(cmp_input_order)
  return unused_names

def show_unused_imports(file_name):
  unused_imports = inspect(
    py_lines=open(file_name).read().splitlines())
  if (len(unused_imports) != 0):
    print "%s: %s" % (file_name, ", ".join(unused_imports))
    print
  return unused_imports

def walk_func(counter, dirname, names):
  for name in names:
    if (not name.endswith(".py")):
      continue
    file_name = op.join(dirname, name)
    if (op.isfile(file_name)):
      if (len(show_unused_imports(file_name)) != 0):
        counter[0] += 1

def run(args):
  if (len(args) == 0):
    args = ["."]
  counter = [0]
  for arg in args:
    if (op.isdir(arg)):
      op.walk(arg, walk_func, counter)
    elif (op.isfile(arg)):
      if (len(show_unused_imports(file_name=arg)) != 0):
        counter[0] += 1
  if (counter[0] != 0):
    print """\
HINT:
  To suppress flagging of unused imports follow these examples:
    import scitbx.array_family.flex # import dependency
    import something.related # implicit import
    import wingdbstub # special import
"""

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
