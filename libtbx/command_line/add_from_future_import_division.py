from __future__ import division
import ast
import os
import sys
import libtbx.load_env

def fix(module_path, imports=None):
  if not imports:
    imports = ['division']
  want = set(imports)

  if os.path.getsize(module_path) == 0: # Do not add __future__ imports to empty files
    return

  with file(module_path) as fh:
    module_lines = fh.readlines()
  module_node = ast.parse(''.join(module_lines), filename=module_path)
  # the attribute lineno is the index of the line *following* that node
  if ast.get_docstring(module_node):
    insertion_lineno = module_node.body[0].lineno
  else:
    insertion_lineno = 0

  for node in ast.iter_child_nodes(module_node):
    if isinstance(node, ast.ImportFrom) and node.module == '__future__':
      insertion_lineno = node.lineno
      module = node.module.split('.')
      for n in node.names:
        if n.name in want:
          want.remove(n.name)

  if not want:
    return

  missing_future = ", ".join(sorted(want))
  print("Adding {stmt} to {file}".format(stmt=missing_future, file=module_path))
  module_lines.insert(insertion_lineno, "from __future__ import %s\n" % missing_future)
  with file(module_path, mode='w') as fh:
    fh.writelines(module_lines)

def run(locations):
  check_for_imports = ['division']
  if '--absolute_import' in locations:
    check_for_imports.append('absolute_import')
    locations.remove('--absolute_import')
  if '--print_function' in locations:
    check_for_imports.append('print_function')
    locations.remove('--print_function')
  if '-3' in locations:
    check_for_imports.append('absolute_import')
    check_for_imports.append('print_function')
    locations.remove('-3')

  if not locations:
    locations = [ os.path.dirname(libtbx.env.dist_path('libtbx')) ]
  for l in locations:
    if os.path.isfile(l) and l.endswith('.py'):
      print(l)
      fix(l, imports=check_for_imports)
    else:
      for dirpath, dirnames, filenames in os.walk(l):
        if dirpath.endswith('jinja2'): continue
        print(dirpath)
        for f in filenames:
          if f.endswith('.py'):
            fix(os.path.join(dirpath, f), imports=check_for_imports)

if __name__ == '__main__':
  run(sys.argv[1:])
