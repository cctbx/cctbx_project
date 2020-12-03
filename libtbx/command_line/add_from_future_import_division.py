from __future__ import absolute_import, division, print_function
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

  with open(module_path) as fh:
    module_lines = fh.readlines()
  syntax_tree = ast.parse(''.join(module_lines), filename=module_path)
  # the attribute lineno is the index of the line *following* that node
  if ast.get_docstring(syntax_tree):
    insertion_lineno = syntax_tree.body[0].lineno
  else:
    insertion_lineno = 0

  for node in ast.iter_child_nodes(syntax_tree):
    if isinstance(node, ast.ImportFrom) and node.module == '__future__':
      insertion_lineno = node.lineno
      for n in node.names:
        want.discard(n.name)

  if not want:
    return

  missing_future = ", ".join(sorted(want))
  module_lines.insert(insertion_lineno, "from __future__ import %s\n" % missing_future)

  # Check syntax before making change
  try:
    modified_syntax_tree = ast.parse(''.join(module_lines))
  except SyntaxError:
    print("Cannot add {stmt} to {file}, causes syntax error".format(stmt=missing_future, file=module_path))
    return

  # Compare the modified syntax tree to the original tree
  for n, node in enumerate(modified_syntax_tree.body):
    if isinstance(node, ast.ImportFrom) and node.module == '__future__' and node.names[0].name in want:
      del(modified_syntax_tree.body[n])
      break
  if ast.dump(modified_syntax_tree) != ast.dump(syntax_tree):
    print("Cannot add {stmt} to {file}, changes parsing tree".format(stmt=missing_future, file=module_path))
    return

  print("Adding {stmt} to {file}".format(stmt=missing_future, file=module_path))
  with open(module_path, mode='w') as fh:
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
      fix(l, imports=check_for_imports)
    else:
      for dirpath, dirs, filenames in os.walk(l):
        dirs[:] = [d for d in dirs if not d.endswith('jinja2') and d not in ('.git', '__pycache__')]
        for f in filenames:
          if f.endswith('.py'):
            fix(os.path.join(dirpath, f), imports=check_for_imports)

if __name__ == '__main__':
  run(sys.argv[1:])
