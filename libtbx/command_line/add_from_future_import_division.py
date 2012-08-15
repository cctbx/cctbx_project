from __future__ import division
import libtbx.load_env
import os, sys
import ast

def fix(module_path):
  module_lines = file(module_path).readlines()
  module_node = ast.parse(''.join(module_lines), filename=module_path)
  # the attribute lineno is the index of the line *following* that node
  if ast.get_docstring(module_node):
    insertion_lineno = module_node.body[0].lineno
    index_to_check = 1
  else:
    insertion_lineno = 0
    index_to_check = 0
  if index_to_check < len(module_node.body):
    node_to_check = module_node.body[index_to_check]
    if isinstance(node_to_check, ast.ImportFrom) and node_to_check.module == '__future__':
      for name in node_to_check.names:
        if name.name == 'division':
          return
  module_lines.insert(insertion_lineno, "from __future__ import division\n")
  fixed_source = ''.join(module_lines)
  file(module_path, mode='w').writelines(module_lines)

def run(locations):
  if not locations:
    locations = [ os.path.dirname(libtbx.env.dist_path('libtbx')) ]
  modules = []
  for l in locations:
    if os.path.isfile(l) and l.endswith('.py'):
      fix(l)
    else:
      for dirpath, dirnames, filenames in os.walk(l):
        if dirpath.endswith('jinja2'): continue
        for f in filenames:
          if f.endswith('.py'): fix(os.path.join(dirpath, f))



if __name__ == '__main__':
  run(sys.argv[1:])
