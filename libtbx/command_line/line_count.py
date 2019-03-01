from __future__ import absolute_import, division, print_function
import sys, os
import re
import libtbx.load_env

boost_python_include_pat = re.compile(r"#include\s*<boost(?:/|_)python");

def run(modules):
  directory_paths = [ libtbx.env.dist_path(m) for m in modules ]
  line_counts_in_files_of_type = {}
  for d in directory_paths:
    for root, dirs, files in os.walk(d):
      for f in files:
        if f.startswith('.'): continue
        _, ext = os.path.splitext(f)
        if ext in ('.pyo', '.pyc'): continue
        boost_python_binding = False
        n_lines = 0
        with open(os.path.join(root,f)) as fo:
          for li in fo:
            n_lines += 1
            if (not boost_python_binding
                and boost_python_include_pat.search(li)):
              boost_python_binding = True
        if boost_python_binding:
          file_type = "Boost.Python"
        elif not ext:
          file_type = "unknown"
        else:
          file_type = ext[1:]
        line_counts_in_files_of_type.setdefault(file_type, []).append(n_lines)
  print("Lines of code in %s" % ', '.join(modules))
  print("%-15s%8s" % ('extension', '#lines'))
  output = []
  for file_type, line_counts in line_counts_in_files_of_type.items():
    cnt = sum(line_counts)
    output.append((cnt, "%-15s%8d" % (file_type, cnt)))
  output.sort(reverse=True)
  output = [ entry[1] for entry in output ]
  print('\n'.join(output))


if __name__ == '__main__':
  run(sys.argv[1:])
