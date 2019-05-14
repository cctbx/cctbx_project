"""\
Simple script for post-processing output of TAPENADE at:

  http://tapenade.inria.fr:8080/tapenade/paste.jsp

TAPENADE is an Automatic Differentiation Engine developed at INRIA
Sophia-Antipolis by the TROPICS team.

http://www-sop.inria.fr/tropics/tapenade.html
"""
from __future__ import absolute_import, division, print_function

def run(args):
  file_names = []
  no_comments = False
  for arg in args:
    if (arg == "--no-comments"):
      no_comments = True
    else:
      file_names.append(arg)
  lines_to_ignore = set("""\
      INCLUDE 'DIFFSIZES.inc'
C  Hint: nbdirsmax should be the maximum number of differentiation directions
""".splitlines())
  default_replacement_dict = {"nbdirsmax": "nbdirs"}
  output_lines = []
  for file_name in file_names:
    replacement_dict = dict(default_replacement_dict)
    for line in open(file_name).read().splitlines():
      line = line.rstrip()
      if (line in lines_to_ignore):
        continue
      if (line.startswith("C  Hint: ")):
        flds = line.split()
        assert len(flds) == 9
        assert " ".join(flds[3:8]) == "should be the value of"
        replacement_dict[flds[2]] = flds[8]
        continue
      if (no_comments and line.startswith("C")):
        continue
      if (line == "      END"):
        replacement_dict = dict(default_replacement_dict)
      else:
        for old,new in replacement_dict.items():
          line = line.replace(old, new)
      output_lines.append(line)
  while (len(output_lines) != 0 and len(output_lines[-1]) == 0):
    output_lines.pop()
  print("\n".join(output_lines))

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
