from __future__ import absolute_import, division, print_function
from itertools import count
import sys, os
from six.moves import zip

def run(args, command_name="libtbx.extract_code_from_txt"):
  markup_extract_begin = "  ## extract code begin: "
  markup_extract_end = "  ## extract code end"
  outputs = {}
  for file_name in args:
    rst = []
    i_line_lines = iter(zip(count(), open(file_name).read().splitlines()))
    for i_line,line in i_line_lines:
      if (not line.startswith(markup_extract_begin)):
        rst.append(line)
      else:
        destination = line[len(markup_extract_begin):].strip()
        buffer = [("  # ---- line %d " % (i_line+1) + "-"*75)[:79]]
        def check_next_is_empty_line():
          i_line,line = next(i_line_lines)
          if (len(line) != 0):
            raise RuntimeError(
              "Markup must be followed by an empty line (line %d)" %
                (i_line+1))
        check_next_is_empty_line()
        buffer.append("")
        for i_line,line in i_line_lines:
          if (line.startswith(markup_extract_begin)):
            raise RuntimeError(
              'Unexpected markup: "%s" (line %d)' % (
                markup_extract_begin, i_line+1))
          if (line.rstrip() == markup_extract_end):
            check_next_is_empty_line()
            break
          buffer.append(line)
          rst.append(line)
        else:
          raise RuntimeError(
            'Unexpected end of file: missing "%s"' % markup_extract_end)
        outputs.setdefault(destination, []).extend(buffer)
    file_rst = os.path.basename(file_name)
    if (file_rst.endswith(".txt")):
      file_rst = file_rst[:-4]
    file_rst += ".rst"
    print("Writing: %s (%d lines)" % (file_rst, len(rst)))
    print("\n".join(rst), file=open(file_rst, "w"))
  for destination,lines in outputs.items():
    lines.insert(0, 'if (__name__ == "__main__"):')
    lines.insert(1, "")
    print("Writing: %s (%d lines)" % (destination, len(lines)))
    print("\n".join(lines), file=open(destination, "w"))

if (__name__ == "__main__"):
  run(sys.argv[1:])
