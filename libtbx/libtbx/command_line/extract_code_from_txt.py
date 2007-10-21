import sys

def run(args, command_name="libtbx.extract_code_from_txt"):
  markup_extract_begin = "  ## extract code begin: "
  markup_extract_end = "  ## extract code end"
  outputs = {}
  for file_name in args:
    lines = open(file_name).read().splitlines()
    for i_line,line in enumerate(lines):
      if (line.startswith(markup_extract_begin)):
        destination = line[len(markup_extract_begin):].strip()
        buffer = [("  # ---- line %d " % (i_line+1) + "-"*75)[:79]]
        for line in lines[i_line+1:]:
          if (line.startswith(markup_extract_begin)):
            raise RuntimeError(
              'Unexpected markup: "%s" (line %d)' % (
                markup_extract_begin, i_line+1))
          if (line.rstrip() == markup_extract_end):
            break
          buffer.append(line)
        else:
          raise RuntimeError(
            'Unexpected end of file: missing "%s"' % markup_extract_end)
        outputs.setdefault(destination, []).extend(buffer)
  for destination,lines in outputs.items():
    lines.insert(0, 'if (__name__ == "__main__"):')
    lines.insert(1, "")
    print "Writing: %s (%d lines)" % (destination, len(lines))
    print >> open(destination, "w"), "\n".join(lines)

if (__name__ == "__main__"):
  run(sys.argv[1:])
