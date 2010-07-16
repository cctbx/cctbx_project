def run(args):
  import fable.read
  out_names_used = set()
  for file_name in args:
    all_units = fable.read.process(
      file_names=[file_name],
      basic_only=True,
      skip_load_includes=True)
    for unit in all_units.all_in_input_order:
      out_name = unit.name.value
      i = 2
      while (out_name in out_names_used):
        out_name = "%s_%d" % (unit.name.value, i)
        i += 1
      out_names_used.add(out_name)
      out = open(out_name+".f", "w")
      print out.name
      first_line = True
      empty_lines = []
      for ssl in unit.all_ssl():
        for sl in ssl.source_line_cluster:
          line = sl.text
          if (len(line.strip()) == 0):
            empty_lines.append(line)
          else:
            if (not first_line):
              for prev_line in empty_lines:
                print >> out, prev_line
            print >> out, line
            first_line = False
            empty_lines = []
      del out

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
