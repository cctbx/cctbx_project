"Minimalistic, pragmatic implementation. Absolutely no bells and whistles."

def run(args):
  assert args[0] == "--yes-i-know-this-overwrites-the-original-files"
  import fable.read
  import os
  op = os.path
  n_files_changed = 0
  for file_name in args[1:]:
    all_prcds = fable.read.process(file_names=[op.abspath(file_name)])
    insert_info = []
    for prcd in all_prcds.all_in_input_order:
      if (prcd.prcd_type == "blockdata"): continue
      if (len(prcd.executable) == 0):
        print "WARNING: no executable statements in %s" % prcd.name.value
      else:
        sl0 = prcd.executable[0].ssl.source_line_cluster[0]
        assert sl0.file_name == file_name
        insert_info.append((prcd.name.value, sl0.line_number))
    if (len(insert_info) != 0):
      insert_info.reverse()
      lines = open(file_name).read().splitlines()
      for name,line_number in insert_info:
        lines.insert(line_number-1,
          "      write(6, '(a)') 'PROCEDURE_START: %s'" % name)
      print >> open(file_name, "w"), "\n".join(lines)
      n_files_changed += 1
  print "Number of files changed:", n_files_changed

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
