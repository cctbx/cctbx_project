"Minimalistic, pragmatic implementation. Absolutely no bells and whistles."

def run(args):
  assert args[0] == "--yes-i-know-this-overwrites-the-original-files"
  import fable.read
  import os
  op = os.path
  n_files_changed = 0
  for file_name in args[1:]:
    file_name = op.abspath(file_name)
    all_fprocs = fable.read.process(file_names=[file_name])
    insert_info = []
    for fproc in all_fprocs.all_in_input_order:
      if (fproc.fproc_type == "blockdata"): continue
      if (len(fproc.executable) == 0):
        print "WARNING: no executable statements in %s" % fproc.name.value
      else:
        sl0 = fproc.executable[0].ssl.source_line_cluster[0]
        assert sl0.file_name == file_name
        insert_info.append((fproc.name.value, sl0.line_number))
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
