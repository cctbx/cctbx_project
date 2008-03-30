import iotbx.pdb
from libtbx.utils import Usage, format_cpu_times
import sys, os

def run(args,  command_name="iotbx.pdb.join_fragment_files"):
  def usage():
    raise Usage("""\
%s file1.pdb file2.pdb [...]

or define the environment variable
  PDB_MIRROR_PDB
to join all fragment files in the PDB.""" % command_name)
  if (len(args) == 0 or args == ["--exercise"]):
    pdb_mirror_pdb = os.environ.get("PDB_MIRROR_PDB")
    if (pdb_mirror_pdb is None):
      if (len(args) == 0): usage()
    else:
      for line in iotbx.pdb.pdb_codes_fragment_files.splitlines():
        print "PDB code group:", line
        codes = line.split()
        joined = iotbx.pdb.join_fragment_files(
          file_names = [
            os.path.join(pdb_mirror_pdb, code[1:3], "pdb%s.ent.gz" % code)
              for code in codes])
        file_name_out = "%s_%s.pdb" % (codes[0], codes[-1])
        print "  writing:", file_name_out
        out = open(file_name_out, "w")
        print >> out, "\n".join(joined.info)
        out.write(joined.as_pdb_string(append_end=True))
        if (len(args) != 0): break
    print format_cpu_times()
  else:
    if (len(args) < 2): usage()
    joined = iotbx.pdb.join_fragment_files(file_names=args)
    print "\n".join(joined.info)
    sys.stdout.write(joined.as_pdb_string(append_end=True))

if (__name__ == "__main__"):
  run(sys.argv[1:])
