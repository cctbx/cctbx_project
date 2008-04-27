import iotbx.pdb
from cctbx import uctbx
import sys

def run(args):
  for file_name in args:
    for pdb_str in open(file_name):
      sys.stdout.write(pdb_str)
      if (pdb_str.startswith("CRYST1")):
        cryst1_record = iotbx.pdb.records.cryst1(pdb_str=pdb_str)
        if (cryst1_record.ucparams is not None):
          print iotbx.pdb.format_scale_records(
            unit_cell=uctbx.unit_cell(cryst1_record.ucparams))

if (__name__ == "__main__"):
  run(sys.argv[1:])
