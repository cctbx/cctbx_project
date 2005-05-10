from iotbx import pdb
import iotbx.pdb.parser
from cctbx import uctbx
import sys

def run(args):
  for file_name in args:
    for raw_record in open(file_name):
      sys.stdout.write(raw_record)
      if (raw_record.startswith("CRYST1")):
        cryst1_record = pdb.parser.pdb_record(raw_record=raw_record)
        if (cryst1_record.ucparams is not None):
          print pdb.format_scale_records(
            unit_cell=uctbx.unit_cell(cryst1_record.ucparams))

if (__name__ == "__main__"):
  run(sys.argv[1:])
