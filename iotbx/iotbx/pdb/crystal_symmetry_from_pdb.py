from iotbx.pdb import cryst1_interpretation
from iotbx.misc import detect_binary_file
from iotbx.cns import pdb_remarks as cns_pdb_remarks

def extract_from(file_name=None, file=None, monitor_initial=None):
  assert [file_name, file].count(None) == 1
  if (file is None):
    file = open(file_name)
  detect_binary = detect_binary_file(monitor_initial=monitor_initial)
  line_number = 0
  for line in file:
    line_number += 1
    if (detect_binary is not None):
      is_binary = detect_binary.is_binary_file(line)
      if (is_binary is not None):
        if (is_binary): break
        detect_binary = None
    if (line.startswith("CRYST1")):
      return cryst1_interpretation.crystal_symmetry(
        cryst1_record=line,
        line_number=line_number)
    crystal_symmetry = cns_pdb_remarks.extract_symmetry(
      pdb_record=line)
    if (crystal_symmetry is not None):
      return crystal_symmetry
  raise RuntimeError, "No CRYST1 record."
