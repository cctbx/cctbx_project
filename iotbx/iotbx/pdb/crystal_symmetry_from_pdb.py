from cctbx import crystal
from iotbx.pdb import parser
from iotbx.misc import detect_binary_file

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
    if (not line.startswith("CRYST1")): continue
    cryst1 = parser.pdb_record(raw_record=line, line_number=line_number)
    assert [cryst1.ucparams, cryst1.sGroup].count(None) < 2
    return crystal.symmetry(
      unit_cell=cryst1.ucparams,
      space_group_symbol=cryst1.sGroup)
  raise RuntimeError, "No CRYST1 record."
