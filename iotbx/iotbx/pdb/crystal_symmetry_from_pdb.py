from cctbx import crystal
from iotbx.misc import detect_binary_file

def extract_from(file_name=None, file=None, monitor_initial=None):
  assert [file_name, file].count(None) == 1
  if (file is None):
    file = open(file_name)
  detect_binary = detect_binary_file(monitor_initial=monitor_initial)
  for line in file:
    if (detect_binary is not None):
      is_binary = detect_binary.is_binary_file(line)
      if (is_binary is not None):
        if (is_binary): break
        detect_binary = None
    if (not line.startswith("CRYST1")): continue
    #  7 - 15       Real(9.3)      a             a (Angstroms).
    # 16 - 24       Real(9.3)      b             b (Angstroms).
    # 25 - 33       Real(9.3)      c             c (Angstroms).
    # 34 - 40       Real(7.2)      alpha         alpha (degrees).
    # 41 - 47       Real(7.2)      beta          beta (degrees).
    # 48 - 54       Real(7.2)      gamma         gamma (degrees).
    # 56 - 66       LString        sGroup        Space group.
    # 67 - 70       Integer        z             Z value.
    raw_record = (line + " " * 80)[:80]
    ucparams = " ".join(
      [raw_record[ 6:15], raw_record[15:24], raw_record[24:33],
      raw_record[33:40], raw_record[40:47], raw_record[47:54]]).strip()
    sGroup = raw_record[55:66].strip()
    if (len(ucparams) == ""): ucparams = None
    if (len(sGroup) == ""): sGroup = None
    assert [ucparams, sGroup].count(None) < 2
    return crystal.symmetry(
      unit_cell=ucparams,
      space_group_symbol=sGroup)
  raise RuntimeError, "No CRYST1 record."
