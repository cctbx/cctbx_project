from iotbx.cns import sdb_reader

def extract_from(file_name=None, file=None):
  assert [file_name, file].count(None) == 1
  if (file is None):
    file = open(file_name)
  try: sdb_files = sdb_reader.multi_sdb_parser(file)
  except: return None
  assert len(sdb_files) > 0
  return sdb_files[0].crystal_symmetry()
