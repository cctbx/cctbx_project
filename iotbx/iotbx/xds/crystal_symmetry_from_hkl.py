from iotbx.xds import read_ascii

def extract_from(file_name):
  return read_ascii.reader(
    open(file_name), header_only=0001).crystal_symmetry()
