from iotbx.dtrek import reflnlist_reader

def extract_from(file_name):
  return reflnlist_reader.reflnlist(
    open(file_name), header_only=0001).crystal_symmetry()
