import os, pickle

class env:

  def __init__(self):
    libtbx_build = os.environ["LIBTBX_BUILD"]
    file_name = os.path.join(libtbx_build, "libtbx_env")
    try:
      f = open(file_name)
    except:
      raise RuntimeError("Cannot open libtbx environment file: " + file_name)
    try:
      e = pickle.load(f)
    except:
      raise RuntimeError("Corrupt libtbx environment file: " + file_name)
    self.__dict__.update(e)
