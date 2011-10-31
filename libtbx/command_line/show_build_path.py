import libtbx.load_env

def run():
  print abs(libtbx.env.build_path)

if (__name__ == "__main__"):
  run()
