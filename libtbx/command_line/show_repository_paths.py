import libtbx.load_env

def run():
  for path in libtbx.env.repository_paths:
    print abs(path)

if (__name__ == "__main__"):
  run()
