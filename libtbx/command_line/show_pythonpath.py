import libtbx.load_env

def run():
  print ":".join(libtbx.env.pythonpath)

if (__name__ == "__main__"):
  run()
