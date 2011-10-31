import libtbx.load_env

def run():
  print ":".join([ abs(p) for p in libtbx.env.pythonpath ])

if (__name__ == "__main__"):
  run()
