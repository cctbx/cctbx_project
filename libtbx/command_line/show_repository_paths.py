import libtbx.load_env

def run():
  for path in libtbx.env.repository_paths:
    print path

if (__name__ == "__main__"):
  run()
