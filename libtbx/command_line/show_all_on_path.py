import os

def run():
  done = set()
  for directory in os.environ["PATH"].split(os.pathsep):
    directory = os.path.normcase(os.path.abspath(directory))
    if (directory in done): continue
    done.add(directory)
    if (not os.path.isdir(directory)): continue
    if (not os.access(directory, os.X_OK)): continue
    for file_name in os.listdir(directory):
      if (file_name.lower().startswith("libtbx.find_clutter")):
        print directory
        break

if (__name__ == "__main__"):
  run()
