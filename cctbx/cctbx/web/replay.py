import sys, pickle

def run(file_names):
  for file_name in file_names:
    print "Replaying:", file_name
    f = open(file_name, "rb")
    cctbx_url, target_module, inp = pickle.load(f)
    f.close()
    exec "import " + target_module + " as target"
    target.run(cctbx_url, inp)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
