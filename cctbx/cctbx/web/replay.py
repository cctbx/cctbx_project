import sys, pickle

class empty: pass

def run(file_names):
  for file_name in file_names:
    print "Replaying:", file_name
    f = open(file_name, "rb")
    server_info, target_module, inp = pickle.load(f)
    f.close()
    exec "import " + target_module + " as target"
    status = empty()
    status.in_table = False
    target.run(server_info, inp, status)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
