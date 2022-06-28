from __future__ import absolute_import, division, print_function
import sys, pickle
import importlib

class empty: pass

def run(file_names):
  for file_name in file_names:
    print("Replaying:", file_name)
    f = open(file_name, "rb")
    server_info, target_module, inp = pickle.load(f)
    f.close()
    target = importlib.import_module(target_module)
    status = empty()
    status.in_table = False
    target.run(server_info, inp, status)
    print()

if (__name__ == "__main__"):
  run(sys.argv[1:])
