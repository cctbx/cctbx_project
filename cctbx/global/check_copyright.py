from os.path import walk

def visitfunc(arg, dirname, names):
  if (dirname.endswith("CVS")): return
  if (not "COPYRIGHT.txt" in names):
    print dirname

def run():
  walk(".", visitfunc, None)

if (__name__ == "__main__"):
  run()
