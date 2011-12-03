def run(args):
  import libtbx.load_env
  remaining_args = []
  dirname_count = 0
  for arg in args:
    if (arg == "--dirname"):
      dirname_count += 1
    else:
      remaining_args.append(arg)
  def show(path):
    if (path is not None):
      from os.path import dirname
      for _ in xrange(dirname_count):
        path = dirname(path)
    print path
  if (len(remaining_args) == 0):
    for path in libtbx.env.dist_paths():
      show(path)
  else:
    for arg in remaining_args:
      show(libtbx.env.dist_path(module_name=arg, default=None))

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
