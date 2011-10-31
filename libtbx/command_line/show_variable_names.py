import libtbx.load_env

def run():
  pairs = dict(libtbx.env.var_name_and_build_or_dist_path_pairs())
  var_names = pairs.keys()
  var_names.sort()
  for var_name in var_names:
    print "%s=%s" % (var_name, abs(pairs[var_name]))

if (__name__ == "__main__"):
  run()
