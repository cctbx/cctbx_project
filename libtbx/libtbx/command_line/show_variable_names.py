import libtbx.load_env
pairs = dict(libtbx.env.var_name_and_build_or_dist_path_pairs())
var_names = pairs.keys()
var_names.sort()
for var_name in var_names:
  print "%s=%s" % (var_name, pairs[var_name])
