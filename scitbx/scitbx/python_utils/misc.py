def adopt_init_args(obj, args):
  del args["self"]
  obj.__dict__.update(args)

def import_regular_symbols(dict_target, dict_source):
  for key, value in dict_source.items():
    if (not key.startswith("_") and not key in dict_target):
      dict_target[key] = value
