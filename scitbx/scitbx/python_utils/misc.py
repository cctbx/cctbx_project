def adopt_init_args(obj, args, exclude=()):
  del args["self"]
  for param in exclude:
    del args[param]
  for key in args.keys():
    assert not hasattr(obj.__dict__, key)
  obj.__dict__.update(args)

def import_regular_symbols(dict_target, dict_source):
  for key, value in dict_source.items():
    if (not key.startswith("_") and not key in dict_target):
      dict_target[key] = value
