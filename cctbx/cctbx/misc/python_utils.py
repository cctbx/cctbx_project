def adopt_init_args(obj, args):
  del args["self"]
  obj.__dict__.update(args)
