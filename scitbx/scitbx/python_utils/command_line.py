class parse_options:

  def __init__(self, argv, keywords):
    self.keywords = keywords
    self.n = 0
    for keyword in keywords:
      setattr(self, keyword, 0)
    self.regular_args = []
    for arg in argv:
      if (not arg.startswith("--")):
        self.regular_args.append(arg)
        continue
      if (not arg[2:] in keywords):
        raise RuntimeError, "Unknown option: " + arg
      setattr(self, arg[2:], 1)
      self.n += 1
