class parse_options:

  def __init__(self, argv, keywords):
    self.keywords = keywords
    self.n = 0
    for keyword in keywords:
      setattr(self, keyword, 00000)
    self.regular_args = []
    for arg in argv:
      if (not arg.startswith("--")):
        self.regular_args.append(arg)
        continue
      flds = arg[2:].split("=")
      assert len(flds) in (1,2)
      if (not flds[0] in keywords):
        raise RuntimeError, "Unknown option: " + arg
      if (len(flds) == 1):
        setattr(self, flds[0], 0001)
      else:
        setattr(self, flds[0], flds[1])
      self.n += 1
