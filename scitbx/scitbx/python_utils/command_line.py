class parse_options:

  def __init__(self, argv, keywords):
    self.keywords = keywords
    self.n = 0
    for keyword in keywords:
      setattr(self, keyword, False)
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
        setattr(self, flds[0], True)
      else:
        setattr(self, flds[0], flds[1])
      self.n += 1

def parse_options_with_chunk(argv, keywords=[]):
  flags = parse_options(argv, ["ChunkSize", "ChunkMember"] + keywords)
  if (flags.ChunkSize == False):
    flags.ChunkSize = 1
  else:
     flags.ChunkSize = int(flags.ChunkSize)
  if (flags.ChunkMember == False):
    flags.ChunkMember = 0
  else:
     flags.ChunkMember = int(flags.ChunkMember)
  assert flags.ChunkSize > 0
  assert flags.ChunkMember >= 0
  assert flags.ChunkMember < flags.ChunkSize
  return flags
