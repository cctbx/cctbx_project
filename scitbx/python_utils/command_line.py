from __future__ import absolute_import, division, print_function
class parse_options(object):

  def __init__(self, argv, keywords, case_sensitive=True):
    self.keywords = keywords
    self.n = 0
    keywords_lower = []
    for keyword in keywords:
      setattr(self, keyword, False)
      if (not case_sensitive):
        keywords_lower.append(keyword.lower())
    self.regular_args = []
    for arg in argv:
      if (not arg.startswith("--")):
        self.regular_args.append(arg)
        continue
      flds = arg[2:].split("=")
      assert len(flds) in (1,2)
      try:
        if (case_sensitive):
          i = keywords.index(flds[0])
        else:
          i = keywords_lower.index(flds[0].lower())
      except ValueError:
        raise RuntimeError("Unknown option: " + arg)
      if (len(flds) == 1):
        setattr(self, keywords[i], True)
      else:
        setattr(self, keywords[i], flds[1])
      self.n += 1

def parse_options_with_chunk(argv, keywords=[], case_sensitive=True):
  flags = parse_options(
    argv=argv,
    keywords=["ChunkSize", "ChunkMember"] + keywords,
    case_sensitive=case_sensitive)
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
