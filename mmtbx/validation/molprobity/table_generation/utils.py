def make_reskey(model_id, chain_id, resseq, icode):
  """
  Creates a unique string residue key with mmCIF-like placeholders for blank fields.
  The residue key is space-delimited.
  """
  return ' '.join([
      model_id.strip() or '.',
      chain_id.strip() or '.',
      resseq.strip(),
      icode.strip() or '?'
  ])
