
from __future__ import division
import os, sys
import libtbx.load_env
from libtbx.utils import Sorry

curated_repository_dir = libtbx.env.find_in_repositories(
  relative_path="chem_data/curated_datasets",
  test=os.path.isdir)

def imp_dataset(name, verbose=False):
  import imp
  try:
    fp, pathname, description = imp.find_module(name, [curated_repository_dir])
  except Exception :
    assert 0
  if verbose:
      print name
      print fp
      print pathname
      print description
  m = imp.load_module(name, fp, pathname, description)
  return m

def generate_top8000():
  top8000_data = imp_dataset("top8000").data
  for row in top8000_data["data"]:
    yield row[1],row[2]

def run():
  if (curated_repository_dir is None) :
    raise Sorry("chem_data/curated_datasets is not available.")
  print 'curated_repository_dir',curated_repository_dir
  for i, (pdb_code, chain_id) in enumerate(generate_top8000()):
    print i, pdb_code, chain_id
    if i>10: break

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run()
