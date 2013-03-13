import os, sys
import libtbx.load_env

curated_repository_dir = os.path.dirname(libtbx.env.dist_path("elbow"))
curated_repository_dir = os.path.join(curated_repository_dir,
                                      "chem_data",
                                      "curated_datasets",
                                      )
def imp_dataset(name, verbose=False):
  import imp
  try:
    fp, pathname, description = imp.find_module(name, [curated_repository_dir])
  except:
    assert 0
  if verbose:
      print name
      print fp
      print pathname
      print description
  m = imp.load_module(name, fp, pathname, description)
  return m

top8000_data = imp_dataset("top8000").data

def generate_top8000():
  for row in top8000_data["data"]:
    yield row[1],row[2]

def run():
  print 'curated_repository_dir',curated_repository_dir
  for i, (pdb_code, chain_id) in enumerate(generate_top8000()):
    print i, pdb_code, chain_id
    if i>10: break

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run()
