
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

top8000_data = imp_dataset("top8000").data
rna11_data = imp_dataset("rna11").data
hiq54_data = imp_dataset("hiq54").data

hiq54_data["info"] = """
The HiQ54 dataset for tests in methods development: small (60-200 residues), 
monomeric, non-redundant, high-quality (both resolution and MolProbity score 
both >=1.4), with no tightly-bound ligands.  

Literature reference: Leaver-Fay et al (2013) Meth Enzymol 523: 109-143."""

class dataset_item(object):
  def __init__(self, attrs, items):
    for i, attr in enumerate(attrs):
      setattr(self, attr, items[i])

  def __repr__(self):
    outl = "dataset_item"
    for attr in self.__dict__:
      outl += "\n  %-20s : %s" % (attr, getattr(self, attr))
    return outl

def generate_top8000():
  top8000_data = imp_dataset("top8000").data
  for row in top8000_data["data"]:
    yield dataset_item(top8000_data["headers"], row)

def generate_rna11():
  for row in rna11_data["data"]:
    yield dataset_item(rna11_data["headers"], row)

def generate_hiq54():
  print hiq54_data.get("info", "")
  for row in hiq54_data["data"]:
    yield dataset_item(hiq54_data["headers"], row)

def run():
  if (curated_repository_dir is None) :
    raise Sorry("chem_data/curated_datasets is not available.")
  print '\nCurated Repository Dir',curated_repository_dir

  print '\n','Top 8000 '*10
  for i, item in enumerate(generate_top8000()):
    print '-'*80
    print i, item
    if i>0: break
  print '\n','RNA 11 '*10
  for i, item in enumerate(generate_rna11()):
    print '-'*80
    print i, item
    if i>0: break
  print '\n','HiQ54 '*10
  for i, item in enumerate(generate_hiq54()):
    print '-'*80
    print i, item
    if i>0: break

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run()
