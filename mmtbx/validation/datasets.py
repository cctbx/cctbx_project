
from __future__ import absolute_import, division, print_function
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
      print(name)
      print(fp)
      print(pathname)
      print(description)
  m = imp.load_module(name, fp, pathname, description)
  return m

top8000_data = imp_dataset("top8000").data
rna11_data = imp_dataset("rna11").data
hiq54_data = imp_dataset("hiq54").data
iridium_data = imp_dataset("iridium").data

hiq54_data["info"] = """
The HiQ54 dataset for tests in methods development: small (60-200 residues),
monomeric, non-redundant, high-quality (both resolution and MolProbity score
both >=1.4), with no tightly-bound ligands.

Literature reference: Leaver-Fay et al (2013) Meth Enzymol 523: 109-143."""

iridium_data["info"] = """
The Iridium dataset for ligands as curated by Nat Echols to have the original
(as best as possible) search model.
"""

class dataset_item(object):
  def __init__(self, attrs, items):
    for i, attr in enumerate(attrs):
      setattr(self, attr, items[i])

  def __repr__(self):
    outl = "dataset_item"
    for attr in self.__dict__:
      outl += "\n  %-20s : %s" % (attr, getattr(self, attr))
    return outl

def generate_data_items(data):
  if data.get("info", ""): print(data.get("info", ""))
  for row in data["data"]:
    yield dataset_item(data["headers"], row)

data_lookups = {}
for key in locals().keys():
  if not key.endswith("_data"): continue
  cmd = "generate_%s = generate_data_items(%s)" % (key[:-5],
                                                   key,
                                                 )
  exec(cmd)
  if key=="top8000_data": data_lookups[key]="Top 8000"
  elif key=="rna11_data": data_lookups[key]="RNA 11"
  elif key=="hiq54_data": data_lookups[key]="HiQ54"
  elif key=="iridium_data": data_lookups[key]="Iridium"
  else: assert 0


def run(local):
  if (curated_repository_dir is None):
    raise Sorry("chem_data/curated_datasets is not available.")
  print('\nCurated Repository Dir',curated_repository_dir)
  for key in local.keys():
    if not key.endswith("_data"): continue
    print('\n',('%s ' % data_lookups[key])*10)
    for i, item in enumerate(local.get("generate_%s" % key[:-5], [])):
      print('-'*80)
      print(i, item)
      if i>0: break

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(locals())
