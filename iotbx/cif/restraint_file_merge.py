from __future__ import division, print_function
import os, sys

import iotbx
from iotbx import cif

def merge_block(b1, b2):
  for key, item in b1.items():
    for i in b2[key]:
      item.append(i)

def run(filenames,
        filename,
        no_file_access=False,
        no_file_but_save=False,
        verbose=False,
        ):
  outl = "  Joining"
  for f in filenames:
    if f.find("\n")==-1:
      outl += "\n\t%s" % f
    else:
      outl += "\n\t%s lines" % (len(f.split("\n")))
  outl += "\n  into\n\t%s" % filename
  if not no_file_access:
    if verbose: print(outl)

  cifs = []
  if not no_file_access:
    for cif_file_name in filenames:
      if not os.path.exists(cif_file_name): continue
      cif_object = iotbx.cif.reader(file_path=cif_file_name).model()
      cifs.append(cif_object)
  else:
    for lines in filenames:
      cif_object = iotbx.cif.reader(input_string=lines).model()
      cifs.append(cif_object)

  comp_list = None
  for i, cif_object in enumerate(cifs):
    if comp_list is None:
      comp_list = cif_object['comp_list']
    else:
      merge_block(comp_list, cif_object['comp_list'])
    del cif_object['comp_list']

  outl = 'data_comp_list\n'
  outl += str(comp_list)
  for i, cif_object in enumerate(cifs):
    outl += str(cif_object)

  if verbose:
    print(outl)

  if not no_file_access or no_file_but_save:
    f = open(filename, "w")
    f.write(outl)
    f.close()
  else:
    return outl

if __name__=="__main__":
  if len(sys.argv)<4:
    print("""
  usage:
    iotbx.cif.restraints_file_merge.py merge1 merge2 ... target
""")
    sys.exit()
  target = sys.argv[-1]
  filenames=[]
  for i in range(1,len(sys.argv)-1):
    filenames.append(sys.argv[i])
  run(filenames, target)
