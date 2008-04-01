import os, sys

cif_keyword_dictionary = {
  "_chem_comp_atom" : {"comp_id": str,
                       "atom_id": str,
                       "alt_atom_id": str,
                       "type_symbol": str,
                       "charge" : int,
                       "pdbx_align": int,
                       "pdbx_aromatic_flag": str,
                       "pdbx_leaving_atom_flag": str,
                       "pdbx_stereo_config": str,
                       "model_Cartn_x": float,
                       "model_Cartn_y": float,
                       "model_Cartn_z": float,
                       "pdbx_model_Cartn_x_ideal": float,
                       "pdbx_model_Cartn_y_ideal": float,
                       "pdbx_model_Cartn_z_ideal": float,
                       "pdbx_ordinal" : int,
                       },
  "_chem_comp_bond" : {"comp_id": str,
                       "atom_id_1": str,
                       "atom_id_2": str,
                       "value_order": str,
                       "pdbx_aromatic_flag": str,
                       "pdbx_stereo_config": str,
                       "pdbx_ordinal": int,
                       },
  }

class empty(object):
  def __repr__(self):
    outl = "\nObject"
    for attr in self.__dict__:
      outl += "\n  %s : %s" % (attr, getattr(self, attr))
    return outl

def run(filename):
  if not os.path.exists(filename): return None
  lines = open(filename).read().splitlines()
  line_iter = iter(lines)

  complete_cif_data = {}
  loop_list = []
  for line in line_iter:
    start = line.find('"')
    finish = line.find('"', start+1)
    line = line[:start] + line[start:finish].replace(" ", "_") + line[finish:]
    if line.find("#")==0: continue
    if line.find("loop_")==0:
      loop_list = []
      continue
    if line.find("_")==0:
      for cif_key in cif_keyword_dictionary:
        if line.find(cif_key)==-1: continue
        for attr in cif_keyword_dictionary[cif_key]:
          test = "%s.%s" % (cif_key, attr)
          if line.find(test)>-1:
            loop_list.append(attr)
            break
        if line.find(test)>-1: break
    else:
      if loop_list:
        obj = empty()
        for i, item in enumerate(line.split()):
          if item not in ["?", "."]:
            item = cif_keyword_dictionary[cif_key][loop_list[i]](item)
          if type(item)==type(""):
            if item[0]=='"' and item[-1]=='"':
              item = item[1:-1]
              item = item.replace("_", " ")
          setattr(obj, loop_list[i], item)
        complete_cif_data.setdefault(cif_key, [])
        complete_cif_data[cif_key].append(obj)

  return complete_cif_data

if __name__=="__main__":
  import os, sys
  run(sys.argv[1])
