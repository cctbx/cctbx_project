import os, sys

cif_keyword_dictionary = {
  "_chem_comp" : {"id" : str,
                  "name" : str,
                  "type" : str,
                  "pdbx_type" : str,
                  "formula" : str,
                  "mon_nstd_parent_comp_id" : str,
                  "pdbx_synonyms" : str,
                  "pdbx_formal_charge" : int,
                  "pdbx_initial_date" : str,
                  "pdbx_modified_date" : str,
                  "pdbx_ambiguous_flag" : str,
                  "pdbx_release_status" : str,
                  "pdbx_replaced_by" : str,
                  "pdbx_replaces" : str,
                  "formula_weight" : float,
                  "one_letter_code" : str,
                  "three_letter_code" : str,
                  "pdbx_model_coordinates_details" : str,
                  "pdbx_model_coordinates_missing_flag" : str,
                  "pdbx_ideal_coordinates_details" : str,
                  "pdbx_ideal_coordinates_missing_flag" : str,
                  "pdbx_model_coordinates_db_code" : str,
                  "pdbx_processing_site" : str,
                  # added 11/2008
                  "pdbx_subcomponent_list" : str,
                  },
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
                       # added 11/2008
                       "pdbx_component_atom_id" : str,
                       "pdbx_component_comp_id" : str,
                       },
  "_chem_comp_bond" : {"comp_id": str,
                       "atom_id_1": str,
                       "atom_id_2": str,
                       "value_order": str,
                       "pdbx_aromatic_flag": str,
                       "pdbx_stereo_config": str,
                       "pdbx_ordinal": int,
                       },
  "_pdbx_chem_comp_descriptor" : {"comp_id" : str,
                                  "type" : str,
                                  "program" : str,
                                  "program_version" : str,
                                  "descriptor" : str,
                                  },
  "_pdbx_chem_comp_identifier" : {"comp_id" : str,
                                  "type" : str,
                                  "program" : str,
                                  "program_version" : str,
                                  "identifier" : str,
                                  },
  }

class empty(object):
  def __repr__(self):
    outl = "\nObject"
    for attr in self.__dict__:
      outl += "\n  %s : %s" % (attr, getattr(self, attr))
    return outl

  def __len__(self):
    return len(self.__dict__.keys())

def smart_split_cif_line(line):
  line = " %s " % line
  tmp = []
  delimiters = [" ", "'", '"']
  while line:
    delimiter=" "
    if line[0] in delimiters:
      delimiter=line[0]
    start = line.find(delimiter)
    finish = line.find(delimiter, start+1)
    items = line.split(delimiter)
    items = filter(None, items)
    if not items: break
    item=items[0]
    tmp.append(item)
    if finish==-1:
      if len(items)==2:
        tmp.append(items[1])
      break
    if delimiter==" ":
      line = line[len(item)+1:].strip()
    else:
      line = line[len(item)+3:].strip()
  return tmp

def run(filename):
  if not os.path.exists(filename): return None
  lines = open(filename).read().splitlines()
  merge = True
  while merge:
    merge = []
    for i, line in enumerate(lines):
      if line.find(";")==0:
        if i not in merge:
          merge.append(i)
      elif line.find('"')==0:
        if i not in merge:
          merge.append(i)
      elif line.find("_chem_comp.name")==0:
        if len(line.split())==1:
          if i+1 not in merge:
            merge.append(i+1)
      elif line.find("_chem_comp.pdbx_synonyms")==0:
        if len(line.split())==1:
          if i+1 not in merge:
            merge.append(i+1)
    if merge:
      merge.reverse()
      for i in merge:
        if lines[i][1:].strip():
          if lines[i][0] in [";"]:
            lines[i-1] += ' "%s"' % lines[i][1:].strip()
          else:
            lines[i-1] += ' "%s"' % lines[i].strip()
        del lines[i]

  line_iter = iter(lines)

  complete_cif_data = {}
  loop_list = []
  non_loop = {}
  code = ""
  pdbx_reading = False
  remove_loop_fields = {}
  loop_index = 0
  for line in line_iter:
    line = line.replace("\t", " "*8)
    if line.find("#")==0: continue
    if line.find("_pdbx")==0: pdbx_reading = True
    line = "%s  " % line
    if False:
      while line.find('"')>-1:
        start = line.find('"')
        finish = line.find('"', start+1)
        if finish==-1: break
        line = "%s%s%s" % (line[:start],
                           line[start+1:finish].replace(" ", "_space_"),
                           line[finish+1:])
      if pdbx_reading:
        while line.find("'")>-1:
          start = line.find("'")
          finish = line.find("'", start+1)
          if finish==-1: break
          line = "%s%s%s" % (line[:start],
                             line[start+1:finish].replace(" ", "_space_"),
                             line[finish+1:])
    if line.find("loop_")==0:
      loop_list = []
      loop_index = 0
      continue
    if line.find(";")==0: continue
    if line.find("_")==0:
      for cif_key in cif_keyword_dictionary:
        if line.find(cif_key)==-1: continue
        for attr in cif_keyword_dictionary[cif_key]:
          test = "%s.%s " % (cif_key, attr)
          if line.find(test)>-1:
            loop_list.append(attr)
            if len(line.split())>1:
              non_loop[test] = line.split()[1:]
            break
        if line.find(test)>-1: break
      else:
        tmp = line.split()[0]
        cif_key, attr = tmp.split(".")
        remove_loop_fields.setdefault(cif_key, [])
        remove_loop_fields[cif_key].append(loop_index)
      loop_index+=1
    else:
      if loop_list:
        if code and line.find(code)!=0: continue
        obj = empty()
        i=0
        #for ptr, item in enumerate(line.split()):
        for ptr, item in enumerate(smart_split_cif_line(line)):
          if ptr in remove_loop_fields.get(cif_key, []): continue
          if len(loop_list)<=i:
            print 'Problem with CIF line parsing'
            print line
            print loop_list
            continue
          if loop_list[i]=="comp_id": code = item
          if item not in ["?", "."]:
            item = cif_keyword_dictionary[cif_key][loop_list[i]](item)
          if type(item)==type(""):
            if item[0]=='"' and item[-1]=='"':
              item = item[1:-1]
            item = item.replace("_space_", " ")
          setattr(obj, loop_list[i], item)
          i+=1
        complete_cif_data.setdefault(cif_key, [])
        if hasattr(obj, "alt_atom_id"):
          if len(obj.alt_atom_id)>5: assert 0
        complete_cif_data[cif_key].append(obj)
  # non loop parsing
  for cif_key in cif_keyword_dictionary:
    obj = empty()
    for key in non_loop:
      if key.find("%s." % cif_key)==0:
        if len(non_loop[key])==1:
          item = non_loop[key][0]
        else:
          item = " ".join(non_loop[key])
        sk = key.split(".")[1].strip()
        if item not in ["?", "."]:
          try: item = cif_keyword_dictionary[cif_key][sk](item)
          except Exception:
            print key
            print sk
            print cif_key
            print cif_keyword_dictionary[cif_key].keys()
            raise
        setattr(obj, sk, item)
    if len(obj):
      complete_cif_data[cif_key] = [obj]

  return complete_cif_data

if __name__=="__main__":
  import os, sys
  cif = run(sys.argv[1])
  for key in cif:
    print '_'*80
    print key
    for item in cif[key]:
      print item
