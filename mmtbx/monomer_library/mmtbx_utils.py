from __future__ import division

import os
import libtbx.load_env
from mmtbx.monomer_library.server import find_mon_lib_file

phenix_repository_dir = os.path.dirname(libtbx.env.dist_path("phenix"))
lib_paths = [
  #os.path.join("chem_data", "monomers"),
  #os.path.join("chem_data", "geostd"),
  os.path.join("chem_data", "mon_lib"),
  "mon_lib",
  os.path.join("ext_ref_files", "mon_lib"),
  ]
mon_lib_directory = os.path.join(phenix_repository_dir, lib_paths[0])

def get_residue_synonyms():
  synonyms = {}
  for relative_path in lib_paths:
    mon_lib_path=libtbx.env.find_in_repositories(
      relative_path=relative_path)
    if mon_lib_path: break
  else:
    return {}
  list_cif = os.path.join(mon_lib_path,
                          "list",
                          "mon_lib_list.cif",
                          )
  f = open(list_cif, "r")
  lines = f.readlines()
  f.close()

  import iotbx.cif
  cif = iotbx.cif.reader(input_string="\n".join(lines),
                         strict=False).model()
  if cif.get("comp_synonym_list", None):
    for a, b in zip(cif["comp_synonym_list"]["_chem_comp_synonym.comp_id"],
                    cif["comp_synonym_list"]["_chem_comp_synonym.comp_alternative_id"],
                    ):
      synonyms[b]=a
      if a.find("-"):
        synonyms[b] = a.split("-")[0]

  if 0:
    syn_read = False
    for line in lines:
      if line.find("#")>-1:
        line = line[:line.find("#")]
      if syn_read:
        tmp = line.split()
        if len(tmp)<3:
          syn_read=False
          continue
        synonyms[tmp[1]] = tmp[0]
        if tmp[0].find("-"):
          synonyms[tmp[1]] = tmp[0].split("-")[0]
      if line.find("_chem_comp_synonym.")>-1:
        syn_read=True
  return synonyms

def _get_monomer_cif_file_basic(code, old_fashion=False):
  #code = monomer_library_filename_renames.get(code, code)
  if old_fashion:
    for relative_path in lib_paths:
      mon_lib_path=libtbx.env.find_in_repositories(
        relative_path=relative_path)
      if mon_lib_path: break
    else:
      return None
    for file_name in os.listdir(os.path.join(mon_lib_path,
                                             code[0].lower(),
                                             )
                                ):
      if file_name.startswith(code):
        file_name = os.path.join(mon_lib_path,
                                 code[0].lower(),
                                 file_name,
                                 )
        break
    else:
      return None
  else:
    file_name = find_mon_lib_file(
      relative_path_components=[code[0].lower(),
                                "data_%s.cif" % code.upper(),
                                ]
      )
    if file_name is None:
      file_name = find_mon_lib_file(
        relative_path_components=[code[0].lower(),
                                  "%s.cif" % code.upper(),
                                  ]
        )
    if file_name is None:
      for d in os.listdir(mon_lib_directory):
        if not os.path.isdir(os.path.join(mon_lib_directory, d)): continue
        for tf in os.listdir(os.path.join(mon_lib_directory, d)):

          ptr = tf.lower()[:-3].find("%s" % code.lower())
          if ptr>-1:
            ptr += len(code)
            if tf[ptr] in ['-', '_']:
              file_name = os.path.join(mon_lib_directory, d, tf)
              break
  return file_name

def get_monomer_cif_file(code, old_fashion=False):
  synonyms = get_residue_synonyms()
  file_name = _get_monomer_cif_file_basic(code, old_fashion=old_fashion)
  if file_name and os.path.exists(file_name): return file_name
  if code in synonyms:
    code = synonyms[code]
  file_name = _get_monomer_cif_file_basic(code, old_fashion=old_fashion)
  return file_name

