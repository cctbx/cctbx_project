from mmtbx.monomer_library import server
import string
import os

def exercise():
  list_cif = server.mon_lib_list_cif()
  srv = server.server(list_cif=list_cif)
  print "srv.root_path:", srv.root_path
  table_of_contents = []
  n_get_comp_comp_id_successes = 0
  for first_char in string.lowercase+string.digits:
    sub_dir = os.path.join(srv.root_path, first_char)
    if (not os.path.isdir(sub_dir)): continue
    for node in os.listdir(sub_dir):
      if (not node.lower().endswith(".cif")): continue
      comp_id = node[:-4]
      if (comp_id.endswith("_EL")): continue
      if (comp_id in ["CON_CON", "PRN_PRN"]):
        comp_id = comp_id[:3]
      if (comp_id.upper() != comp_id):
        print "Mixed case:", os.path.join(first_char, node)
      comp_comp_id = srv.get_comp_comp_id(
        comp_id=comp_id, hide_mon_lib_dna_rna_cif=False)
      if (comp_comp_id is None):
        print "Error instantiating comp_comp_id %s (%s)" % (
          comp_id, os.path.join(sub_dir, node))
      else:
        n_get_comp_comp_id_successes += 1
        table_of_contents.append(
          " ".join([comp_id.upper(), os.path.join(first_char, node)]))
  print "number of cif files read successfully:", n_get_comp_comp_id_successes
  print "writing file table_of_contents"
  open("table_of_contents", "w").write("\n".join(table_of_contents)+"\n")

if (__name__ == "__main__"):
  exercise()
