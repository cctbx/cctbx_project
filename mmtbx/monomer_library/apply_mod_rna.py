import mmtbx.monomer_library.server
import sys

def run(args):
  if (len(args) == 0):
    args = ["mod_rnaC3.cif", "mod_rnaC2.cif", "mod_rnaEsd.cif"]
  for file_name in args:
    mon_lib_srv = mmtbx.monomer_library.server.server()
    cif_object = mmtbx.monomer_library.server.read_cif(file_name=file_name)
    mon_lib_srv.process_cif_object(cif_object=cif_object, file_name=file_name)
    for code in ["AR", "CR", "GR", "UR"]:
      gr = mon_lib_srv.get_comp_comp_id_direct(comp_id=code)
      file_name_orig = "original_%s_show" % code
      print "writing", file_name_orig
      gr.show(f=open(file_name_orig, "w"))
      for key in ["rnaC3", "rnaC2", "rnaEsd"]:
        mod = mon_lib_srv.mod_mod_id_dict.get(key)
        if (mod is not None):
          gr_mod = gr.apply_mod(mod)
          file_name_mod = "modified_%s_%s_show" % (code, key)
          print "writing", file_name_mod
          gr_mod.show(f=open(file_name_mod, "w"))

if (__name__ == "__main__"):
  run(sys.argv[1:])
