from __future__ import absolute_import, division, print_function
import os
import mmtbx.model
import libtbx.load_env
from mmtbx import monomer_library
import mmtbx.monomer_library.server
from mmtbx import ias
import iotbx.pdb
import copy

def exercise():
  # initial setup
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  cif_file = libtbx.env.find_in_repositories(
            relative_path="phenix_regression/pdb/tyr.cif", test=os.path.isfile)
  cif_objects = []
  mon_lib_srv.process_cif(file_name= cif_file)
  ener_lib.process_cif(file_name= cif_file)
  pdb_file = libtbx.env.find_in_repositories(
            relative_path="phenix_regression/pdb/ygg.pdb", test=os.path.isfile)
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)

  # # get ias manager
  params = ias.ias_master_params.extract()
  all = ias.ias_master_params.extract().build_ias_types
  for opt in [["L"], ["R"], ["B"], ["BH"], all]:
    mol = mmtbx.model.manager(
        model_input=copy.deepcopy(pdb_inp),
        restraint_objects=cif_objects)
    mol.process(make_restraints=True)
    if(opt == ["L"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      f = open("zz.pdb","w")
      pdb_str = mol.model_as_pdb()
      f.write(pdb_str)
      f.close()
      ias_selection = mol.ias_manager.get_ias_selection()
      assert ias_selection.size() == 47, ias_selection.size()
      assert ias_selection.count(True) == 6,   ias_selection.count(True)
      assert ias_selection.count(False) == 41, ias_selection.count(False)
    if(opt == all):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      ias_selection = mol.ias_manager.get_ias_selection()
      assert ias_selection.size() == 88, ias_selection.size()
      assert ias_selection.count(True) == 47
      assert ias_selection.count(False) == 41
    if(opt == ["R"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      ias_selection = mol.ias_manager.get_ias_selection()
      assert ias_selection.size() == 42
      assert ias_selection.count(True) == 1
      assert ias_selection.count(False) == 41
    if(opt == ["B"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      ias_selection = mol.ias_manager.get_ias_selection()
      assert ias_selection.size() == 62
      assert ias_selection.count(True) == 21
      assert ias_selection.count(False) == 41
    if(opt == ["BH"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      ias_selection = mol.ias_manager.get_ias_selection()
      assert ias_selection.size() == 60
      assert ias_selection.count(True) == 19
      assert ias_selection.count(False) == 41

def run():
  exercise()

if (__name__ == "__main__"):
  run()
  print("OK")
