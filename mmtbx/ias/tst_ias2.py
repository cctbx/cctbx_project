from __future__ import division
import os
import mmtbx.model
import libtbx.load_env
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import ias
import iotbx.pdb
from iotbx import file_reader
import copy

def exercise():
  """
  This is copy-paste of tst_ias.py with use of new way of creating
  mmtbx.model.manager
  Later original test could be removed.
  """
  # initial setup
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  cif_file = libtbx.env.find_in_repositories(
            relative_path="phenix_regression/pdb/tyr.cif", test=os.path.isfile)
  # cif_objects = [(cif_file, file_reader.any_file(cif_file))]
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
        restraint_objects=cif_objects,
        build_grm=True)
    if(opt == ["L"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      mol.write_pdb_file(out=open("zz.pdb","w"))
      assert mol.ias_selection.size() == 47, mol.ias_selection.size()
      assert mol.ias_selection.count(True) == 6, mol.ias_selection.count(True)
      assert mol.ias_selection.count(False) == 41, mol.ias_selection.count(False)
    if(opt == all):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      assert mol.ias_selection.size() == 88,mol.ias_selection.size()
      assert mol.ias_selection.count(True) == 47
      assert mol.ias_selection.count(False) == 41
    if(opt == ["R"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      assert mol.ias_selection.size() == 42
      assert mol.ias_selection.count(True) == 1
      assert mol.ias_selection.count(False) == 41
    if(opt == ["B"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      assert mol.ias_selection.size() == 62
      assert mol.ias_selection.count(True) == 21
      assert mol.ias_selection.count(False) == 41
    if(opt == ["BH"]):
      params.build_ias_types = opt
      mol.add_ias(fmodel = None, ias_params = params)
      assert mol.ias_selection.size() == 60
      assert mol.ias_selection.count(True) == 19
      assert mol.ias_selection.count(False) == 41

def run():
  exercise()

if (__name__ == "__main__"):
  run()
  print "OK"
