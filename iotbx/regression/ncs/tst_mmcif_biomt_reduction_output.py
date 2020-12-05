from __future__ import absolute_import, division, print_function

import iotbx.pdb
import mmtbx.model
from libtbx.test_utils import assert_lines_in_text
from iotbx.regression.ncs.tst_mtrix_biomt_cmdl import pdb_str_0
from cctbx.maptbx.box import shift_and_box_model

def test_1():
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str_0)
  model = mmtbx.model.manager(model_input=inp)
  model.expand_with_BIOMT_records()
  model = shift_and_box_model(model)

  sel = model.selection("chain '0' or chain 'C' or chain 'F' or chain 'I' or chain 'L' or chain 'O' or chain 'R' or chain 'U' or chain 'X'")
  model = model.select(sel)

  model.search_for_ncs()
  model.setup_ncs_constraints_groups(filter_groups=True)
  n1 = model.get_number_of_atoms()
  assert n1 == 648, n1
  assert model.ncs_constraints_present()
  nrgl = model.get_ncs_groups()
  assert len(nrgl[0].master_iselection) == 72
  assert len(nrgl[0].copies) == 8
  # nrgl._show()
  # print (model.can_be_reduced_with_biomt())
  cif_txt = model.model_as_mmcif(try_reduce_with_biomt=True)
  # print (cif_txt)
  assert_lines_in_text(cif_txt, """
loop_
  _pdbx_struct_assembly_gen.assembly_id
  _pdbx_struct_assembly_gen.oper_expression
  _pdbx_struct_assembly_gen.asym_id_list
   1 (1-9) A""")
  assert_lines_in_text(cif_txt, """
loop_
  _pdbx_struct_assembly.id
  _pdbx_struct_assembly.details
  _pdbx_struct_assembly.method_details
  _pdbx_struct_assembly.oligomeric_details
  _pdbx_struct_assembly.oligomeric_count
   1 'Symmetry assembly' ? ? ? """)
  assert_lines_in_text(cif_txt, """
  _pdbx_struct_oper_list.vector[1]
  _pdbx_struct_oper_list.vector[2]
  _pdbx_struct_oper_list.vector[3]
   1 'point symmetry operation' ? ? 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0""")


  inp = iotbx.pdb.input(source_info=None, lines=cif_txt)
  m2 = mmtbx.model.manager(model_input=inp)
  n2_1 = m2.get_number_of_atoms()
  assert n2_1 == 72
  m2.expand_with_BIOMT_records()
  n2_2 = m2.get_number_of_atoms()
  # print (n1, n2)
  assert n1 == n2_2, "%d, %d" % (n1, n2)

def test_2():
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str_0)
  model = mmtbx.model.manager(model_input=inp)
  model.expand_with_BIOMT_records()
  model = shift_and_box_model(model)
  model.search_for_ncs()
  model.setup_ncs_constraints_groups(filter_groups=True)
  assert model.ncs_constraints_present()
  assert not model.can_be_reduced_with_biomt()
  assert "" == model.model_as_mmcif(try_reduce_with_biomt=True)

if __name__ == '__main__':
  import libtbx.load_env
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_1()
    test_2()
    print('OK')
  else:
    print('chem_data is not available, skipping')
