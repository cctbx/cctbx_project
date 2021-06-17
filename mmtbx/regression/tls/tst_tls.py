from __future__ import absolute_import, division, print_function
from mmtbx.tls import tools
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import iotbx.pdb.remark_3_interpretation
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
import os
from six.moves import zip

params = monomer_library.pdb_interpretation.master_params.extract()
params.flip_symmetric_amino_acids = False

def uaniso_from_tls_and_back():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1OC2_tst.pdb",
    test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       params = params,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  xray_structure = processed_pdb_file.xray_structure()
  selections = []
  for string in ["chain A", "chain B"]:
      selections.append(processed_pdb_file.all_chain_proxies.selection(
                                                              string = string))
  input_tls_data = iotbx.pdb.remark_3_interpretation.extract_tls_parameters(
    remark_3_records=processed_pdb_file.all_chain_proxies.pdb_inp
      .extract_remark_iii_records(iii=3),
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    chain_ids=[]).tls_params
  tls_params = []
  for item in input_tls_data:
      tls_params.append(tools.tlso(t      = item.t,
                                   l      = item.l,
                                   s      = item.s,
                                   origin = item.origin))
  tools.show_tls(tlsos = tls_params)
  u_cart_from_tls = tools.u_cart_from_tls(
                        sites_cart = xray_structure.sites_cart(),
                        selections = selections,
                        tlsos      = tls_params)

  i = 0
  for utls,atom in zip(u_cart_from_tls,
                       processed_pdb_file.all_chain_proxies.pdb_atoms):
    updb = atom.uij
    #i += 1
    #print "      ", i
    #print "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f "% \
    #                         (utls[0],utls[1],utls[2],utls[3],utls[4],utls[5])
    #print "%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f "% \
    #                         (updb[0],updb[1],updb[2],updb[3],updb[4],updb[5])
    assert approx_equal(utls,updb, 1.e-4)

  tlsos_initial = []
  for input_tls_data_ in input_tls_data:
      tlsos_initial.append(tools.tlso(
                                   t = ([0.0,0.0,0.0,0.0,0.0,0.0]),
                                   l = ([0.0,0.0,0.0,0.0,0.0,0.0]),
                                   s = ([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),
                                   origin = input_tls_data_.origin))
  tls_from_uanisos = tools.tls_from_uanisos(
    number_of_macro_cycles       = 300,
    max_iterations               = 1000,
    xray_structure = xray_structure,
    selections     = selections,
    tlsos_initial  = tlsos_initial)
  print("\nTLS from Uaniso:\n")
  tools.show_tls(tlsos = tls_from_uanisos)

  for input_tls_data_,tls_from_uanisos_ in zip(tls_params,tls_from_uanisos):
    assert approx_equal(input_tls_data_.t,      tls_from_uanisos_.t, 1.e-4)
    assert approx_equal(input_tls_data_.l,      tls_from_uanisos_.l, 1.e-4)
    assert approx_equal(input_tls_data_.s,      tls_from_uanisos_.s, 1.e-4)
    assert approx_equal(input_tls_data_.origin, tls_from_uanisos_.origin, 1.e-3)
  #
  print(format_cpu_times())

if (__name__ == "__main__"):
  uaniso_from_tls_and_back()
