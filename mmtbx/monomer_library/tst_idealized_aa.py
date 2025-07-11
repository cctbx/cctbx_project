from __future__ import absolute_import, division, print_function
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx.monomer_library import idealized_aa as iaa
import iotbx.pdb.amino_acid_codes
from scitbx.array_family import flex
import time

def exercise_00():
  d = iaa.residue_dict()
  assert len(d) == 48
  for aac in iotbx.pdb.amino_acid_codes.one_letter_given_three_letter:
    assert aac.lower() in d
  #
  mon_lib_srv = monomer_library.server.server()
  ener_lib    = monomer_library.server.ener_lib()
  for aac in iotbx.pdb.amino_acid_codes.one_letter_given_three_letter:
    aac = aac.lower()
    for aac_ in [aac, aac+"_h"]:
      residue_as_string = iaa.__dict__[aac_]
      rs = flex.std_string(residue_as_string.splitlines())
      params = monomer_library.pdb_interpretation.master_params.extract()
      params.flip_symmetric_amino_acids=False
      processed_pdb_file = monomer_library.pdb_interpretation.process(
         mon_lib_srv = mon_lib_srv,
         ener_lib    = ener_lib,
         params      = params,
         raw_records = rs)
      sites_cart = processed_pdb_file.xray_structure().sites_cart()
      grm = processed_pdb_file.geometry_restraints_manager(
        show_energies = False, plain_pairs_radius = 5.0)
      es = grm.energies_sites(
        sites_cart = sites_cart)
      b = es.bond_deviations()
      a = es.angle_deviations()
      b_z = es.bond_deviations_z()
      a_z = es.angle_deviations_z()
      print("%5s"%aac_, "bonds     : %5.3f %5.3f %5.3f %d"%b, \
        "angles     : %5.3f %5.3f %5.3f %d"%a)
      assert a[2] < 1.2, a[2]
      assert b[2] < 0.006, b[2]
      print("%5s"%aac_, "bonds rmsZ: %5.3f %5.3f %5.3f %d"%b_z, \
        "angles rmsZ: %5.3f %5.3f %5.3f %d"%a_z)
      assert a_z[2] < 0.7, a_z[2]
      assert b_z[2] < 0.7, b_z[2]


if (__name__ == "__main__"):
  t0 = time.time()
  exercise_00()
  print("Time: %6.3f"%(time.time()-t0))
