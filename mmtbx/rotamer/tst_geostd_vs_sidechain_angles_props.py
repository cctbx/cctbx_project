from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library.server
import mmtbx.rotamer.sidechain_angles
import sys, os
op = os.path

def run(args):
  assert len(args) == 0
  mon_lib_srv = mmtbx.monomer_library.server.server()
  SidechainAngles = mmtbx.rotamer.sidechain_angles.SidechainAngles(
    show_errs=True)
  assert SidechainAngles.get_rotamer_expectation_frequencies('GLY', 't') is None
  assert SidechainAngles.get_rotamer_expectation_frequencies('CYS', 't')=='26.33%'
  assert SidechainAngles.get_rotamer_expectation_frequencies('CYS', 'w') is None
  for resname_dot_tor_id,atom_names_raw \
        in SidechainAngles.atomsForAngle.items():
    resname, tor_id = resname_dot_tor_id.split(".")
    resname = resname.upper()
    atom_names = tuple([atom_name.split(";")[-1].strip()
      for atom_name in atom_names_raw])
    comp = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
    for tor in comp.tor_list:
      if (tor.id == tor_id):
        if (tor.atom_ids() == atom_names):
          annotation = "OK"
        else:
          annotation = "MISMATCH"
        break
    else:
      annotation = "MISSING"
    print(resname, tor_id, atom_names, annotation)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
