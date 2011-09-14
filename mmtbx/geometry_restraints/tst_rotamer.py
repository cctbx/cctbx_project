from mmtbx.monomer_library import server, pdb_interpretation
import mmtbx.geometry_restraints.rotamer
from iotbx import file_reader
import libtbx.load_env
from libtbx.utils import null_out
import cStringIO
import os

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/enk.pdb",
    test=os.path.isfile)
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  params = pdb_interpretation.master_params.extract()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    pdb_inp=file_reader.any_file(pdb_file, force_type="pdb").file_object,
    log=null_out())
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  proxies = mmtbx.geometry_restraints.rotamer.extract_proxies(
    pdb_hierarchy=pdb_hierarchy,
    log=null_out())
  assert (len(proxies) == 3)
  m = mmtbx.geometry_restraints.rotamer.manager(pdb_hierarchy, log=null_out())
  grm = processed_pdb_file.geometry_restraints_manager()
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  tmp_log = cStringIO.StringIO()
  old_angles = []
  for proxy in grm.dihedral_proxies :
    old_angles.append(proxy.angle_ideal)
  n_updates = m.update_dihedral_proxies(
    sites_cart=sites_cart,
    dihedral_proxies=grm.dihedral_proxies,
    verbose=True,
    log=tmp_log)
#  print tmp_log.getvalue()
  assert (tmp_log.getvalue() == """\
 Updating dihedral proxies from rotamer library
   TYR A   1  : Sm30 (rmsd=1.395)
   PHE A   4  : Sm30 (rmsd=0.684)
   LEU A   5  : tp (rmsd=0.004)
 6 proxies modified.
""")
  assert (n_updates == 6)
  for i, proxy in enumerate(grm.dihedral_proxies) :
    if (i == 1) : assert (proxy.angle_ideal == 77.0)
    elif (i == 9) : assert (proxy.angle_ideal == -85.0)

if (__name__ == "__main__") :
  exercise()
  print "OK"
