
from __future__ import absolute_import, division, print_function

def run():
  from mmtbx.monomer_library import server
  from libtbx.utils import Sorry
  try :
    cif1 = server.mon_lib_list_cif()
    cif2 = server.geostd_list_cif()
    cif3 = server.mon_lib_ener_lib_cif()
    assert (not None in [cif1, cif2, cif3])
  except server.MonomerLibraryServerError :
    raise Sorry("The monomer library installation could not be found.  If "+
      "you are using a Phenix installer, we recommend downloading the "+
      "installer again and reinstalling.")
  else :
    print("OK")

if (__name__ == "__main__"):
  run()
