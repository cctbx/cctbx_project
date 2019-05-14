
# XXX depends on internet connectivity, so not run as part of standard tests

from __future__ import absolute_import, division, print_function

def exercise():
  from mmtbx.wwpdb import rcsb_web_services
  lysozyme = """KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"""
  homologs = rcsb_web_services.sequence_search(lysozyme, d_max=2.0)
  assert (len(homologs) > 500)
  atp_binding = rcsb_web_services.chemical_id_search("ATP", protein_only=True)
  assert (len(atp_binding) > 650)
  report = rcsb_web_services.get_high_resolution_for_structures(atp_binding)
  assert (len(report) == len(atp_binding)) and (len(report[0]) == 2)
  ligand_info = rcsb_web_services.get_ligand_info_for_structures(['1mru'])
  assert (len(ligand_info) == 4)
  print("OK")

if (__name__ == "__main__"):
  exercise()
