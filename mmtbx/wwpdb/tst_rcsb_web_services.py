from __future__ import absolute_import, division, print_function
from mmtbx.wwpdb import rcsb_web_services
import requests

def thorough_exercise():
  """
  This exercises all currently available options, so makes 64 requests.
  Took 20 minutes, therefore disabled from the general testing.
  """
  lysozyme = """KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"""
  n_good = 0
  for po in [True, False]:
    for do in [True, False]:
      for sbr in [True, False]:
        for xro in [True, False]:
          for d_max in [None, 3.0]:
            for d_min in [None, 1.0]:
              print ("ARGUMENTS:", po, do, sbr, xro, d_max, d_min)
              homologs = rcsb_web_services.sequence_search(
                  lysozyme,
                  xray_only=xro,
                  d_max=d_max,
                  d_min=d_min,
                  protein_only=po,
                  data_only=do,
                  sort_by_resolution=sbr)
              if len(homologs) > 200:
                n_good += 1
                print ("GOOD:", n_good)
  print (n_good)
  assert n_good == 64

def exercise():
  lysozyme = """KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"""
  homologs = rcsb_web_services.sequence_search(lysozyme, d_max=2.0)
  assert (len(homologs) > 500)
  atp_binding = rcsb_web_services.chemical_id_search("ATP", protein_only=True)
  assert (len(atp_binding) > 650)
  atp_binding = rcsb_web_services.chemical_id_search("ATP", xray_only=True, protein_only=True)
  assert (len(atp_binding) > 650)
  report = rcsb_web_services.get_high_resolution_for_structures(atp_binding)
  assert (len(report) == len(atp_binding)) and (len(report[0]) == 2)
  # print (report)
  report = rcsb_web_services.get_high_resolution_and_residue_count_for_structures(atp_binding)
  assert (len(report) == len(atp_binding)) and (len(report[0]) == 3)
  # print (report)
  ligand_info = rcsb_web_services.get_ligand_info_for_structures(['1mru'])
  # print (ligand_info)
  assert ligand_info == [
    ['1MRU', 'A', 'MG', 24.305, 'Mg', 'MAGNESIUM ION', '[Mg+2]'],
    ['1MRU', 'B', 'MG', 24.305, 'Mg', 'MAGNESIUM ION', '[Mg+2]'],
    ['1MRU', 'A', 'AGS', 523.247, 'C10 H16 N5 O12 P3 S', 'PHOSPHOTHIOPHOSPHORIC ACID-ADENYLATE ESTER',
        'c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=S)(O)O)O)O)N'],
    ['1MRU', 'B', 'AGS', 523.247, 'C10 H16 N5 O12 P3 S', 'PHOSPHOTHIOPHOSPHORIC ACID-ADENYLATE ESTER',
        'c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=S)(O)O)O)O)N']
  ], ligand_info

def exercise_2():
  fes_binding = rcsb_web_services.chemical_id_search(
      "FES",
      xray_only=False,
      # log=sys.stdout
      )
  assert len(fes_binding) > 765, len(fes_binding)
  # checking correct work when no results found
  di1_examples = rcsb_web_services.chemical_id_search(
      'Di1',
      data_only=True,
      sort_by_resolution=True,
      # log=sys.stdout,
      )
  assert len(di1_examples) == 0

if (__name__ == "__main__"):
  # thorough_exercise()
  # check if internet and rcsb are available
  exception_occured = False
  try:
    r = requests.get('https://search.rcsb.org/')
  except Exception:
    print("OK but exception.")
    exception_occured = True
  if not exception_occured and r.ok and len(r.text) > 100:
    exercise()
    exercise_2()
    print("OK")
  else:
    print("OK but skipped.")
