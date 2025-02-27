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
  # lysozyme = """KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"""
  lysozyme = """KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSR"""
  homologs = rcsb_web_services.sequence_search(lysozyme, d_max=2.0)
  assert (len(homologs) > 500)
  atp_binding = rcsb_web_services.chemical_id_search("ATP", protein_only=True)
  # print("len(atp_binding)", len(atp_binding), atp_binding) # 1389
  assert (len(atp_binding) > 1000)
  atp_binding = rcsb_web_services.chemical_id_search("ATP", xray_only=True, protein_only=True)
  # print("len(atp_binding)", len(atp_binding)) # 1389
  assert (len(atp_binding) > 1000)
  report = rcsb_web_services.get_high_resolution_for_structures(atp_binding)
  assert (len(report) == len(atp_binding)) and (len(report[0]) == 2)
  # print (report)
  report = rcsb_web_services.get_high_resolution_and_residue_count_for_structures(atp_binding)
  assert (len(report) == len(atp_binding)) and (len(report[0]) == 3)
  # print (report)
  report = rcsb_web_services.get_r_work_rfree_for_structures(['1ucs', '1yjp'])
  assert report == [['1UCS', 0.133, 0.155], ['1YJP', 0.18086, 0.19014]]
  ligand_info = rcsb_web_services.get_ligand_info_for_structures(['1mru'])
  # print (ligand_info)
  reference_ligand_info = [
    ['1MRU', 'A', 'MG', 24.305, 'Mg', 'MAGNESIUM ION', '[Mg+2]'],
    ['1MRU', 'B', 'MG', 24.305, 'Mg', 'MAGNESIUM ION', '[Mg+2]'],
    ['1MRU', 'A', 'AGS', 523.247, 'C10 H16 N5 O12 P3 S', 'PHOSPHOTHIOPHOSPHORIC ACID-ADENYLATE ESTER',
        'c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=S)(O)O)O)O)N'],
    ['1MRU', 'B', 'AGS', 523.247, 'C10 H16 N5 O12 P3 S', 'PHOSPHOTHIOPHOSPHORIC ACID-ADENYLATE ESTER',
        'c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=S)(O)O)O)O)N']]
  for answer in reference_ligand_info:
    assert answer in ligand_info, "%s\n%s" % (answer, ligand_info)

def exercise_2():
  fes_binding = rcsb_web_services.chemical_id_search(
      "FES",
      xray_only=False,
      # log=sys.stdout
      )
  assert len(fes_binding) > 765, len(fes_binding)
  # checking correct work when no results found
  di1_examples = rcsb_web_services.chemical_id_search(
      'Di2',
      data_only=True,
      sort_by_resolution=True,
      # log=sys.stdout,
      )
  assert len(di1_examples) == 0, len(di1_examples)

def exercise_3():
  # Direct test of post_query
  r = rcsb_web_services.post_query(xray_only=False)
  assert len(r) > 1
  r = rcsb_web_services.post_query()
  assert len(r) > 1
  r = rcsb_web_services.post_query(data_only=True)
  assert len(r) > 1
  r = rcsb_web_services.post_query(d_max=1)
  assert len(r) > 1
  r = rcsb_web_services.post_query(d_min=5)
  assert len(r) > 1
  r = rcsb_web_services.post_query(sort_by_resolution=True)
  assert len(r) > 180000
  r = rcsb_web_services.post_query(clashscore_range=(1,5))
  print('n clashscore filter:', len(r))
  assert len(r) > 80000
  r = rcsb_web_services.post_query(rama_outliers_range=(1,2))
  print('n rama filter:', len(r))
  assert len(r) > 10000
  r = rcsb_web_services.post_query(rota_outliers_range=(0,5))
  print('n rota filter:', len(r))
  assert len(r) > 135000

def exercise_get_emdb_id():
  emdb_ids = rcsb_web_services.get_emdb_id_for_pdb_id('8wcc')
  assert emdb_ids == ['EMD-37438']
  emdb_ids = rcsb_web_services.get_emdb_id_for_pdb_id('1yjp')
  assert emdb_ids == None

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
    exercise_3()
    exercise_get_emdb_id()
    print("OK")
  else:
    print("OK but skipped.")
