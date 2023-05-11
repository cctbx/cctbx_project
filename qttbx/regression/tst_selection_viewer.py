from __future__ import absolute_import, division, print_function

import json

from qttbx.sel_convert_chimera import translate_phenix_selection_string



def test_known_successes():
  # test a list of (phenix,chimera) pairs that should succeed
  with open("selection_pairs_test_success.json","r") as fh:
    selection_pairs = json.load(fh)
  #print("Testing selection pairs:",len(selection_pairs))
  for selection_phenix,selection_chimera in selection_pairs:
    selection_test = translate_phenix_selection_string(selection_phenix)
    if selection_test!=selection_chimera:
      print("Selection mismatch: ")
      print("Phenix:",selection_phenix)
      print("Expected Chimera:",selection_chimera)
      print("Translated Chimera: ",selection_test)
      assert False, "Selection mismatch"



if __name__=='__main__':
  test_known_successes()
  print("OK")
