from __future__ import absolute_import, division, print_function

import os
import subprocess
import hashlib
import numpy as np
from pdb_strs import (
    pdb_str_normal,
    pdb_str_normal_misordered,
    pdb_str_normal_missing_res,
    pdb_str_normal_missing_atoms,
    pdb_str_icode,
    pdb_str_icode_reorder,
    pdb_str_alts,
    pdb_str_alts_misorderd,
    pdb_str_non_numeric,
    pdb_str_non_ascending
)

from iotbx.data_manager import DataManager
from qttbx.viewers.chimerax import ChimeraXViewer
from qttbx.sel_convert_chimera import translate_phenix_selection_string

"""
This file defines a "test" as a list or tuple of: 

  pdb_str_key, operator, keyword ,start,stop = test

For example: 
  ("pdb_str_normal", "through", "resid", 3, 5)

Where "pdb_str_normal" is the string key referring to a pdb string/model object used for testing.
This is just a structured alternative to something like:
  
  model = pdb_models[pdb_str_key]
  model.selection("resid 3 through 5")

"""



#####################################
# Set up test models
#####################################
pdb_strs = {'pdb_str_normal':pdb_str_normal,
            'pdb_str_normal_misordered':pdb_str_normal_misordered,
 'pdb_str_normal_missing_res':pdb_str_normal_missing_res,
 'pdb_str_normal_missing_atoms':pdb_str_normal_missing_atoms,
 'pdb_str_icode':pdb_str_icode,
 'pdb_str_icode_reorder':pdb_str_icode_reorder,
 'pdb_str_alts':pdb_str_alts,
 'pdb_str_alts_misorderd':pdb_str_alts_misorderd,
 'pdb_str_non_numeric':pdb_str_non_numeric,
 'pdb_str_non_ascending':pdb_str_non_ascending
           }

pdb_models = {}
for key,str_ in pdb_strs.items():
  dm = DataManager()
  dm.process_model_str(key,str_)
  model = dm.get_model()
  model.add_crystal_symmetry_if_necessary()
  pdb_models[key] = model


viewer = ChimeraXViewer()

#####################################
# Utility functions
#####################################


def eval_test(test1):
   # evaluate a test. 
    try:
        
        pdb_str_key1, operator1, keyword1 ,start1,stop1 = test1
        
        pdb_str1 = pdb_strs[pdb_str_key1]
        
        model1 = pdb_models[pdb_str_key1]
        sel_str1 =f"{keyword1} {start1} {operator1} {stop1}"
        selection1 = model1.selection(sel_str1)
        result1 = model1.select(selection1).model_as_pdb()
        return True,list(selection1.as_numpy_array().nonzero()[0])
    except:
        return False,None



def compare_tests(t1,t2,sel_str1=None,sel_str2=None):
    # compare two tests to verify they return the same content
    # First checks string, and if that is not equal, checks cartesian coordinates
    pdb_str_key1, operator1, keyword1 ,start1,stop1 = t1
    pdb_str_key2, operator2, keyword2 ,start2,stop2 = t2
    model1 = pdb_models[pdb_str_key1]


    # test 1
    if sel_str1 is None:
        sel_str1 =f"{keyword1} {start1} {operator1} {stop1}"
    try:
        selection1 = model1.selection(sel_str1)
        result1 = model1.select(selection1)
    except:
        result1 = None


    model2 = pdb_models[pdb_str_key2]
    try:
        if sel_str2 is None:
            sel_str2 =f"{keyword2} {start2} {operator2} {stop2}"
        selection2 = model2.selection(sel_str2)
        result2 = model2.select(selection2)
    except:
        result2 = None
    if result1.get_number_of_atoms() != result2.get_number_of_atoms():
        same = False
    else:
        same =np.all(np.isclose(result1.get_sites_cart().as_numpy_array(),
                            result2.get_sites_cart().as_numpy_array()))
        if not same:
            s1 = set([frozenset([a,b,c]) for a,b,c in result1.get_sites_cart().as_numpy_array()])
            s2 = set([frozenset([a,b,c]) for a,b,c in result2.get_sites_cart().as_numpy_array()])
            same = s1==s2
    return same#,result1.model_as_pdb(),result2.model_as_pdb()

def eval_chimerax_conversion(test):
    """
    Use subprocess to run ChimeraX. For a test input:

    1. A test (cctbx) is translated to ChimeraX
    2. The ChimeraX selection is used to save a new pdb file with only selected atoms.
    3. The new file is opened in cctbx
    4. Verification of identical atoms selected in cctbx/ChimeraX


    Note: Tests which select nothing in cctbx return None, since ChimeraX will not open a file without any atoms.
    
    """  
  
  
    pdb_str_key, operator, keyword ,start,stop = test
    model = pdb_models[pdb_str_key]
    sel_str_cctbx =f"{keyword} {start} {operator} {stop}"
    sel_string_cctbx = sel_str_cctbx
    sel_str_chimera = translate_phenix_selection_string(sel_str_cctbx)


    
    test_folder = os.getcwd()+"/chimera_tests"
    filename = os.getcwd()+"/../devel/devel/"+pdb_str_key+".pdb"
    if not os.path.exists(filename):
      dm = DataManager()
      dm.write_model_file(model,filename=filename)
    file = filename
    if not os.path.exists(test_folder):
        os.mkdir(test_folder)
    test_file = test_folder+"/chimerax_selection_test_"+str(hash((filename,sel_str_cctbx)))+".pdb"

    
    
    try:
        # load the model and make selection
        
        selection = model.selection(sel_string_cctbx)
        #print(selection.count(True))
        if selection.count(True)==0:
            return None
        model_sel = model.select(selection)
        #print(model_sel.get_number_of_atoms())
        #print(model_sel.model_as_pdb())
        #print("cctbx atoms: ",model1_sel.get_number_of_atoms())
        coord_set1 = set([frozenset(coord) for coord in model_sel.get_sites_cart().as_numpy_array()])
        
        
        chimerax_commands = f"""
        open {file}
        sel {sel_str_chimera}
        save {test_file} sel true
        exit
        """
        # Write the commands to a temporary ChimeraX script file.
        with open('temp.cxc', 'w') as f:
            f.write(chimerax_commands)
        
        # Use subprocess to run ChimeraX with the script file.
        _ = subprocess.run([viewer.find_viewer(), "--nogui", "--silent","temp.cxc"])
        
        # Remember to clean up the temporary script file when you're done.
        os.remove('temp.cxc')
        

        # reopen as text to check for atoms
        with open(test_file,"r") as fh:
            string = fh.read()

        if "ATOM" not in string:
            return None
        # reopen with cctbx
        dm2 = DataManager()
        dm2.process_model_file(test_file)
        os.remove(test_file)
        model2 = dm2.get_model()
        #print(model2.model_as_pdb())
        #print("chimerax atoms: ",model2.get_number_of_atoms())
        model2.add_crystal_symmetry_if_necessary()
        coord_set2 = set([frozenset(coord) for coord in model2.get_sites_cart().as_numpy_array()])
        #assert coord_set1==coord_set2, "Translation from cctbx to chimerax yields two different sets of cartesian coordinates"
        passed = coord_set1==coord_set2
        #print("\t","Passed:",passed)
        return passed
    except:
        raise
        return False


#####################################
# Tests
#####################################
  
def tst_normal(chimera=False):
  # Using a "normal" file
  
  # Resseq
  
  #   Colon
  
  # resseq,colon takes all the alts
  t1 = ("pdb_str_normal",":","resseq",1,3)
  works, int_sel = eval_test(t1)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  if chimera:
      assert eval_chimerax_conversion(t1)
  
  # will cross chain boundry
  t2 = ("pdb_str_normal",":","resseq",1,8)
  works, int_sel = eval_test(t2)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  # Resid
  
  # Colon
  
  # resid,colon works the same as resseq
  t3 = ("pdb_str_normal",":","resid",1,3)
  works, int_sel = eval_test(t3)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t1,t3)
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  # Will cross chain boundry
  t4 = ("pdb_str_normal",":","resid",1,8)
  works,int_sel = eval_test(t4)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]
  assert compare_tests(t2,t4)
  if chimera:
      assert eval_chimerax_conversion(t4)
  
  #through
  
  # resid will take all alts
  t5 = ("pdb_str_normal","through","resid",1,3)
  works,int_sel = eval_test(t5)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t3,t5)
  if chimera:
      assert eval_chimerax_conversion(t5)
  
  # take all alts
  t6 = ("pdb_str_normal","through","resid",1,3)
  works_int_sel = eval_test(t6)
  assert int_sel ==[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t3,t6)
  if chimera:
      assert eval_chimerax_conversion(t6)
  
  # **will NOT cross chain boundry
  t7 = ("pdb_str_normal","through","resid",1,8)
  works, int_sel = eval_test(t7)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
  assert not compare_tests(t4,t7)
  if chimera:
      assert not eval_chimerax_conversion(t7)


def tst_missing_residues(chimera=False):
  # Using a file missing residues (missing 3)

  # Resseq
  
  #   Colon
  
  # resseq,colon takes the residues in numerical range only, just omits the missing ones
  t1 = ("pdb_str_normal_missing_res",":","resseq",1,3)
  works, int_sel = eval_test(t1)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
  if chimera:
      assert eval_chimerax_conversion(t1)
  
  # same
  t2 = ("pdb_str_normal_missing_res",":","resseq",1,4)
  works, int_sel = eval_test(t2)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  # Resid
  
  # Colon
  
  # resid,colon works the same as resseq
  t3 = ("pdb_str_normal_missing_res",":","resid",1,3)
  works,int_sel = eval_test(t3)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
  assert compare_tests(t1,t3)
  if chimera:
      assert eval_chimerax_conversion(t3)
  
  t4 = ("pdb_str_normal_missing_res",":","resid",1,4)
  works,int_sel = eval_test(t4)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
  assert compare_tests(t2,t4)
  if chimera:
      assert eval_chimerax_conversion(t4)
  
  # through
  
  # ** ending on a missing residue will extend to end of chain
  t5 = ("pdb_str_normal_missing_res","through","resid",1,3)
  works, int_sel = eval_test(t5)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
  if chimera:
      assert not eval_chimerax_conversion(t5)
  
  # ** starting on missing residue returns nothing
  t6 = ("pdb_str_normal_missing_res","through","resid",3,4)
  works,int_sel = eval_test(t6)
  assert int_sel == []
  if chimera:
      assert eval_chimerax_conversion(t6) is None

def tst_non_ascending(chimera=False):
  # Using a file with integer but inserting a high non-ascending seq id (777)
  
  # Resseq
  
  #   Colon
  
  # resseq,colon takes the residues in numerical range only
  t1 = ("pdb_str_non_ascending",":","resseq",1,3)
  works, int_sel = eval_test(t1)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22]
  if chimera:
      assert eval_chimerax_conversion(t1)
  
  # takes the full numerical range
  t2 = ("pdb_str_non_ascending",":","resseq",1,777)
  works, int_sel = eval_test(t2)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  # Resid
  
  # Colon
  
  # resid,colon works the same as resseq
  t3 = ("pdb_str_non_ascending",":","resid",1,3)
  works,int_sel = eval_test(t3)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t1,t3)
  if chimera:
      assert eval_chimerax_conversion(t3)
  
  # This will take anything in the full numerical range
  t4 = ("pdb_str_non_ascending",":","resid",1,777)
  works, int_sel = eval_test(t4)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]
  assert compare_tests(t2,t4)
  if chimera:
      assert eval_chimerax_conversion(t4)
  
  # through
  
  # resid with through will take a high range value if it is within file range
  t5 = ("pdb_str_non_ascending","through","resid",1,3)
  works, int_sel = eval_test(t5)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  if chimera:
      assert not eval_chimerax_conversion(t5)
  
  
  # specifying the value will terminate early, not include full numerical range
  t6 = ("pdb_str_non_ascending","through","resid",1,777)
  works, int_sel = eval_test(t6)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
  if chimera:
      assert not eval_chimerax_conversion(t6)
  
  # going backwards doesn't work, 1 is treated as missing
  t7 = ("pdb_str_non_ascending","through","resid",777,1)
  works, int_sel = eval_test(t7)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
  if chimera:
      assert not eval_chimerax_conversion(t7)

    
def tst_non_numeric(chimera=False):
    
  # Using non-numeric sequence ids (like “A01”) 
  
  
  # Resseq
  # Cannot parse as a specifier with resseq, colon
  t1 = ("pdb_str_non_numeric",":","resseq","A01",7) # expect eval fail
  works, int_sel = eval_test(t1)
  assert not works
  
  # If within range it parses, but they are also not selected.
  t2 = ("pdb_str_non_numeric",":","resseq",1,7) # expect all numeric within range to be selected
  works, int_sel = eval_test(t2)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 29, 30, 31, 32, 33]
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  #Resid
  
  # colon
  # Cannot parse as a specifier with resseq, colon
  t3 = ("pdb_str_non_numeric",":","resseq","A01",7) # expect eval fail
  works, int_sel = eval_test(t3)
  assert not works
  
  # (Same as resseq) If within range it parses, but they are also not selected.
  t4 = ("pdb_str_non_numeric",":","resid",1,7) # expect all numeric within range to be selected
  works, int_sel = eval_test(t4)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 29, 30, 31, 32, 33]
  assert compare_tests(t2,t4)
  if chimera:
      eval_chimerax_conversion(t4)
  
  
  # through
  # Cannot parse as a specifier with resid, through
  t5 = ("pdb_str_non_numeric","through","resseq","A01",7)
  works, int_sel = eval_test(t5)
  assert not works
  
  # **Will select these values with resid,through due to them being in the file-based range
  t6 = ("pdb_str_non_numeric","through","resid",1,7)
  works, int_sel = eval_test(t6)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28]
  if chimera:
      assert not eval_chimerax_conversion(t6)


def tst_insertion_codes(chimera=False):
  # Using a file with insertion codes
  
  
  # Resseq
  
  #   Colon
  
  # resseq,colon takes the icodes if in range
  t1 = ("pdb_str_icode",":","resseq",1,3)
  works, int_sel = eval_test(t1)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  if chimera:
      assert eval_chimerax_conversion(t1)
  
  # ending on insertion integer takes all insertions
  t2 =  ("pdb_str_icode",":","resseq",1,2)
  works, int_sel = eval_test(t2)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  # trying to specify icode takes to end of structure. Note, this is Note the same behavior as "A01" which fails
  t3 =  ("pdb_str_icode",":","resseq",1,"2A")
  works, int_sel = eval_test(t3)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
  if chimera:
      assert not eval_chimerax_conversion(t3)
  
  #  starting with the insertion integer takes the icodes too
  t4 =  ("pdb_str_icode",":","resseq",2,5)
  works, int_sel = eval_test(t4)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  if chimera:
      assert eval_chimerax_conversion(t4)
  
  # starting with the icode notation fails to parse
  t5 = ("pdb_str_icode",":","resseq","2A",5)
  works, int_sel = eval_test(t5)
  assert not works
  
  # Resid
  
  # colon 
  
  # Same as resseq
  t6 = ("pdb_str_icode",":","resid",1,3)
  works, int_sel = eval_test(t6)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  assert compare_tests(t1,t6)
  if chimera:
      assert eval_chimerax_conversion(t6)
  
  # ending on insertion integer does NOT take insertions, which come after
  t7 =  ("pdb_str_icode",":","resseq",1,2)
  works, int_sel = eval_test(t7)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t2,t7)
  if chimera:
      assert eval_chimerax_conversion(t7)
  
  # **trying to specify icode works as expected, includes 2,2A
  t8 =  ("pdb_str_icode",":","resid",1,"2A")
  works, int_sel = eval_test(t8)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert not compare_tests(t3,t8)
  if chimera:
      assert eval_chimerax_conversion(t8)
  
  #  starting with the insertion integer takes the icodes too.
  t9 =  ("pdb_str_icode",":","resid",2,5)
  works, int_sel = eval_test(t9)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  assert compare_tests(t4,t9)
  if chimera:
      assert eval_chimerax_conversion(t9)
  
  # ** starting with the icode notation ONLY takes 2A - 5, not 2. 
  t10 = ("pdb_str_icode",":","resid","2A",5)
  works, int_sel = eval_test(t10)
  assert int_sel == [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  if chimera:
      assert eval_chimerax_conversion(t10)
  
  # Probably Chimera will need to restrict specifying icodes to single residue selections
  
  # through
  
  # Same as resseq
  t11 = ("pdb_str_icode","through","resid",1,3)
  works, int_sel = eval_test(t11)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  assert compare_tests(t6,t11)
  if chimera:
      assert eval_chimerax_conversion(t11)
  
  # ** same as colon, ending on insertion integer does NOT take insertions, which come after
  t12 =  ("pdb_str_icode","through","resid",1,2)
  works, int_sel = eval_test(t12)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
  assert not compare_tests(t7,t12)
  if chimera:
      assert eval_chimerax_conversion(t12)
  
  
  # trying to specify icode works as expected, includes 2,2A
  t13 =  ("pdb_str_icode","through","resid",1,"2A")
  works, int_sel = eval_test(t13)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t8,t13)
  if chimera:
      assert eval_chimerax_conversion(t13)
  
  #  starting with the insertion integer takes the icodes too.
  t14 =  ("pdb_str_icode","through","resid",2,5)
  works, int_sel = eval_test(t14)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  assert compare_tests(t9,t14)
  if chimera:
      assert eval_chimerax_conversion(t14)
  
  # starting with the icode notation ONLY takes 2A - 5, not 2. 
  t15 = ("pdb_str_icode","through","resid","2A",5)
  works, int_sel = eval_test(t15)
  assert int_sel == [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  assert compare_tests(t10,t15)
  if chimera:
      assert eval_chimerax_conversion(t15)
  
  # Here, through and colon are identical

def tst_misordered_insertions(chimera=False):
  # Using a file with misordered icodes

  
  # Resseq
  
  #   Colon
  
  # resseq,colon takes the icodes if in range
  t1 = ("pdb_str_icode_reorder",":","resseq",1,3)
  works, int_sel = eval_test(t1)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  if chimera:
      assert eval_chimerax_conversion(t1)
  
  # ending on insertion integer takes all insertions
  t2 =  ("pdb_str_icode_reorder",":","resseq",1,2)
  works, int_sel = eval_test(t2)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  # trying to specify icode takes to end of structure.
  t3 =  ("pdb_str_icode_reorder",":","resseq",1,"2A")
  works, int_sel = eval_test(t3)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44]
  if chimera:
      assert not eval_chimerax_conversion(t3)
  
  #  starting with the insertion integer takes the icodes too
  t4 =  ("pdb_str_icode_reorder",":","resseq",2,5)
  works, int_sel = eval_test(t4)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  if chimera:
      assert eval_chimerax_conversion(t4)
  
  # starting with the icode notation fails to parse
  t5 = ("pdb_str_icode_reorder",":","resseq","2A",5)
  works, int_sel = eval_test(t5)
  assert not works
  
  # Resid
  
  # colon 
  
  # Same as resseq
  t6 = ("pdb_str_icode_reorder",":","resid",1,3)
  works, int_sel = eval_test(t6)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  assert compare_tests(t1,t6)
  if chimera:
      assert eval_chimerax_conversion(t6)
  
  # ending on insertion integer does NOT take insertions, which come after
  t7 =  ("pdb_str_icode_reorder",":","resseq",1,2)
  works, int_sel = eval_test(t7)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t2,t7)
  if chimera:
      assert eval_chimerax_conversion(t7)
  
  # ** trying to specify icode works as expected, includes 2,2A
  t8 =  ("pdb_str_icode_reorder",":","resid",1,"2A")
  works, int_sel = eval_test(t8)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert not compare_tests(t3,t8)
  if chimera:
      assert eval_chimerax_conversion(t8)
  
  #  starting with the insertion integer takes the icodes too.
  t9 =  ("pdb_str_icode_reorder",":","resid",2,5)
  works, int_sel = eval_test(t9)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  assert compare_tests(t4,t9)
  if chimera:
      assert eval_chimerax_conversion(t9)
  
  # starting with the icode notation ONLY takes 2A - 5, not 2. An implicit ordering of insertions is followed
  t10 = ("pdb_str_icode_reorder",":","resid","2A",5)
  works, int_sel = eval_test(t10)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  if chimera:
      assert eval_chimerax_conversion(t10)
  
  # through
  
  # same as resseq
  t11 = ("pdb_str_icode_reorder","through","resid",1,3)
  works, int_sel = eval_test(t11)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  assert compare_tests(t6,t11)
  if chimera:
      assert eval_chimerax_conversion(t11)
  
  # same as colon, ending on insertion integer DOES take insertions, because inserted residues were before
  t12 =  ("pdb_str_icode_reorder","through","resid",1,2)
  works, int_sel = eval_test(t12)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  assert compare_tests(t7,t12)
  if chimera:
      assert eval_chimerax_conversion(t12)
  
  # **trying to specify icode ONLY takes 2A, because it stops in file order before 2
  t13 =  ("pdb_str_icode_reorder","through","resid",1,"2A")
  works, int_sel = eval_test(t13)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
  assert not compare_tests(t8,t13)
  if chimera:
      assert not eval_chimerax_conversion(t13)
  
  #  starting with the insertion integer does NOT take icodes, because they are earlier in file order
  t14 =  ("pdb_str_icode_reorder","through","resid",2,5)
  works, int_sel = eval_test(t14)
  assert int_sel == [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  assert not compare_tests(t9,t14)
  if chimera:
      assert eval_chimerax_conversion(t14)
      
  # starting with the icode notation DOES 2A,2 - 5. Again, file order 
  t15 = ("pdb_str_icode_reorder","through","resid","2A",5)
  works, int_sel = eval_test(t15)
  assert int_sel == [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  assert not compare_tests(t10,t15)
  if chimera:
      assert not eval_chimerax_conversion(t15)


def tst_alt_locs(chimera=False):
  # Using a file with alt locs (3 has A,B)

  # Resseq
  
  #   Colon
  
  # resseq,colon takes all the alts
  t1 = ("pdb_str_alts",":","resseq",1,3)
  works, int_sel = eval_test(t1)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  if chimera:
      assert eval_chimerax_conversion(t1)
  
  # takes the full numerical range
  t2 = ("pdb_str_alts",":","resseq",1,4)
  works, int_sel = eval_test(t2)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  if chimera:
      assert eval_chimerax_conversion(t2)
  
  # Resid
  
  # Colon
  
  # resid,colon works the same as resseq
  t3 = ("pdb_str_alts",":","resid",1,3)
  works, int_sel = eval_test(t3)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  assert compare_tests(t1,t3)
  if chimera:
      assert eval_chimerax_conversion(t3)
  
  # This will take anything in the full numerical range, all alts
  t4 = ("pdb_str_alts",":","resid",1,4)
  works, int_sel = eval_test(t4)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
  assert compare_tests(t2,t4)
  if chimera:
      assert eval_chimerax_conversion(t4)
  
  # through
  
  # resid will take all alts
  t5 = ("pdb_str_alts","through","resid",1,3)
  works, int_sel = eval_test(t5)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  assert compare_tests(t3,t5)
  if chimera:
      assert eval_chimerax_conversion(t5)
  
  # take all alts
  t6 = ("pdb_str_alts","through","resid",1,3)
  works, int_sel = eval_test(t6)
  assert int_sel == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
  assert compare_tests(t3,t6)
  if chimera:
      assert eval_chimerax_conversion(t6)
  
  # alts don't really have any changes among the tests. 
  # NOTE: The eval chimerax commands work via disk, but in the gui only 1 conformer will be visible


if __name__=='__main__':
  chimera = True
  tst_normal(chimera=chimera)
  tst_missing_residues(chimera=chimera)
  tst_non_ascending(chimera=chimera)
  tst_non_numeric(chimera=chimera)
  tst_insertion_codes(chimera=chimera)
  tst_misordered_insertions(chimera=chimera)
  tst_alt_locs(chimera=chimera)
  
  print("OK")
