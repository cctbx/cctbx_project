#from lastmol import LastMol # TODO: Remove
from atom_sites import AtomSites
from iotbx.data_manager import DataManager
import random
import numpy as np
try:
  from tqdm import tqdm
except ImportError:
  def tqdm(iterable, *args, **kwargs):
    return iterable

    
# A fundamental concept here is that integer selections are NEVER trustworthy between objects.
# One integer selection is likely to produce two different string selections for two objects.
# When testing integer ranges between objects, go:
#     int_sel0 -> str_sel1
#     int_sel0 -> str_sel2
#     df1 = obj1.select(str_sel1)
#     df2 = obj2.select(str_sel2)
#     


def generate_random_ranges(k, n_min, n_max, N_max):
  generated_ranges = []
  for _ in range(k):
    range_len = random.randint(n_min, n_max)
    range_start = random.randint(0, N_max - range_len)
    range_end = range_start + range_len
    generated_ranges.append((range_start, range_end))

  return generated_ranges

def generate_random_selections(N_max):
  n_selections = 100 # number of selections to return
  k_max = 10 # max number of ranges per selection
  k_min = 1  # min number of ranges per selection
  N_max = N_max # max range for any selection
  n_min = 1 # min length for any range
  n_max = 30 # max length for any range
  frac_remove = 0.33 # turn off these to make some ranges less continuous
  selection_ranges = []
  for i in range(n_selections):
    k = random.randint(k_min, k_max)
    ranges = generate_random_ranges(k,n_min,n_max,N_max)
    selection_ranges.append(ranges)
  selections = []
  for range_list in selection_ranges:
    selection_int = []
    for r in range_list:
      sel_int = list(range(*r))
      selection_int+=sel_int
    selection_int = np.array(sorted(list(set(selection_int))))
    selection_bool = np.full(N_max+1,False)
    selection_bool[selection_int] = True
    # remove some
    remove_ints = [random.randint(0,N_max) for i in range(int(N_max*frac_remove))]
    selection_bool[remove_ints] = False
    selections.append(selection_bool)
  return selections

def test_two_objects(obj1,obj2,sel_str):
  df1 = obj1.select_from_phenix_str(sel_str)
  df2 = obj2.select_from_phenix_str(sel_str)
  df1 = df1.drop(columns="id")
  df2 = df2.drop(columns="id")
  df1 = df1.sort_values(by=list(df1.columns)).reset_index(drop=True)
  df2 = df2.sort_values(by=list(df2.columns)).reset_index(drop=True)
  
  assert (df1==df2).all().all(), "A string selection returned different data for two objects"

# Start
# read from disk to df with hierarchy graph
mol = LastMol.from_file_cif("6cvm.cif")
df = mol.atom_sites

# read using iotbx
dm = DataManager()
dm.process_model_file("6cvm.cif")
model = dm.get_model()

# Build sites df without a mmtbx.model
sites = AtomSites(df)

# Build with an mmtbx.model
sites_model = AtomSites.from_mmtbx_model(model)

# test
selections = generate_random_selections(len(sites_model))
for sel in tqdm(selections):
  sel_bool = sel
  sel_int = np.where(sel)[0]

  # testing round trip int on obj1
  sel_str = sites._str_phenix_from_selection(sel_int)
  check_int1 = sites._select_from_phenix_str(sel_str).index.values
  assert np.all(sel_int==check_int1), "Round trip integer selection failed on obj1"

  # test round trip int on object2
  sel_str = sites_model._str_phenix_from_selection(sel_int)
  check_int2 = sites_model._select(sel_str,str_format="phenix").index.values
  assert np.all(sel_int==check_int2), "Round trip integer selection failed on obj2"

  # test round trip bool on object2
  sel_str = sites_model._str_phenix_from_selection(sel_bool)
  check_bool1 = np.full(len(sites_model),False)
  check_bool1[check_int1] = True
  assert np.all(sel_bool==check_bool1), "Round trip bool selection failed on obj2"

  # test round trip bool on object2
  sel_str = sites_model._str_phenix_from_selection(sel_bool)
  check_bool2 = np.full(len(sites_model),False)
  check_bool2[check_int2] = True
  assert np.all(sel_bool==check_bool2), "Round trip bool selection failed on obj2"


  # testing between two sites objects
  obj1 = sites
  obj2 = sites_model
  sel_str1 = sites._str_phenix_from_selection(sel_int)
  sel_str2 = sites_model._str_phenix_from_selection(sel_int)

  test_two_objects(obj1,obj2,sel_str1)
  test_two_objects(obj1,obj2,sel_str2)

print("OK")
