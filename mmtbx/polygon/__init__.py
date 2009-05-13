from scitbx.array_family import flex
import iotbx.phil
import mmtbx.polygon
import libtbx, os
from libtbx.utils import Sorry
from libtbx import easy_pickle


master_params = iotbx.phil.parse("""\
polygon {
  database_file_name = None
    .type = str
  keys_to_show = rama_favored *adp_mean cmpl_in_range cmpl_d_min_inf \
                 rama_allowed *r_work_pdb rama_outliers *k_sol *chirality_rmsd \
                 rama_general *angles_rmsd *r_free_cutoff r_free_re_computed \
                 *dihedrals_rmsd cmpl_6A_inf *r_free_pdb *b_sol matthews_coeff \
                 solvent_cont wilson_b r_work_re_computed angles_max \
                 dihedrals_max adp_min chirality_max *planarity_rmsd adp_max \
                 bonds_max *bonds_rmsd planarity_max *r_work_cutoff

    .type = choice(multi=True)
  number_of_histogram_slots = None
    .type = int
    .help = Number of histogram slots for the final histrogram to be used to \
            draw the POLYGON's rays.
  max_reject_fraction = 0.1
    .type = float
    .help = Fraction of models allowed to be rejected as outliers.
  max_models_for_default_filter = 1000
    .type = int
  filter
    .multiple = True
    .help = Selection keys.
  {
    key = twinned rama_proline n_fobs_outl rama_favored adp_mean rna_dna \
          cmpl_in_range cmpl_d_min_inf n_npd unit_cell n_aniso space_group \
          rama_prepro non_bonded_min number_of_models rama_allowed r_work_pdb \
          n_refl rama_outliers k_sol occ_min occ_max chirality_rmsd \
          sigma_cutoff rama_general name angles_rmsd d_max_pdb r_free_cutoff \
          r_free_re_computed test_flag_value test_set_size dihedrals_rmsd \
          atom_counts d_max cmpl_6A_inf r_free_pdb year b_sol matthews_coeff \
          anom_flag n_refl_cutoff small_molecule d_min_pdb solvent_cont \
          wilson_b r_work_re_computed n_tls_groups angles_max dihedrals_max \
          adp_min d_min chirality_max planarity_rmsd adp_max n_altloc \
          data_label rama_glycine bonds_max bonds_rmsd planarity_max water \
          program_name occ_mean tls unit_cell_volume amino_acid element \
          n_atoms r_work_cutoff b_cart other
      .type = choice(multi=False)
    value_min = None
      .type = float
    value_max = None
      .type = float
    target_value= None
      .type = str
  }
}
""")

def select_database_dict_by_keys(select_keys, database_dict):
  result = {}
  for key in select_keys:
    result.setdefault(key, database_dict[key])
  return result

def convert_to_numeric(values):
  is_digit_selection = flex.bool()
  for kv in values:
    flag = True
    try:
      tmp = float(kv)
    except: flag = False
    is_digit_selection.append(flag)
  if(is_digit_selection.count(False) > 0):
    raise RuntimeError("Non-numerical values.")
  points = 0
  for kv in values:
    points += kv.count(".")
  if(points == 0):
    new_values = flex.int()
    for kv in values:
      new_values.append(int(kv))
  else:
    new_values = flex.double()
    for kv in values:
      new_values.append(float(kv))
  return new_values

def select_dict(database_dict, selection):
  result = {}
  for key in database_dict.keys():
    result.setdefault(key, database_dict[key].select(selection))
  return result

def order_by_value(database_dict, key, reverse = False):
  values = database_dict[key]
  new_values = convert_to_numeric(values = values)
  selection = flex.sort_permutation(new_values, reverse = reverse)
  return select_dict(database_dict = database_dict, selection = selection)

def leave_all_available_entries(database_dict, keys):
  selection = flex.bool(database_dict[keys[0]].size(), True)
  for key in keys:
    values = database_dict[key]
    selection &= values != "none"
  return select_dict(database_dict = database_dict, selection = selection)

def filter_database(database_dict, key, value_min = None, value_max = None,
                    target_value = None):
  values = database_dict[key]
  if(target_value is not None):
    assert [value_min,value_max].count(None) == 2
    selection = flex.bool(values.size(), False)
    for i, v in enumerate(values):
      if(v.count(target_value)>0): selection[i] = True

  else:
    assert [value_min,value_max].count(None) < 2
    if([value_min,value_max].count(None) == 0): assert value_min < value_max
    new_values = convert_to_numeric(values = values)
    selection = flex.bool(new_values.size(), True)
    if(value_min is not None):
      if(isinstance(new_values, flex.int)): value_min = int(value_min)
      selection &= new_values > value_min
    if(value_max is not None):
      if(isinstance(new_values, flex.int)): value_max = int(value_max)
      selection &= new_values < value_max
  return select_dict(database_dict = database_dict, selection = selection)

def filter_histogram_of_key_value(database_dict, key, max_reject_fraction,
                                  edge_tolerance_small = 1.e-4, n_slots = 3):
  values = database_dict[key]
  new_values = convert_to_numeric(values = values)
  #print flex.min(new_values), flex.max(new_values)
  size = new_values.size()
  #print size
  if(size == 0): return
  while True:
    values = database_dict[key]
    new_values = convert_to_numeric(values = values)
    selection = flex.bool(new_values.size(), True)
    histogram = flex.histogram(data = new_values, n_slots = n_slots)
    l = histogram.data_min()
    for i, s in enumerate(histogram.slots()):
      r = histogram.data_min() + histogram.slot_width() * (i+1)
      r = r+edge_tolerance_small
      l = max(0, l-edge_tolerance_small)
      #print "%8.4f %8.4f %d" % (l, r, s)
      if(s < size * max_reject_fraction):
         selection &= ~((new_values >= l) & (new_values <= r))
      l = r
    #print
    leave, remove = selection.count(True), selection.count(False)
    #print leave, remove
    if(remove == 0): break
    if(size - leave > int(size * max_reject_fraction)): break
    database_dict = select_dict(database_dict = database_dict,
                                selection     = selection)
  return database_dict

def show_histogram(data, n_slots, smooth = True):
  triplets = []
  histogram = flex.histogram(data = data, n_slots = n_slots)
  l = histogram.data_min()
  for i, s in enumerate(histogram.slots()):
    r = histogram.data_min() + histogram.slot_width() * (i+1)
    triplets.append( [l, r, s] )
    print "%8.4f %8.4f %d" % (l, r, s)
    l = r
  if(smooth):
    print "... smooth histogram"
    triplets_smooth = []
    for i, t in enumerate(triplets):
      values = flex.double()
      for j in [-1,0,1]:
        if(i+j >=0 and i+j < len(triplets)):
          values.append(float(triplets[i+j][2]))
      triplets_smooth.append((t[0],t[1],flex.mean(values)))
    for t in triplets_smooth:
      print "%8.4f %8.4f %d" % (t[0], t[1], int("%.0f"%t[2]))
  return histogram

def convert_to_histogram(data, n_slots) :
  histogram = flex.histogram(data=data, n_slots=n_slots)
  return histogram

def apply_default_filter(database_dict, d_min, max_models_for_default_filter,
                         key = "d_min"):
  database_dict = order_by_value(database_dict = database_dict, key = key)
  values = flex.double()
  for v in database_dict[key]: values.append(float(v))
  diff = flex.abs(values-d_min)
  min_val = flex.min(diff)
  i_min_sel = (diff == min_val).iselection()
  assert i_min_sel.size() > 0
  i_min = i_min_sel[i_min_sel.size()/2]
  i_l = max(0, i_min-max_models_for_default_filter/2)
  i_r = min(values.size()-1, i_min+max_models_for_default_filter/2)
  #
  print "apply_default_filter:"
  print "  found data points dmin->higher =", abs(i_l-i_min)
  print "  found data points dmin->lower  =", abs(i_r-i_min)
  imm = min(abs(i_l-i_min), abs(i_r-i_min))
  i_l, i_r = i_min-imm, i_min+imm
  print "  used data points dmin->higher =", imm
  print "  used data points dmin->lower  =", imm
  #
  selection = flex.bool(values.size(), False)
  for i in xrange(i_l,i_r): selection[i] = True
  return select_dict(database_dict = database_dict, selection = selection)

def polygon(params = master_params.extract(), d_min = None,
            show_histograms = True):
  if(params.polygon.database_file_name is None):
    file_name = libtbx.env.find_in_repositories(
      relative_path = "chem_data/polygon_data/phenix_mvd_2009_APR_14_21h23.pickle",
      test = os.path.isfile)
  else:
    file_name = params.polygon.database_file_name
    assert os.path.isfile(file_name)
  database_dict = easy_pickle.load(file_name)
  result = leave_all_available_entries(
    database_dict = database_dict,
    keys          = params.polygon.keys_to_show)
  filters = params.polygon.filter
  if((d_min is not None and len(filters) == 1 and filters[0].key is None) or
     len(filters) == 0):
    result = apply_default_filter(database_dict = result, d_min = d_min,
      max_models_for_default_filter =
        params.polygon.max_models_for_default_filter)
  else:
    for filter in filters:
      result = filter_database(
        database_dict = result,
        key           = filter.key,
        value_min     = filter.value_min,
        value_max     = filter.value_max,
        target_value  = filter.target_value)
  result = select_database_dict_by_keys(
    select_keys   = params.polygon.keys_to_show,
    database_dict = result)
  for key_to_show in params.polygon.keys_to_show:
    result = filter_histogram_of_key_value(
      database_dict       = result,
      key                 = key_to_show,
      max_reject_fraction = params.polygon.max_reject_fraction)
  histograms = []
  if(show_histograms):
    for selected_key in params.polygon.keys_to_show:
      data = convert_to_numeric(values=result[selected_key])
      print "%s data_points=%d" % (selected_key, data.size()), \
        "min/max/mean= %12.4f %12.4f %12.4f"%data.min_max_mean().as_tuple()
      n_slots = params.polygon.number_of_histogram_slots
      if(n_slots is None):
        n_slots = data.size()/50
        if(n_slots < 5):
          for scale in range(25,10,-1):
            n_slots = data.size()/scale
            if(n_slots >= 10): break
      if(n_slots == 0):
        raise Sorry("Not enough data selected.")
      h = show_histogram(data = data, n_slots = n_slots)
      histograms.append([selected_key,h])
  return histograms
