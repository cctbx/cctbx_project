from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.phil
import mmtbx.polygon
import libtbx, os, re, sys
from libtbx.utils import Sorry
from libtbx import easy_pickle
from six.moves import range

keys_to_show = ["r_work", "r_free",
  "pdb_header_r_work", "pdb_header_r_free",
  "r_work_cutoffs", "r_free_cutoffs",
  "completeness_in_range", "completeness_d_min_inf", "completeness_6A_inf",
  "adp_mean_all", "adp_min_all", "adp_max_all",
  "wilson_b", "solvent_content_via_mask",
  "bond_rmsd", "bond_max_deviation", "angle_rmsd", "angle_max_deviation",
  "dihedral_rmsd", "dihedral_max_deviation",
  "planarity_rmsd", "planarity_max_deviation",
  "chirality_rmsd", "chirality_max_deviation",
  "rama_favored", "rama_allowed", "rama_outliers",
  "rotamer_outliers", "clashscore"]

other_numerical_keys = ["high_resolution", "low_resolution",
  'adp_min_backbone', 'adp_min_sidechain', 'adp_min_solvent',
  'adp_max_backbone', 'adp_max_sidechain', 'adp_max_solvent',
  'adp_mean_backbone', 'adp_mean_sidechain', 'adp_mean_solvent',
  'unit_cell_volume']

default_keys = ["r_work", "r_free", "adp_mean_all", "bond_rmsd", "angle_rmsd",
  "clashscore"]

key_captions = ["R-work", "R-free", "R-work (PDB)", "R-free (PDB)",
  "R-work (after cutoff)", "R-free (after cutoff)",
  "Completeness in range", "Completeness", "Completeness to 6A",
  "Average B", "Minimum B", "Maximum B",
  "Wilson B", "Solvent content",
  "RMSD(bonds)", "Bonds max.", "RMSD(angles)", "Angles max.",
  "RMSD(dihedrals)", "Dihedrals max.", "RMSD(planarity)", "Planarity max",
  "RMSD(chirality)", "Chirality max.",
  "Ramachandran favored", "Ramachandran allowed", "Ramachandran outliers",
  "Rotamer outliers", "Clashscore"]

other_captions = ["High resolution", "Low resolution",
  "Min. backbone ADP", "Min. sidechain ADP", "Min. solvent ADP",
  "Max. backbone ADP", "Max. sidechain ADP", "Max. solvent ADP",
  "Mean backbone ADP", "Mean sidechain ADP", "Mean solvent ADP",
  "Unit cell volume"]

assert len(keys_to_show) == len(key_captions)
_selected = []
for key_name in keys_to_show :
  if key_name in default_keys :
    _selected.append("*%s" % key_name)
  else :
    _selected.append(key_name)
key_params_str = " ".join(_selected)
captions_str = " ".join([ re.sub(" ", "_", txt) for txt in key_captions ])

polygon_params_str = """\
  database_file_name = None
    .type = str
    .style = noauto
  keys_to_show = %s
    .type = choice(multi=True)
    .short_caption = Statistics to display
    .caption = %s
    .style = bold hide_label
  number_of_histogram_slots = 10
    .type = int
    .help = Number of histogram slots for the final histogram to be used to \
            draw the POLYGON's rays.
    .input_size = 64
    .style = noauto bold
""" % (key_params_str, captions_str)

all_params_str = """
polygon {
  %s
}""" % polygon_params_str

master_params = iotbx.phil.parse(all_params_str)

def select_dict(database_dict, selection):
  result = {}
  for key in database_dict.keys():
    result.setdefault(key, database_dict[key].select(selection))
  return result

def filter_and_convert(database_dict, keys):
  selection = flex.bool(database_dict[keys[0]].size(), True)
  for key in keys+["high_resolution"]:
    values = database_dict[key]
    selection &= (values != "none")
  tmp = select_dict(database_dict = database_dict, selection = selection)
  result = {}
  for key in keys+["high_resolution"]:
    vals = flex.double([float(v) for v in tmp[key]])
    result.setdefault(key, vals)
  return result

def show_histogram(data, n_slots, smooth = True):
  triplets = []
  histogram = flex.histogram(data = data, n_slots = n_slots)
  l = histogram.data_min()
  for i, s in enumerate(histogram.slots()):
    r = histogram.data_min() + histogram.slot_width() * (i+1)
    triplets.append( [l, r, s] )
    print("%8.4f %8.4f %d" % (l, r, s))
    l = r
  if(smooth):
    print("... smooth histogram")
    triplets_smooth = []
    for i, t in enumerate(triplets):
      values = flex.double()
      for j in [-1,0,1]:
        if(i+j >=0 and i+j < len(triplets)):
          values.append(float(triplets[i+j][2]))
      triplets_smooth.append((t[0],t[1],flex.mean(values)))
    for t in triplets_smooth:
      print("%8.4f %8.4f %d" % (t[0], t[1], int("%.0f"%t[2])))
  return histogram

def convert_to_histogram(data, n_slots):
  histogram = flex.histogram(data=data, n_slots=n_slots)
  return histogram

def apply_default_filter(database_dict, d_min, key = "high_resolution"):
  d_mins = database_dict["high_resolution"]
  offset = 0.1
  if(d_min>=3 and d_min<4): offset = 0.2
  if(d_min>=4 and d_min<6): offset = 0.5
  if(d_min>=6):             offset = 1.0
  sel  = (d_mins>(d_min-offset))
  sel &= (d_mins<(d_min+offset))
  result = select_dict(database_dict = database_dict, selection = sel)
  # Totally ad-hoc manipulation for histograms to make sense and format nicely.
  # Perhaps needs to be revised at some point.
  sel = flex.bool(sel.count(True), True)
  for key in result.keys():
    if(key in ["high_resolution"]): continue
    vals = result[key]
    if vals.size()==0: continue
    if(key == "bond_rmsd"):
      sel &= vals < 0.05
    elif(key == "angle_rmsd"):
      sel &= vals < 5.
    else:
      mean = flex.mean(vals)
      sel &= vals > mean/2
      sel &= vals < mean*2
      if(key == "r_work" or key == "r_free"):
        sel &= vals < 0.45
  result = select_dict(database_dict=result, selection=sel)
  #
  return result

def load_db(file_name=None):
  if(file_name is None):
    file_name = libtbx.env.find_in_repositories(
      relative_path = "chem_data/polygon_data/all_mvd.pickle",
      test = os.path.isfile)
  assert os.path.isfile(file_name)
  database_dict = easy_pickle.load(file_name)

  # Python 3 pickle fix
  # =========================================================================
  if sys.version_info.major == 3:
    database_dict = easy_pickle.fix_py2_pickle(database_dict)
  # =========================================================================

  return database_dict

def polygon(params = master_params.extract(), d_min = None,
            show_histograms = True, extract_gui_data=False):
  database_dict = load_db(file_name=params.polygon.database_file_name)
  result = filter_and_convert(
    database_dict = database_dict,
    keys          = params.polygon.keys_to_show)
  if(d_min is not None):
    result = apply_default_filter(database_dict = result, d_min = d_min)
  histograms = []
  if extract_gui_data :
    for selected_key in params.polygon.keys_to_show:
      data = result[selected_key]
      histograms.append([selected_key, data]) # XXX: not really histograms!
  elif(show_histograms):
    for selected_key in params.polygon.keys_to_show:
      data = result[selected_key]
      print("%s data_points=%d" % (selected_key, data.size()), \
        "min/max/mean= %12.4f %12.4f %12.4f"%data.min_max_mean().as_tuple())
      n_slots = params.polygon.number_of_histogram_slots
      if(n_slots is None):
        assert 0
        n_slots = data.size()//50
        if(n_slots < 5):
          for scale in range(25,10,-1):
            n_slots = data.size()//scale
            if(n_slots >= 10): break
      if(n_slots == 0):
        raise Sorry("Not enough data selected.")
      h = show_histogram(data = data, n_slots = n_slots)
      histograms.append([selected_key,h])
  return histograms

def get_statistics_percentiles(d_min, stats):
  """
  For a given set of statistics, determine their percentile ranking compared
  to other crystal structures at similar resolution.
  """
  if (d_min is None):
    return dict([ (s, None) for s in stats.keys()  ])
  try :
    db = load_db()
  except Exception as e :
    return {}
  d_min_mvd = flex.double([ float(x) for x in db['high_resolution'] ])
  sel_perm = flex.sort_permutation(d_min_mvd)
  d_min_mvd = d_min_mvd.select(sel_perm)
  def find_value_in_list(values, value):
    i = 0
    j = len(values) - 1
    while (i != j):
      k = i + (j - i) // 2
      if (value and value <= values[k]):
        j = k
      else :
        i = k + 1
    return i
  index = find_value_in_list(d_min_mvd, d_min)
  sel_around = flex.bool(d_min_mvd.size(), False)
  index_tmp = index
  while (index_tmp > 0):
    d_min_other = d_min_mvd[index_tmp]
    if (d_min_other < d_min - 0.1):
      break
    sel_around[index_tmp] = True
    index_tmp -= 1
  index_tmp = index
  while (index_tmp < d_min_mvd.size()):
    d_min_other = d_min_mvd[index_tmp]
    if (d_min_other > d_min + 0.1):
      break
    sel_around[index_tmp] = True
    index_tmp += 1
  #print "%d structures around %g" % (sel_around.count(True), d_min)
  percentiles = {}
  for stat_name in stats.keys():
    stat = stats[stat_name]
    if (not stat_name in db):
      percentiles[stat_name] = None
      continue
    values = db[stat_name].select(sel_perm).select(sel_around)
    fvalues = flex.double()
    for value in values :
      try :
        fvalues.append(float(value))
      except ValueError :
        pass
    assert fvalues.size() != 0
    fvalues_sorted = fvalues.select(flex.sort_permutation(fvalues))
    stat_index = find_value_in_list(fvalues_sorted, stat)
    # FIXME I think for some of these statistics we need to reverse this -
    # i.e. if higher statistics are better
    stat_perc = 100 * (1 - (stat_index / fvalues.size()))
    percentiles[stat_name] = stat_perc
    #print stat_name, stat_index, fvalues.size(), stat_perc
    #flex.histogram(fvalues, n_slots=10).show(prefix="  ")
  return percentiles
