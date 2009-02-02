from scitbx.array_family import flex
from libtbx import group_args, smart_open
import iotbx.phil
import mmtbx.polygon

master_params = iotbx.phil.parse("""\
polygon {
  database_file_name = None
    .type = str
  keys_to_show = nseq npdb code dhigh dlow nrefl *rfact *rfree *bndav bndmx \
                 *angav angmx dihav dihmx plnav plnmx chrav chrmx nbdis adpmn \
                 *adpav adpmx occmn occav occmx natoms natiso nonpos nmod nsg \
                 acell bcell ccell alpha beta gamma data twin anom bwils \
                 *compl cmp-t cmp-l bulkk bulkb noutl progr year rhgpdb rlwpdb \
                 sigma t-set t-flag rpdb rfpdb vol/1000 mw(kda) cmatthw rmout \
                 rmfav rmgen mol:pno wat/res
    .type = choice(multi=True)
  number_of_histogram_slots = 10
    .type = int
    .help = Number of histogram slots for the final histrogram to be used to \
            draw the POLYGON's rays.
  max_reject_fraction = 0.1
    .type = float
    .help = Fraction of models allowed to be rejected as outliers.
  filter
    .multiple = True
    .help = Selection keys.
  {
    value_min = None
      .type = float
    value_max = None
      .type = float
    target_value= None
      .type = str
  }
}
""")

def extract_database_from_file(file_name): # XXX use pickle to speed up?
  database_lines = smart_open.for_reading(file_name =
    file_name).read().splitlines()
  keys_dict = {}
  for i_seq, key in enumerate(database_lines[0].split()):
    keys_dict.setdefault(i_seq, key)
    assert key.islower()
  number_of_keys = len(keys_dict.values())
  database_dict = {}
  for line in database_lines[1:]:
    line_split = line.split()
    assert len(line_split) == number_of_keys
    fs = flex.std_string
    for i_seq, line_split_element in enumerate(line_split):
      database_dict.setdefault(
        keys_dict[i_seq], fs()).append(line_split_element.lower())
  return database_dict

def select_database_dict_by_keys(select_keys, database_dict):
  result = {}
  for key in select_keys:
    result.setdefault(key, database_dict[key])
  return result

def _convert_to_numeric(values):
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

def _select_dict(database_dict, selection):
  result = {}
  for key in database_dict.keys():
    result.setdefault(key, database_dict[key].select(selection))
  return result

def order_by_value(database_dict, key, reverse = False):
  values = database_dict[key]
  new_values = _convert_to_numeric(values = values)
  selection = flex.sort_permutation(new_values, reverse = reverse)
  return _select_dict(database_dict = database_dict, selection = selection)

def leave_all_available_entries(database_dict, keys):
  selection = flex.bool(database_dict[keys[0]].size(), True)
  for key in keys:
    values = database_dict[key]
    selection &= values != "none"
  return _select_dict(database_dict = database_dict, selection = selection)

def filter_database(database_dict, key, value_min = None, value_max = None,
                    target_value = None):
  values = database_dict[key]
  if(target_value is not None):
    assert [value_min,value_max].count(None) == 2
    selection = values == target_value
  else:
    assert [value_min,value_max].count(None) < 2
    if([value_min,value_max].count(None) == 0): assert value_min < value_max
    new_values = _convert_to_numeric(values = values)
    selection = flex.bool(new_values.size(), True)
    if(value_min is not None):
      selection &= new_values > value_min
    if(value_max is not None):
      selection &= new_values < value_max
  return _select_dict(database_dict = database_dict, selection = selection)

def filter_histogram_of_key_value(database_dict, key, max_reject_fraction,
                                  edge_tolerance_small = 1.e-4):
  values = database_dict[key]
  new_values = _convert_to_numeric(values = values)
  #print flex.min(new_values), flex.max(new_values)
  size = new_values.size()
  #print size
  if(size == 0): return
  while True:
    values = database_dict[key]
    new_values = _convert_to_numeric(values = values)
    selection = flex.bool(new_values.size(), True)
    histogram = flex.histogram(data = new_values, n_slots = 3)
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
    database_dict = _select_dict(database_dict = database_dict, selection = selection)
  return database_dict

def show_histogram(data, n_slots):
  histogram = flex.histogram(data = data, n_slots = n_slots)
  l = histogram.data_min()
  for i, s in enumerate(histogram.slots()):
    r = histogram.data_min() + histogram.slot_width() * (i+1)
    print "%8.4f %8.4f %d" % (l, r, s)
    l = r

def polygon(file_name,
            selected_keys,
            filter_keys_and_targets,
            max_reject_fraction,
            n_histogram_slots):
  database_dict = extract_database_from_file(file_name = file_name)
  result = leave_all_available_entries(
    database_dict = database_dict,
    keys          = selected_keys)
  for filter_keys_and_target in filter_keys_and_targets:
    result = filter_database(
      database_dict = result,
      key           = filter_keys_and_target.key,
      value_min     = filter_keys_and_target.value_min,
      value_max     = filter_keys_and_target.value_max,
      target_value  = filter_keys_and_target.target_value)
  result = select_database_dict_by_keys(
    select_keys   = selected_keys,
    database_dict = result)
  for selected_key in selected_keys:
    result = filter_histogram_of_key_value(
      database_dict       = result,
      key                 = selected_key,
      max_reject_fraction = max_reject_fraction)
  for selected_key in selected_keys:
    print selected_key
    show_histogram(data = _convert_to_numeric(values = result[selected_key]),
                   n_slots = n_histogram_slots)

