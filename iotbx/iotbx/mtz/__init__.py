from cctbx.uctbx import unit_cell

from scitbx.python_utils.misc import import_regular_symbols
from iotbx_boost import mtz as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols

from iotbx.mtz import writer
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from scitbx.python_utils.misc import adopt_init_args
from scitbx.python_utils.str_utils import overwrite_at

column_type_legend_source = \
  "http://www.ccp4.ac.uk/dist/html/mtzlib.html#fileformat"
column_type_legend = {
  "H": "index h,k,l",
  "J": "intensity",
  "F": "amplitude",
  "D": "anomalous difference",
  "Q": "standard deviation",
  "G": "F(+) or F(-)",
  "L": "standard deviation",
  "K": "I(+) or I(-)",
  "M": "standard deviation",
  "P": "phase angle in degrees",
  "W": "weight (of some sort)",
  "A": "phase probability coefficients (Hendrickson/Lattmann)",
  "B": "BATCH number",
  "Y": "M/ISYM, packed partial/reject flag and symmetry number",
  "I": "integer",
  "R": "real",
}

class Mtz (ext.Mtz):
  def __init__(self,s):
    ext.Mtz.__init__(self,s)

  def __getattr__(self,s):
    assert type(s) == type(str())
    return self.getColumn(s)

  def label_to_crystal(self, label):
    assert label in self.columns()
    for i in xrange(self.ncrystals()):
      cryst = self.getCrystal(i)
      for j in xrange(cryst.ndatasets()):
        data = cryst.getDataset(j)
        for k in xrange(data.ncolumns()):
          if data.getColumn(k).label() == label:
            return cryst

  def get_space_group_info(self):
    if (self.nsym() == 0):
      return None
    return sgtbx.space_group_info(group=self.getSgtbxSpaceGroup())

  def as_miller_arrays(self, crystal_symmetry=None, force_symmetry=00000,
                             info_prefix=""):
    other_symmetry = crystal_symmetry
    result = []
    for i_crystal in xrange(self.ncrystals()):
      cryst = self.getCrystal(i_crystal)
      crystal_symmetry = crystal.symmetry(
        unit_cell=cryst.get_unit_cell(),
        space_group_info=self.get_space_group_info()).join_symmetry(
          other_symmetry=other_symmetry,
          force=force_symmetry)
      for i_dataset in xrange(cryst.ndatasets()):
        dataset = cryst.getDataset(i_dataset)
        column_groups = self.group_columns(crystal_symmetry, dataset)
        for column_group in column_groups:
          info = info_prefix + column_group.info()
          result.append(column_group.set_info(info))
    return result

  def group_columns(self, crystal_symmetry, dataset):
    known_mtz_column_types = "".join(column_type_legend.keys())
    assert len(known_mtz_column_types) == 16 # safety guard
    all_column_labels = dataset.all_column_labels()
    all_column_types = mend_non_conforming_anomalous_column_types(
      dataset.all_column_types(), all_column_labels)
    groups = []
    i_column = -1
    while 1:
      i_column += 1
      if (i_column == dataset.ncolumns()): break
      assert dataset.getColumn(i_column).type() in known_mtz_column_types
      l0 = all_column_labels[i_column]
      t0 = all_column_types[i_column]
      if (t0 == "H"): continue # skip h,k,l
      if (t0 in "BYI"): # integer columns
        groups.append(column_group(
          crystal_symmetry=crystal_symmetry,
          primary_column_type=t0,
          labels=[l0],
          indices=self.valid_indices(l0),
          anomalous_flag=False,
          data=self.valid_integers(l0)))
      elif (t0 in "R"): # general real column
        groups.append(column_group(
          crystal_symmetry=crystal_symmetry,
          primary_column_type=t0,
          labels=[l0],
          indices=self.valid_indices(l0),
          anomalous_flag=False,
          data=self.valid_values(l0)))
      elif (t0 in "A"): # Hendrickson-Lattman coefficients
        assert all_column_types[i_column:i_column+4] == "AAAA"
        labels = all_column_labels[i_column:i_column+4]
        i_column += 3
        groups.append(column_group(
          crystal_symmetry=crystal_symmetry,
          primary_column_type=t0,
          labels=labels,
          indices=self.valid_indices(l0),
          anomalous_flag=False,
          data=self.valid_hl(labels[0], labels[1], labels[2], labels[3])))
      elif (t0 in "JFD"):
        # "J": "intensity"
        # "F": "amplitude"
        # "D": "anomalous difference"
        # "Q": "standard deviation"
        # "P": "phase angle in degrees"
        labels = [l0]
        data = None
        sigmas = None
        if (    i_column+1 < len(all_column_types)
            and all_column_types[i_column+1] in "QP"):
          labels = all_column_labels[i_column:i_column+2]
          i_column += 1
          if (all_column_types[i_column] == "Q"):
            sigmas = self.valid_values(labels[-1])
          else:
            if (t0 == "J"):
              raise AssertionError, "Invalid column combination."
            data=self.valid_complex(l0, labels[-1])
        if (data is None):
          data=self.valid_values(l0)
        groups.append(column_group(
          crystal_symmetry=crystal_symmetry,
          primary_column_type=t0,
          labels=labels,
          indices=self.valid_indices(l0),
          anomalous_flag=False,
          data=data,
          sigmas=sigmas))
      elif (t0 in "GK"):
        # "G": "F(+) or F(-)"
        # "L": "standard deviation"
        # "P": "phase angle in degrees"
        # "K": "I(+) or I(-)"
        # "M": "standard deviation"
        perm = None
        remaining_types = all_column_types[i_column:]
        if (remaining_types[:4] in ("GLGL", "GPGP", "KMKM")):
          perm = [0,1,2,3]
        elif (remaining_types[:4] in ("GGLL", "GGPP", "KKMM")):
          perm = [0,2,1,3]
        elif (remaining_types[:2] in ("GG", "KK")):
          perm = [0,1]
        else:
          raise AssertionError, "Invalid column combination."
        labels = [all_column_labels[i_column+i] for i in perm]
        i_column += len(perm)-1
        if (len(perm) == 2):
          groups.append(column_group(
            crystal_symmetry=crystal_symmetry,
            primary_column_type=t0,
            labels=labels,
            indices=self.valid_indices(labels[0], labels[1]),
            anomalous_flag=True,
            data=self.valid_values(labels[0], labels[1])))
        elif ("P" in remaining_types[:4]):
          groups.append(column_group(
            crystal_symmetry=crystal_symmetry,
            primary_column_type=t0,
            labels=labels,
            indices=self.valid_indices(labels[0], labels[2]),
            anomalous_flag=True,
            data=self.valid_complex(labels[0],labels[1],labels[2],labels[3])))
        else:
          groups.append(column_group(
            crystal_symmetry=crystal_symmetry,
            primary_column_type=t0,
            labels=labels,
            indices=self.valid_indices(labels[0], labels[2]),
            anomalous_flag=True,
            data=self.valid_values(labels[0],labels[2]),
            sigmas=self.valid_values(labels[1],labels[3])))
      else:
        groups.append(column_group(
          crystal_symmetry=crystal_symmetry,
          primary_column_type=t0,
          labels=[l0],
          indices=self.valid_indices(l0),
          anomalous_flag=False,
          data=self.valid_values(l0)))
    return groups

def _Crystal_get_unit_cell(self):
  result = self.UnitCell()
  if (result.is_similar_to(uctbx.unit_cell((1,1,1,90,90,90)))):
    return None
  return result

Crystal.get_unit_cell = _Crystal_get_unit_cell

def _Dataset_all_column_types(self):
  result = ""
  for i_column in xrange(self.ncolumns()):
    result += self.getColumn(i_column).type()
  return result

def _Dataset_all_column_labels(self):
  result = []
  for i_column in xrange(self.ncolumns()):
    result.append(self.getColumn(i_column).label())
  return result

Dataset.all_column_types = _Dataset_all_column_types
Dataset.all_column_labels = _Dataset_all_column_labels

def column_group(crystal_symmetry, primary_column_type, labels,
                 indices, anomalous_flag,
                 data, sigmas=None):
  assert data is not None
  if (sigmas is not None): assert sigmas.size() == data.size()
  if (anomalous_flag is 0001):
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=indices,
      anomalous_flag=0001)
  else:
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=indices).auto_anomalous(min_fraction_bijvoet_pairs=2/3.)
  result = (miller.array(
    miller_set=miller_set,
    data=data,
    sigmas=sigmas)
    .set_info(",".join(labels)))
  if (primary_column_type in "FG"):
    result.set_observation_type_xray_amplitude()
  elif (primary_column_type in "JK"):
    result.set_observation_type_xray_intensity()
  return result

for v in writer.__dict__.values():
  if (hasattr(v, "func_name")):
    setattr(MtzWriter, v.func_name, v)

def mend_non_conforming_anomalous_column_types(all_types, all_labels):

  replacements = {
    "JQJQ": "KMKM",
    "FQFQ": "GLGL",
    "JJQQ": "KKMM",
    "FFQQ": "GGLL",
  }

  def contains_one_of(label, patterns):
    for pattern in patterns:
      if (label.find(pattern) >= 0):
        return 0001
    return 00000

  def are_anomalous_labels():
    if (group_types[1] == "Q"):
      permutation = ((0,1),(2,3))
    else:
      permutation = ((0,2),(1,3))
    for offsets,patterns in zip(permutation, (("+", "PLUS"), ("-", "MINU"))):
      for offs in offsets:
        label = all_labels[i_group_start + offs].upper()
        if (not contains_one_of(label, patterns)):
          return 00000
    return 0001

  def find_group():
    for group_types in replacements.keys():
      i_group_start = all_types_x.find(group_types)
      if (i_group_start >= 0):
        return group_types, i_group_start
    return None, None

  if (0): # for debugging only
    all_types = all_types.replace("GLGL", "FQFQ")
    all_types = all_types.replace("KMKM", "JQJQ")

  all_types_x = all_types
  while 1:
    group_types, i_group_start = find_group()
    if (group_types is None):
      break
    if (are_anomalous_labels()):
      replacement = replacements[group_types]
      all_types = overwrite_at(all_types, i_group_start, replacement)
    else:
      replacement = "X"
    all_types_x = overwrite_at(all_types_x, i_group_start, replacement)
  return all_types
