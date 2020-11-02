from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex

class index_fobs_sigma_line(object):

  def __init__(self, raw_line):
    self.is_complete = False
    flds = raw_line.replace("="," ").split()
    if (len(flds) != 8): return
    if (flds[0].lower() not in ("inde", "index")): return
    if (flds[4].lower() != "fobs"): return
    if (flds[6].lower() != "sigma"): return
    self.names = [flds[4], flds[6]]
    try: self.index = tuple([int(i) for i in flds[1:4]])
    except Exception: return
    try: self.fobs = float(flds[5])
    except Exception: return
    try: self.sigma = float(flds[7])
    except Exception: return
    self.is_complete = True

class reader(object):

  def __init__(self, file_name=None, file_object=None, max_header_lines=30):
    assert [file_name, file_object].count(None) == 1
    if (file_object is None):
      file_object = open(file_name)
    lines = file_object.readlines()
    file_object.close()
    self._names = None
    self._indices = flex.miller_index()
    self._data = flex.double()
    self._sigmas = flex.double()
    have_data = False
    self.n_lines = 0
    for raw_line in lines:
      self.n_lines += 1
      ifs = index_fobs_sigma_line(raw_line)
      if (not ifs.is_complete):
        if (raw_line.strip().lower() == "end"):
          break
        if (self.n_lines == max_header_lines or have_data):
          raise RuntimeError("Unkown file format.")
      else:
        if (self._names is None): self._names = ifs.names
        self._indices.append(ifs.index)
        self._data.append(ifs.fobs)
        self._sigmas.append(ifs.sigma)
        have_data = True
    if (not have_data):
      raise RuntimeError("No data found in file.")

  def indices(self):
    return self._indices

  def data(self):
    return self._data

  def sigmas(self):
    return self._sigmas

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None,
        anomalous=None):
    if (crystal_symmetry is None):
      crystal_symmetry = crystal.symmetry()
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="cns_index_fobs_sigma")
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=self.indices(), anomalous_flag=anomalous)
    if anomalous is None:
      miller_set = miller_set.auto_anomalous()
    return [miller.array(
      miller_set=miller_set,
      data=self.data(),
      sigmas=self.sigmas())
      .set_info(base_array_info.customized_copy(labels=self._names))
      .set_observation_type_xray_amplitude()]
