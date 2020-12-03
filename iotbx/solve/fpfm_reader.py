from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
from six.moves import range

class reader(object):

  def __init__(self, file_name=None, file_object=None, max_header_lines=30):
    assert [file_name, file_object].count(None) == 1
    if (file_object is None):
      file_object = open(file_name)
    self.file_object = file_object
    line = self.file_object.readline()
    n_data_columns = int(line)
    assert n_data_columns == 4, "Exactly 4 data columns expected."
    self.file_object.readline() # ignore the FORTRAN format line
    label = None
    for sign in ("+", "-"):
      flds = self.file_object.readline().split()
      assert flds[0] == "F"+sign, "F"+sign+" header line expected."
      if (label is None):
        label = flds[-1]
      else:
        assert flds[-1] == label, "Inconsistent wavelength label."
      flds = self.file_object.readline().split()
      assert flds[:3] == ("sig of F"+sign).split(), \
        "sig of F"+sign+" header line expected."
      assert flds[-1] == label, "Inconsistent wavelength label."
    self._indices = flex.miller_index()
    self._data = flex.double()
    self._sigmas = flex.double()
    self.n_lines = 6
    for raw_line in self.file_object:
      self.n_lines += 1
      try:
        h = [int(raw_line[1+i*4:1+i*4+4]) for i in range(3)]
        d = [float(raw_line[23+i*10:23+i*10+10]) for i in range(4)]
      except Exception as e:
        raise RuntimeError("Line %d: %s" % (self.n_lines, str(e)))
      for i in range(2):
        assert (d[2*i] < 0) == (d[2*i+1] < 0), \
          "Line %d: Inconsistent observation and sigma." % self.n_lines
      if (d[0] >= 0):
        self._indices.append(h)
        self._data.append(d[0])
        self._sigmas.append(d[1])
      if (d[2] >= 0):
        self._indices.append([-e for e in h])
        self._data.append(d[2])
        self._sigmas.append(d[3])

  def __del__(self):
    self.file_object.close()

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
      base_array_info = miller.array_info(source_type="solve_fpfm")
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=self.indices(),
      anomalous_flag=True)
    return [miller.array(
      miller_set=miller_set,
      data=self.data(),
      sigmas=self.sigmas())
      .set_info(base_array_info.customized_copy(
        labels=["fpfm", "sigma_fpfm"]))
      .set_observation_type_xray_amplitude()]
