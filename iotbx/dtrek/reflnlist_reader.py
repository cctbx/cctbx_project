from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
import sys

import boost_adaptbx.boost.python as bp
from six.moves import range
dtrek_ext = bp.import_ext("iotbx_dtrek_ext")

class reflnlist(object):

  def __init__(self, file, header_only=False):
    self.file = file
    self.line_no = 0
    line_1 = self.next_line().split()
    assert len(line_1) in (3,4)
    line_1 = [int(field) for field in line_1]
    if (len(line_1) == 3):
      line_1.append(0)
    self.NumInts = line_1[0]
    assert self.NumInts >= 3
    self.NumFloats = line_1[1]
    self.NumStrings = line_1[2]
    self.NumInfoLines = line_1[3]
    self.InfoLines = []
    for i in range(self.NumInfoLines):
      self.InfoLines.append(self.next_line().strip())
    self.column_info = []
    num_columns = self.NumInts + self.NumFloats + self.NumStrings
    self.column_names = []
    for i in range(num_columns):
      column_name = self.next_line().strip()
      assert len(column_name) > 1
      assert column_name[0] in "nfs"
      if (i < self.NumInts):
        assert column_name[0] == "n"
      elif (i < self.NumInts + self.NumFloats):
        assert column_name[0] == "f"
      else:
        assert column_name[0] == "s"
      self.column_names.append(column_name)
    assert "fIntensity" in self.column_names
    assert "fSigmaI" in self.column_names
    assert self.column_names[:3] == ["nH", "nK", "nL"]
    if (header_only):
      self.miller_indices = None
      self.column_dict = None
    else:
      self.miller_indices = flex.miller_index()
      self.column_dict = {}
      column_list = [None, None, None]
      for i in range(3, num_columns):
        if (i < self.NumInts):
          column_list.append(flex.int())
        elif (i < self.NumInts + self.NumFloats):
          column_list.append(flex.double())
        else:
          column_list.append(flex.std_string())
        self.column_dict[self.column_names[i]] = column_list[-1]
      while 1:
        line = self.next_line()
        if (line is None):
          break
        data_values = line.split()
        assert len(data_values) == num_columns, "Line no. %d" % self.line_no
        h = [int(v) for v in data_values[:3]]
        self.miller_indices.append(h)
        offset = 3
        for i in range(offset, self.NumInts):
          v = int(data_values[i])
          column_list[i].append(v)
        offset = self.NumInts
        for i in range(offset, offset+self.NumFloats):
          v = float(data_values[i])
          column_list[i].append(v)
        offset += self.NumFloats
        for i in range(offset, num_columns):
          column_list[i].append(data_values[i])
    del self.file
    del self.line_no

  def next_line(self):
    line = self.file.readline()
    if (line == ""): return None
    self.line_no += 1
    return line

  def space_group_info(self):
    for info_line in self.InfoLines:
      if (info_line.upper().startswith("CRYSTAL_SPACEGROUP=")):
        assert info_line[-1] == ";"
        space_group_number = info_line[:-1].split("=", 1)[1]
        return sgtbx.space_group_info(symbol=space_group_number)
    return None

  def unit_cell(self):
    for info_line in self.InfoLines:
      if (info_line.upper().startswith("CRYSTAL_UNIT_CELL=")):
        assert info_line[-1] == ";"
        unit_cell = info_line[:-1].split("=", 1)[1]
        return uctbx.unit_cell(unit_cell)
    return None

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=self.unit_cell(),
      space_group_info=self.space_group_info())

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    self.crystal_symmetry().show_summary(f=out, prefix=prefix)
    print(prefix + "Column names:", " ".join(self.column_names[3:]), file=out)
    print(prefix + "Number of reflections:", self.miller_indices.size(), file=out)

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None,
        anomalous=None):
    crystal_symmetry_from_file = self.crystal_symmetry()
    crystal_symmetry = crystal_symmetry_from_file.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=force_symmetry)
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="dtrek_reflnlist")
    miller_arrays = []
    sigmas=self.column_dict["fSigmaI"]
    miller_arrays.append(miller.array(
      miller_set=miller.set(
        crystal_symmetry=crystal_symmetry,
        indices=self.miller_indices,
        anomalous_flag=False),
      data=self.column_dict["fIntensity"],
      sigmas=sigmas)
      .select(sigmas > 0)
      .set_info(base_array_info.customized_copy(
        labels=["Intensity", "SigmaI"],
        crystal_symmetry_from_file=crystal_symmetry_from_file))
      .set_observation_type_xray_intensity())
    if ("fIntensity+" in self.column_dict):
      assert "fSigmaI+" in self.column_dict
      assert "fIntensity-" in self.column_dict
      assert "fSigmaI-" in self.column_dict
      if (crystal_symmetry.space_group_info() is None):
        space_group = sgtbx.space_group() # P 1
      else:
        space_group = crystal_symmetry.space_group()
      ac = dtrek_ext.anomalous_combined(
        space_group,
        self.miller_indices,
        self.column_dict["fIntensity+"],
        self.column_dict["fSigmaI+"],
        self.column_dict["fIntensity-"],
        self.column_dict["fSigmaI-"])
      miller_arrays.append(miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=ac.miller_indices(),
          anomalous_flag=True),
        data=ac.data(),
        sigmas=ac.sigmas())
        .set_info(base_array_info.customized_copy(
          labels=["Intensity+-", "SigmaI+-"],
          crystal_symmetry_from_file=crystal_symmetry_from_file))
        .set_observation_type_xray_intensity())
    for column_name in self.column_names:
      if (column_name in ("nH", "nK", "nL",
                          "fSigmaI", "fIntensity",
                          "fIntensity+", "fSigmaI+",
                          "fIntensity-", "fSigmaI-")):
        continue
      miller_arrays.append(miller.array(
        miller_set=miller.set(
          crystal_symmetry=crystal_symmetry,
          indices=self.miller_indices,
          anomalous_flag=False),
        data=self.column_dict[column_name])
        .set_info(base_array_info.customized_copy(
          labels=[column_name[1:]],
          crystal_symmetry_from_file=crystal_symmetry_from_file)))
    return miller_arrays
