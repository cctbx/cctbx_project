"Reading and writing of scalepack merge reflection files."

# Sample scalepack OUTPUT FILE
#    1
# -987
#    34.698    71.491    54.740    90.000   106.549    90.000 p21
#   0   0   4  3617.6   287.2
#   0   1   6 12951.7  1583.6 12039.2  1665.8
#
# Format: (3I4, 4F8.1)

from cctbx import uctbx
from cctbx import sgtbx
from cctbx import crystal
from cctbx import miller
from cctbx.array_family import flex
from libtbx.math_utils import iround
from libtbx import easy_pickle
import exceptions
import os
import sys

class FormatError(exceptions.Exception): pass

class reader(object):

  def __init__(self, file_handle, header_only=False):
    line = file_handle.readline()
    if (line.rstrip() != "    1"):
      raise FormatError, "line 1: expecting '    1'"
    file_handle.readline() # ignore line 2
    line_error = "line 3: expecting unit cell parameters and space group label"
    line = file_handle.readline()
    if (len(line) < 63 or line[60] != ' '):
      raise FormatError, line_error
    try:
      uc_params = [float(line[i * 10 : (i + 1) * 10]) for i in xrange(6)]
    except KeyboardInterrupt: raise
    except:
      raise FormatError, line_error
    self.unit_cell = uctbx.unit_cell(uc_params)
    self.space_group_symbol = line[61:].strip()
    if (len(self.space_group_symbol) == 0):
      raise FormatError, line_error
    try:
      self.space_group_info = sgtbx.space_group_info(self.space_group_symbol)
    except KeyboardInterrupt: raise
    except:
      self.space_group_info = None
    if (header_only): return
    self.miller_indices = flex.miller_index()
    self.i_obs = flex.double()
    self.sigmas = flex.double()
    self.anomalous = False
    line_count = 3
    while 1:
      line = file_handle.readline()
      line_count += 1
      line_error = "line %d: expecting (3I4, 4F8.1)" % line_count
      if (line == ""): break
      line = line.rstrip()
      if (len(line) == 0): continue
      line += (" " * 44)
      flds = []
      used = 0
      for width in (4,4,4,8,8,8,8):
        next_used = used + width
        flds.append(line[used:next_used].strip())
        used = next_used
      try:
        h = [int(flds[i]) for i in xrange(3)]
      except KeyboardInterrupt: raise
      except:
        raise FormatError, line_error
      for i in (0,1):
        j = 3+2*i
        if (len(flds[j])):
          try:
            i_obs, sigma = (float(flds[j]), float(flds[j+1]))
          except KeyboardInterrupt: raise
          except:
            raise FormatError, line_error
          if (i):
            h = [-e for e in h]
            self.anomalous = True
          if (sigma != 0 or (i_obs != 0 and i_obs != 1)):
            self.miller_indices.append(h)
            self.i_obs.append(i_obs)
            self.sigmas.append(sigma)

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    self.unit_cell.show_parameters(f=out, prefix=prefix+"Unit cell: ")
    print >> out, prefix + "Space group symbol:", self.space_group_symbol
    print >> out, prefix + "Anomalous flag:", self.anomalous
    print >> out, prefix + "Number of reflections:", self.miller_indices.size()

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group_info=self.space_group_info,
      correct_rhombohedral_setting_if_necessary=True)

  def as_miller_array(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="scalepack_merge")
    crystal_symmetry_from_file = self.crystal_symmetry()
    return (miller.array(
      miller_set=miller.set(
        crystal_symmetry=crystal_symmetry_from_file.join_symmetry(
          other_symmetry=crystal_symmetry,
          force=force_symmetry),
        indices=self.miller_indices,
        anomalous_flag=self.anomalous),
      data=self.i_obs,
      sigmas=self.sigmas)
      .set_info(base_array_info.customized_copy(
        labels=["i_obs", "sigma"],
        crystal_symmetry_from_file=crystal_symmetry_from_file))
      .set_observation_type_xray_intensity())

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    return [self.as_miller_array(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      merge_equivalents=merge_equivalents,
      base_array_info=base_array_info)]

def format_f8_1_or_i8(h, label, value):
  if (value < 1.e6): return "%8.1f" % value
  result = "%8d" % iround(value)
  if (len(result) > 8):
    raise ValueError(
      "Value is too large for scalepack merge format: hkl=%s, %s=%.6g" % (
        str(h).replace(" ",""), label, value))
  return result

def write(file_name=None, file_object=None, miller_array=None,
          space_group_symbol=None,
          line_1="    1",
          line_2=" -987"):
  assert [file_name, file_object].count(None) == 1
  assert miller_array.is_xray_intensity_array() or miller_array.is_xray_amplitude_array()
  assert miller_array.sigmas() is not None
  if (file_object is None):
    file_object = open(file_name, "w")
  if (not miller_array.is_xray_intensity_array()):
    miller_array = miller_array.f_as_f_sq()
  if (space_group_symbol is None):
    space_group_symbol = str(miller_array.space_group_info())
    if (not space_group_symbol.endswith(")")): # not universal Hermann-Mauguin
      space_group_symbol = space_group_symbol.replace(" ", "").lower()
  print >> file_object, line_1
  print >> file_object, line_2
  print >> file_object, ("%10.3f"*6) % miller_array.unit_cell().parameters(),
  print >> file_object, space_group_symbol
  if (not miller_array.anomalous_flag()):
    for h,f,s in zip(miller_array.indices(),
                     miller_array.data(),
                     miller_array.sigmas()):
      print >> file_object, ((("%4d"*3) % h)
        + format_f8_1_or_i8(h, "intensity", f)
        + format_f8_1_or_i8(h, "sigma", s))
  else:
    asu, matches = miller_array.match_bijvoet_mates()
    sel_pairs_plus = matches.pairs_hemisphere_selection("+")
    sel_pairs_minus = matches.pairs_hemisphere_selection("-")
    indices = asu.indices().select(sel_pairs_plus)
    data_plus = asu.data().select(sel_pairs_plus)
    data_minus = asu.data().select(sel_pairs_minus)
    sigmas_plus = asu.sigmas().select(sel_pairs_plus)
    sigmas_minus = asu.sigmas().select(sel_pairs_minus)
    for h,fp,sp,fm,sm in zip(indices, data_plus, sigmas_plus,
                                      data_minus, sigmas_minus):
      print >> file_object, ((("%4d"*3) % h)
        + format_f8_1_or_i8(h, "intensity", fp)
        + format_f8_1_or_i8(h, "sigma", sp)
        + format_f8_1_or_i8(h, "intensity", fm)
        + format_f8_1_or_i8(h, "sigma", sm))
    sel_singles = matches.singles_hemisphere_selection("+")
    indices = asu.indices().select(sel_singles)
    data = asu.data().select(sel_singles)
    sigmas = asu.sigmas().select(sel_singles)
    for h,f,s in zip(indices, data, sigmas):
      print >> file_object, ((("%4d"*3) % h)
        + format_f8_1_or_i8(h, "intensity", f)
        + format_f8_1_or_i8(h, "sigma", s))
    sel_singles = matches.singles_hemisphere_selection("-")
    indices = -asu.indices().select(sel_singles)
    data = asu.data().select(sel_singles)
    sigmas = asu.sigmas().select(sel_singles)
    for h,f,s in zip(indices, data, sigmas):
      print >> file_object, ((("%4d"*3) % h) + (" "*16)
        + format_f8_1_or_i8(h, "intensity", f)
        + format_f8_1_or_i8(h, "sigma", s))

def run(args):
  to_pickle = "--pickle" in args
  for file_name in args:
    if (file_name.startswith("--")): continue
    s = reader(open(file_name, "r"))
    miller_array = s.as_miller_array(info="From file: "+file_name)
    miller_array.show_summary()
    if (to_pickle):
      pickle_file_name = os.path.split(file_name)[1] + ".pickle"
      print "Writing:", pickle_file_name
      easy_pickle.dump(pickle_file_name, miller_array)
    print
