from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
iotbx_shelx_ext = bp.import_ext("iotbx_shelx_ext")
import sys
from cctbx.array_family import flex

def miller_array_export_as_shelx_hklf(
      miller_array,
      file_object=None,
      normalise_if_format_overflow=False,
      full_dynamic_range=False,
      scale_range=None):
  """\
  For the full_dynamic_range option, normalise data outside the range that
  fits into 8.6g format, regardless of settings for normalise_if_format_overflow
  or scale_range
  Otherwise, if the maximum data value does not fit into the f8.2/f8.0 format:
    normalise_if_format_overflow=False: RuntimeError is thrown
    normalise_if_format_overflow=True: data is normalised to the largest
      number to fit f8.2/f8.0 format or within specified scale_range
  """
  assert miller_array.is_real_array()
  if (file_object is None): file_object = sys.stdout
  def raise_f8_overflow(v):
    raise RuntimeError("SHELX HKL file F8.2/F8.0 format overflow: %.8g" % v)
  data = miller_array.data()
  sigmas = miller_array.sigmas()
  assert data is not None
  min_val = flex.min(data)
  max_val = flex.max(data)
  if (sigmas is not None):
    min_val = min(min_val, flex.min(sigmas))
    max_val = max(max_val, flex.max(sigmas))
  min_sc = 1
  max_sc = 1
  scale = 1
  if full_dynamic_range:
    max_abs = 999999.
    max_val = max(abs(max_val),abs(min_val))
    if (max_val > max_abs):
      scale = max_abs / max_val
  else:
    if scale_range is None:
      scale_range = (-999999., 9999999.)
    if (min_val < scale_range[0]):
      if not normalise_if_format_overflow:
        raise_f8_overflow(min_val)
      min_sc = scale_range[0] / min_val
    if (max_val > scale_range[1]):
      if (not normalise_if_format_overflow):
        raise_f8_overflow(max_val)
      max_sc = scale_range[1] / max_val
    scale = min(min_sc, max_sc)
  sigmas = miller_array.sigmas()
  s = 0.01
  for i,h in enumerate(miller_array.indices()):
    if (sigmas is not None): s = sigmas[i]
    def fmt_3i4(h):
      result = "%4d%4d%4d" % h
      if (len(result) != 12):
        raise RuntimeError(
          "SHELXL HKL file 3I4 format overflow: %s" % result)
      return result
    def fmt_f8(v):
      result = "%8.2f" % v
      if (len(result) != 8):
        result = "%8.1f" % v
        if (len(result) != 8):
          result = "%7d." % round(v)
          assert len(result) == 8
      return result
    def fmt_fullrange_data(v):
      if (v < 0.):
        if (abs(v) < 1.):
          result = "%8.5f" % v
        else:
          result = "%8.6g" % v
          if ("." not in result): # FORTRAN F8.2 interprets -123456 as -1234.56
            result = "%7d." % round(v)
      else:
        if (v < 1.):
          result = "%8.6f" % v
        else:
          result = "%8.7g" % v
      if ("." not in result): # FORTRAN F8.2 interprets 1234567 as 12345.67
        result = "%7d." % round(v) # Add "."
      assert len(result) == 8
      return result
    def fmt_fullrange_sigma(v):
      if (abs(v) >= 1.):
        result = "%8.7g" % v
        if ("." not in result):
          result = "%7d." % round(v)
      elif (abs(v) < 0.00001):
        result = "%8.5g" % 0.00001
      else:
        result = "%8.6g" % v
      return result
    if (full_dynamic_range):
      line = fmt_3i4(h) + fmt_fullrange_data(data[i]*scale) + fmt_fullrange_sigma(s*scale)
    else:
      line = fmt_3i4(h) + fmt_f8(data[i]*scale) + fmt_f8(s*scale)
    print(line, file=file_object)
  print("   0   0   0    0.00    0.00", file=file_object)

class reader(iotbx_shelx_ext.hklf_reader):

  def __init__(self, file_object=None, file_name=None):
    assert [file_object, file_name].count(None) == 1
    if (file_object is None):
      from libtbx import smart_open
      with smart_open.for_reading(file_name=file_name) as f:
        lines = f.read()
    else:
      lines = file_object.read()
      file_object.close()
    from cctbx.array_family import flex
    super(reader, self).__init__(lines=flex.split_lines(lines))

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None,
        anomalous=None):
    from cctbx import miller
    from cctbx import crystal
    if (crystal_symmetry is None):
      crystal_symmetry = crystal.symmetry()
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="shelx_hklf")
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=self.indices(), anomalous_flag=anomalous)
    if anomalous is None:
      miller_set = miller_set.auto_anomalous()
    miller_arrays = []
    obs = (miller.array(
      miller_set=miller_set,
      data=self.data(),
      sigmas=self.sigmas())
      .set_info(base_array_info.customized_copy(labels=["obs", "sigmas"])))
    miller_arrays.append(obs)
    if (self.alphas() is not None):
      miller_arrays.append(miller.array(
        miller_set=miller_set,
        data=self.alphas())
        .set_info(base_array_info.customized_copy(labels=["alphas"])))
    return miller_arrays

  def _override(method_name):
    def f(self):
      result = getattr(super(reader, self), method_name)()
      if result.size() == 0:
        return None
      return result
    return f
  alphas = _override('alphas')
  batch_numbers = _override('batch_numbers')
  wavelengths = _override('wavelengths')
  del _override
