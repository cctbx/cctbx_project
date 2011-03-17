import boost.python
iotbx_shelx_ext = boost.python.import_ext("iotbx_shelx_ext")
import sys
from cctbx.array_family import flex

def miller_array_export_as_shelx_hklf(
      miller_array,
      file_object=None,
      normalise_if_format_overflow=False):
  """\
  If the maximum data value does not fit into the f8.2/f8.0 format:
  normalise_if_format_overflow=False: RuntimeError is thrown
  normalise_if_format_overflow=True: data is normalised to the largest
  number to fit f8.2/f8.0 format
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
  if (min_val < 0):
    s = "%8.0f" % min_val
    if (len(s) > 8):
      if (not normalise_if_format_overflow):
        raise_f8_overflow(min_val)
      min_sc = -9999999. / min_val
  if (max_val > 0):
    s = "%8.0f" % max_val
    if (len(s) > 8):
      if (not normalise_if_format_overflow):
        raise_f8_overflow(max_val)
      max_sc = 99999999. / max_val
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
          result = "%8.0f" % v
          assert len(result) == 8
      return result
    line = fmt_3i4(h) + fmt_f8(data[i]*scale) + fmt_f8(s*scale)
    print >> file_object, line
  print >> file_object, "   0   0   0    0.00    0.00"

class reader(iotbx_shelx_ext.hklf_reader):

  def __init__(self, file_object=None, file_name=None):
    assert [file_object, file_name].count(None) == 1
    if (file_object is None):
      from libtbx import smart_open
      file_object = smart_open.for_reading(file_name=file_name)
    from cctbx.array_family import flex
    super(reader, self).__init__(lines=flex.split_lines(file_object.read()))

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    from cctbx import miller
    from cctbx import crystal
    if (crystal_symmetry is None):
      crystal_symmetry = crystal.symmetry()
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="shelx_hklf")
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=self.indices()).auto_anomalous()
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
