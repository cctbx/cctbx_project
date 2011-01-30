import boost.python
iotbx_shelx_ext = boost.python.import_ext("iotbx_shelx_ext")
import sys

def miller_array_export_as_shelx_hklf(miller_array, file_object=None):
  assert miller_array.is_real_array()
  if (file_object is None): file_object = sys.stdout
  data = miller_array.data()
  sigmas = miller_array.sigmas()
  s = 0.01
  fmt = "%4d%4d%4d%8.2f%8.2f"
  for i,h in enumerate(miller_array.indices()):
    if (sigmas is not None): s = sigmas[i]
    print >> file_object, fmt % (h + (data[i],s))
  print >> file_object, fmt % (0,0,0,0,0)

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
