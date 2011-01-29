from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
import boost.python
iotbx_shelx_ext = boost.python.import_ext("iotbx_shelx_ext")
import sys

def miller_export_as_shelx_hklf(self, file_object=None):
  assert self.is_real_array()
  if (file_object is None): file_object = sys.stdout
  data = self.data()
  sigmas = self.sigmas()
  s = 0.01
  fmt = "%4d%4d%4d%8.2f%8.2f"
  for i,h in enumerate(self.indices()):
    if (sigmas is not None): s = sigmas[i]
    print >> file_object, fmt % (h + (data[i],s))
  print >> file_object, fmt % (0,0,0,0,0)

miller.array.export_as_shelx_hklf = miller_export_as_shelx_hklf

class reader(iotbx_shelx_ext.hklf_reader):

  def __init__(self, file_object=None, filename=None, strict=True):
    assert [file_object, filename].count(None) == 1
    if (file_object is None):
      file_object = open(filename)
    super(reader, self).__init__(
      lines=flex.split_lines(file_object.read()),
      strict=strict)

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
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
  del _override
