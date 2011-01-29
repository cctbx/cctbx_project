from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
import iotbx_shelx_ext
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

class reader_base(object):

  def __init__(self, file_object=None, filename=None, strict=True):
    assert [file_object, filename].count(None) == 1
    if file_object is None:
      file_object = open(filename)
    super(reader_base, self).__init__(content=file_object.read(),
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
      result = getattr(super(reader_base, self), method_name)()
      if result.size() == 0:
        return None
      return result
    return f
  alphas = _override('alphas')
  batch_numbers = _override('batch_numbers')
  del _override

class fast_reader(reader_base, iotbx_shelx_ext.hklf_reader): pass

class simple_reader(reader_base, iotbx_shelx_ext.simple_hklf_reader): pass


class python_reader(reader_base):

  def __init__(self, file_object=None, filename=None):
    assert [file_object, filename].count(None) == 1
    if file_object is None:
      file_object = open(filename)
    self._indices = flex.miller_index()
    self._data = flex.double()
    self._sigmas = flex.double()
    self._alphas = None
    for line in file_object:
      if (len(line.strip()) == 0): break
      if (len(line.rstrip()) > 32):
        raise RuntimeError("Not a SHELX hklf file.")
      h = [int(line[i*4:(i+1)*4]) for i in xrange(3)]
      if (h == [0,0,0]): break
      fs = [float(line[12+i*8:12+(i+1)*8]) for i in xrange(2)]
      self._indices.append(h)
      self._data.append(fs[0])
      self._sigmas.append(fs[1])
      if (self._indices.size() == 1):
        try: a = int(line[28:32])
        except: pass
        else: self._alphas = flex.int(1, a)
      elif (self._alphas is not None):
        self._alphas.append(int(line[28:32]))
    if (self._indices.size() == 0):
      raise RuntimeError("No data in SHELX hklf file.")

  def indices(self):
    return self._indices

  def data(self):
    return self._data

  def sigmas(self):
    return self._sigmas

  def alphas(self):
    return self._alphas

  def batch_numbers(self):
    return self._alphas

reader = fast_reader
