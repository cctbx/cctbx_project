from cctbx import miller
from cctbx.array_family import flex
import sys

def miller_export_as_shelx_hklf(self, file_object=sys.stdout):
  data = self.data()
  sigmas = self.sigmas()
  s = 0.01
  fmt = "%4d%4d%4d%8.2f%8.2f"
  for i,h in self.indices().items():
    if (sigmas != None): s = sigmas[i]
    print >> file_object, fmt % (h + (data[i],s))
  print >> file_object, fmt % (0,0,0,0,0)

miller.miller_export_as_shelx_hklf = miller_export_as_shelx_hklf

class reader:

  def __init__(self, file_object):
    self._indices = flex.miller_index()
    self._data = flex.double()
    self._sigmas = flex.double()
    self._alphas = flex.int()
    self._count_alphas = 0
    for line in file_object:
      if (len(line.strip()) == 0): break
      h = [int(line[i*4:(i+1)*4]) for i in xrange(3)]
      fs = [float(line[12+i*8:12+(i+1)*8]) for i in xrange(2)]
      try: a = int(line[28:32])
      except: a = 0
      else: self._count_alphas += 1
      if (h == [0,0,0]): break
      self._indices.append(h)
      self._data.append(fs[0])
      self._sigmas.append(fs[1])
      self._alphas.append(a)

  def indices(self):
    return self._indices

  def data(self):
    return self._data

  def sigmas(self):
    return self._sigmas

  def alphas(self):
    return self._alphas

  def count_alphas(self):
    return self._count_alphas

  def as_miller_arrays(self, crystal_symmetry=None, force_symmetry=00000,
                             info_prefix=""):
    if (crystal_symmetry == None):
      crystal_symmetry = crystal.symmetry()
    miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=self.indices()).auto_anomalous()
    miller_arrays = []
    obs = miller.array(
      miller_set=miller_set,
      data=self.data(),
      sigmas=self.sigmas(),
      info=info_prefix+"obs,sigmas")
    miller_arrays.append(obs)
    if (self.count_alphas() > 0):
      miller_arrays.append(miller.array(
        miller_set=miller_set,
        data=self.alphas(),
        info=info_prefix+"alphas"))
    return miller_arrays
