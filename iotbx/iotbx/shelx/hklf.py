from cctbx.array_family import flex
import sys

def write(miller_array, file_object=sys.stdout):
  data = miller_array.data()
  sigmas = miller_array.sigmas()
  s = 0.01
  fmt = "%4d%4d%4d%8.2f%8.2f"
  for i,h in miller_array.indices().items():
    if (sigmas != None): s = sigmas[i]
    print >> file_object, fmt % (h + (data[i],s))
  print >> file_object, fmt % (0,0,0,0,0)

class read:

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
