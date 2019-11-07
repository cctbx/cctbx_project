#
# periodogram.py
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
"""Calculate the periodogram of real evenly-spaced data"""

from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx import fftpack
from six.moves import range

class Kernel(object):
  """A discrete symmetric normalized smoothing kernel for use with kernapply.
  Currently allowed types are 'daniell' or 'modified.daniell', with 'm' giving
  the kernel dimension. Alternatively, provide the upper half of the smoothing
  kernel coefficients as 'coef' directly"""

  def __init__(self, name='daniell', m=2, coef=None):
    if coef is not None:
      self.m = len(coef) - 1
      self.name = None
      self.coef = flex.double(coef)
      if abs(2.0 * flex.sum(self.coef[1:]) + self.coef[0] - 1) > 1.e-10:
        raise ValueError('Coefficients do not add to 1')
      return

    if not name in ['daniell', 'modified.daniell']:
      raise TypeError('Unknown kernel type "{0}"'.format(name))

    self.m = m
    self.name = name
    if name == 'daniell': self._set_daniell()
    if name == 'modified.daniell': self._set_modified_daniell()
    return

  def _set_daniell(self):
    self.coef = flex.double(self.m + 1, 1./(2. * self.m + 1))

  def _set_modified_daniell(self):
    self.coef = flex.double(self.m + 1, 1./(2. * self.m))
    self.coef[self.m] = self.coef[self.m] * 0.5

def kernapply(x, k, circular=False):
  """Convolve a sequence x with a Kernel k"""

  x = flex.double(x).deep_copy()
  lenx = len(x)
  w = flex.double(lenx, 0.0)
  w.set_selected(flex.size_t_range(k.m + 1), k.coef)
  sel = lenx -1 - flex.size_t_range(k.m)
  w.set_selected(sel, k.coef[1:])

  # do convolution in the Fourier domain
  fft = fftpack.real_to_complex(lenx)
  n = fft.n_real()
  m = fft.m_real()
  x.extend(flex.double(m-n, 0.))
  w.extend(flex.double(m-n, 0.))
  conv = fft.forward(x) * fft.forward(w)

  # extend result by the reverse conjugate, omitting the DC offset and Nyquist
  # frequencies. If fft.n_real() is odd there is no Nyquist term.
  end = fft.n_complex() - (fft.n_real() + 1) % 2
  conv.extend(flex.conj(conv[1:end]).reversed())

  # transform back, take real part and scale
  fft = fftpack.complex_to_complex(len(conv))
  result = fft.backward(conv).parts()[0] / n

  if circular:
    return result
  else:
    return result[(k.m):(lenx-k.m)]

class Periodogram(object):
  """Calculate the periodogram of real evenly-spaced data. This class gives
  the same spectrum as R's spec.pgram function when called using
  spec.pgram(x, detrend=F, taper=0, fast=F)"""

  def __init__(self, x, spans=None, demean=True, detrend=True):

    # Ensure x is copied as it will be changed in-place
    x = flex.double(x).deep_copy()
    n = len(x)

    if detrend:
      t = flex.size_t_range(n).as_double() + 1 - (n + 1)/2
      inv_sumt2 = 1./t.dot(t)
      x = x - flex.mean(x) - x.dot(t) * t * inv_sumt2
    elif demean:
      x -= flex.mean(x)

    # determine frequencies
    stop = ((n - (n % 2)) // 2) + 1
    self.freq = flex.double([i / n for i in range(1, stop)])

    fft = fftpack.real_to_complex(n)
    n = fft.n_real()
    m = fft.m_real()
    x.extend(flex.double(m-n, 0.))
    xf = fft.forward(x)

    # get abs length of complex and normalise by n to get the raw periodogram
    spec = flex.norm(xf) / n

    if spans is None:
      # if not doing smoothing, the spectrum is just the raw periodogram with
      # the uninteresting DC offset removed
      self.spec = spec[1:]
      return

    # for smoothing replace the DC offset term and extend the rest of the
    # sequence by its reverse conjugate, omitting the Nyquist term if it is
    # present
    spec[0] = spec[1]
    end = fft.n_complex() - (fft.n_real() + 1) % 2
    spec.extend(spec[1:end].reversed())

    try:
      # multiple integer spans
      nspans = len(spans)
      m = int(spans[0]) // 2
      multiple = True
    except TypeError:
      # single integer span
      m = int(spans) // 2
      multiple = False

    # construct smoothing kernel
    k = Kernel('modified.daniell', m)

    if multiple:
      for i in range(1, nspans):
        # construct kernel for convolution
        m1 = int(spans[i]) // 2
        k1 = Kernel('modified.daniell', m1)

        # expand coefficients of k to full kernel and zero pad for smoothing
        x1 = flex.double(k1.m, 0.0)
        x1.extend(k.coef.reversed())
        x1.extend(k.coef[1:])
        x1.extend(flex.double(k1.m, 0.0))

        # convolve kernels
        coef = kernapply(x1, k1, circular=True)
        m = len(coef)//2
        coef = coef[m:(len(coef))]
        k = Kernel(coef=coef)

    # apply smoothing kernel
    spec = kernapply(spec, k, circular=True)
    self.spec = spec[1:fft.n_complex()]

    return

  def plot(self, sample_interval=1.0, show=True):

    import matplotlib.pyplot as plt

    # convert frequency to natural units given the sample_interval
    sample_freq = 1./sample_interval
    freq = self.freq * sample_freq

    line, = plt.semilogy(freq, self.spec)
    plt.xlabel('frequency')
    plt.ylabel('spectrum')
    if show: plt.show()
    return plt

