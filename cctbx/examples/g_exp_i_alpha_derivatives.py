from __future__ import absolute_import, division, print_function
import cmath
from six.moves import zip

class parameters:

  def __init__(self, alpha, g, ffp, fdp):
    self.alpha = alpha
    self.g = g
    self.ffp = ffp
    self.fdp = fdp

  def as_list(self):
    return [self.alpha, self.g, self.ffp, self.fdp]

class gradients(parameters): pass

class curvatures:

  def __init__(self, alpha_alpha, alpha_g, alpha_ffp, alpha_fdp, g_ffp, g_fdp):
    self.alpha_alpha = alpha_alpha
    self.alpha_g = alpha_g
    self.alpha_ffp = alpha_ffp
    self.alpha_fdp = alpha_fdp
    self.g_ffp = g_ffp
    self.g_fdp = g_fdp

def pack_parameters(params):
  result = []
  for p in params:
    result.extend(p.as_list())
  return result

def pack_gradients(grads):
  return pack_parameters(grads)

class g_exp_i_alpha_sum:

  def __init__(self, params):
    self.params = params

  def f(self):
    "Mathematica: f=g Exp[I alpha]"
    result = 0
    for p in self.params:
      result += p.g * (p.ffp + 1j*p.fdp) * cmath.exp(1j*p.alpha)
    return result

  def d_params(self):
    "Mathematica: D[f,g]; D[f,ffp]; D[f,fdp]; D[f,alpha]"
    result = []
    for p in self.params:
      eja = cmath.exp(1j*p.alpha)
      ffpifdp = p.ffp + 1j*p.fdp
      result.append(gradients(
        alpha=p.g*ffpifdp*1j*eja,
        g=ffpifdp*eja,
        ffp=p.g*eja,
        fdp=p.g*1j*eja))
    return result

  def d2_params(self):
    """Mathematica:
         D[f,alpha,alpha]; D[f,alpha,g]; D[f,alpha,ffp]; D[f,alpha,fdp]
         D[f,g,    alpha]; D[f,g,    g]; D[f,g,    ffp]; D[f,g,    fdp]
         D[f,ffp,  alpha]; D[f,ffp,  g]; D[f,ffp,  ffp]; D[f,ffp,  fdp]
         D[f,fdp,  alpha]; D[f,fdp,  g]; D[f,fdp,  ffp]; D[f,fdp,  fdp]
       Zeros/non-zeros in upper triangle:
         . . . .
           0 . .
             0 0
               0
    """
    result = []
    for p in self.params:
      eja = cmath.exp(1j*p.alpha)
      ffpifdp = p.ffp + 1j*p.fdp
      result.append(curvatures(
        alpha_alpha=-p.g*ffpifdp*eja,
        alpha_g=ffpifdp*1j*eja,
        alpha_ffp=p.g*1j*eja,
        alpha_fdp=-p.g*eja,
        g_ffp=eja,
        g_fdp=1j*eja))
    return result

  def d_target_d_params(self, target):
    "Rule for derivatives of sum of roots of unity."
    result = []
    da, db = target.da(), target.db()
    for d_params in self.d_params():
      result.append(gradients(*[da * d.real + db * d.imag
        for d in d_params.as_list()]))
    return result

  def d2_target_d_params(self, target):
    "Product rule applied to da * d.real + db * d.imag."
    da, db = target.da(), target.db()
    daa, dbb, dab = target.daa(), target.dbb(), target.dab()
    d = self.d_params()
    d2 = self.d2_params()
    i4 = 0
    for di,d2i in zip(d, d2):
      for ixi,dix in enumerate(di.as_list()):
        row = []
        for dj in d:
          for ixj,djx in enumerate(dj.as_list()):
            sum = daa * dix.real * djx.real \
                + dbb * dix.imag * djx.imag \
                + dab * (dix.real * djx.imag + dix.imag * djx.real)
            row.append(sum)
        if (ixi == 0): # (0,0)
          row[i4] += da * d2i.alpha_alpha.real \
                   + db * d2i.alpha_alpha.imag
        if (ixi == 0 or ixi == 1): # (0,1) or (1,0)
          row[i4+1-ixi] += da * d2i.alpha_g.real \
                         + db * d2i.alpha_g.imag
        if (ixi == 0 or ixi == 2): # (0,2) or (2,0)
          row[i4+2-ixi] += da * d2i.alpha_ffp.real \
                         + db * d2i.alpha_ffp.imag
        if (ixi == 0 or ixi == 3): # (0,3) or (3,0)
          row[i4+3-ixi] += da * d2i.alpha_fdp.real \
                         + db * d2i.alpha_fdp.imag
        if (ixi == 1 or ixi == 2): # (1,2) or (2,1)
          row[i4+3-ixi] += da * d2i.g_ffp.real \
                         + db * d2i.g_ffp.imag
        if (ixi == 1 or ixi == 3): # (1,3) or (3,1)
          row[i4+4-ixi] += da * d2i.g_fdp.real \
                         + db * d2i.g_fdp.imag
        yield row
      i4 += 4
