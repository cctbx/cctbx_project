import cmath

class parameters:

  def __init__(self, g, alpha):
    self.g = g
    self.alpha = alpha

  def as_list(self):
    return [self.g, self.alpha]

class gradients(parameters): pass

class curvatures:

  def __init__(self, g_alpha, alpha_alpha):
    self.g_alpha = g_alpha
    self.alpha_alpha = alpha_alpha

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
      result += p.g * cmath.exp(1j*p.alpha)
    return result

  def d_params(self):
    "Mathematica: D[f,g]; D[f,alpha]"
    result = []
    for p in self.params:
      eja = cmath.exp(1j*p.alpha)
      result.append(gradients(g=eja, alpha=p.g*1j*eja))
    return result

  def d2_params(self):
    "Mathematica: D[f,g,g]; D[f,g,alpha]; D[f,alpha,alpha]"
    result = []
    for p in self.params:
      eja = cmath.exp(1j*p.alpha)
      result.append(curvatures(g_alpha=1j*eja, alpha_alpha=-p.g*eja))
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
    result = []
    da, db = target.da(), target.db()
    daa, dbb, dab = target.daa(), target.dbb(), target.dab()
    d = self.d_params()
    d2 = self.d2_params()
    for di,d2i in zip(d, d2):
      for ixi,dix in enumerate(di.as_list()):
        row = []
        for dj in d:
          for ixj,djx in enumerate(dj.as_list()):
            sum = daa * dix.real * djx.real \
                + dbb * dix.imag * djx.imag \
                + dab * (dix.real * djx.imag + dix.imag * djx.real)
            if (di is dj):
              if (ixi == 1 and ixj == 1):
                sum += da * d2i.alpha_alpha.real \
                     + db * d2i.alpha_alpha.imag
              elif (ixi == 1 or ixj == 1):
                sum += da * d2i.g_alpha.real \
                     + db * d2i.g_alpha.imag
            row.append(sum)
        result.append(row)
    return result
