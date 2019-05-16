from __future__ import absolute_import, division, print_function
import cmath
import math
from six.moves import zip

class least_squares:

  def __init__(self, obs, calc):
    self.obs = obs
    self.calc = calc
    a, b = self.calc.real, self.calc.imag
    self.abs_calc = math.sqrt(a**2 + b**2)
    self.delta = self.obs - self.abs_calc

  def f(self):
    "Mathematica: f=(obs-Sqrt[a^2+b^2])^2"
    return self.delta**2

  def da(self):
    "Mathematica: D[f,a]"
    if (self.abs_calc == 0): return 0
    return -2 * self.delta * self.calc.real / self.abs_calc

  def db(self):
    "Mathematica: D[f,b]"
    if (self.abs_calc == 0): return 0
    return -2 * self.delta * self.calc.imag / self.abs_calc

  def daa(self):
    "Mathematica: FortranForm[FullSimplify[D[f,a,a]]]"
    ac = self.abs_calc
    if (ac == 0):
      if (self.obs == 0): return 2
      return -1.e160
    return 2 - (2*self.calc.imag**2*self.obs)/ac/ac/ac

  def dbb(self):
    "Mathematica: FortranForm[FullSimplify[D[f,b,b]]]"
    ac = self.abs_calc
    if (ac == 0):
      if (self.obs == 0): return 2
      return -1.e160
    return 2 - (2*self.calc.real**2*self.obs)/ac/ac/ac

  def dab(self):
    "Mathematica: FortranForm[FullSimplify[D[f,a,b]]]"
    ac = self.abs_calc
    if (ac == 0):
      if (self.obs == 0): return 0
      return 1.e160
    return (2*self.calc.real*self.calc.imag*self.obs)/ac/ac/ac

class exp_i_alpha_sum:

  def __init__(self, alphas):
    self.alphas = alphas

  def f(self):
    "Mathematica: f=Exp[I alpha]"
    result = 0
    for alpha in self.alphas:
      result += cmath.exp(1j*alpha)
    return result

  def d_alphas(self):
    "Mathematica: D[f,alpha]"
    return [1j*cmath.exp(1j*alpha) for alpha in self.alphas]

  def d2_alphas(self):
    "Mathematica: D[f,alpha,alpha]"
    return [-cmath.exp(1j*alpha) for alpha in self.alphas]

  def d_target_d_alphas(self, target):
    "Rule for derivatives of sum of roots of unity."
    da, db = target.da(), target.db()
    return [da * d.real + db * d.imag for d in self.d_alphas()]

  def d2_target_d_alphas(self, target):
    "Product rule applied to da * d.real + db * d.imag."
    result = []
    da, db = target.da(), target.db()
    daa, dbb, dab = target.daa(), target.dbb(), target.dab()
    d = self.d_alphas()
    d2 = self.d2_alphas()
    for di,d2i in zip(d, d2):
      row = []
      for dj in d:
        sum = daa * di.real * dj.real \
            + dbb * di.imag * dj.imag \
            + dab * (di.real * dj.imag + di.imag * dj.real)
        if (di is dj):
          sum += da * d2i.real + db * d2i.imag
        row.append(sum)
      result.append(row)
    return result
