def adopt_init_args(obj, args):
  del args["self"]
  obj.__dict__.update(args)

class dict_with_default_value(dict):

  def __init__(self, default_value):
    self.default_value = default_value

  def __getitem__(self, key):
    try: return dict.__getitem__(self, key)
    except: pass
    val = self.default_value
    dict.__setitem__(self, key, val)
    return val

class dict_with_default_factory(dict):

  def __init__(self, default_factory):
    self.default_factory = default_factory

  def __getitem__(self, key):
    try: return dict.__getitem__(self, key)
    except: pass
    val = self.default_factory()
    dict.__setitem__(self, key, val)
    return val

def list_plus(lhs, rhs):
  return [l + r for l, r in zip(lhs, rhs)]

def list_minus(lhs, rhs):
  return [l - r for l, r in zip(lhs, rhs)]

def list_multiplies(lhs, rhs):
  return [l * r for l, r in zip(lhs, rhs)]

def list_divides(lhs, rhs):
  return [l / r for l, r in zip(lhs, rhs)]

def list_modulus(lhs, rhs):
  return [l % r for l, r in zip(lhs, rhs)]

def list_dot_product(lhs, rhs=0):
  if (rhs == 0): rhs = lhs
  result = 0
  for l, r in zip(lhs, rhs): result += l * r
  return result

def arg(a):
  "conversion of complex number from real, imag to polar angle"
  import math
  if (abs(a) == 0): return 0.
  return math.atan2(a.imag, a.real)

def polar(rho, theta):
  "conversion of complex number from polar representation to real, imag"
  import math
  return complex(rho*math.cos(theta), rho*math.sin(theta))
