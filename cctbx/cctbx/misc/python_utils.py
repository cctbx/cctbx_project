def adopt_init_args(obj, args):
  del args["self"]
  obj.__dict__.update(args)

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
