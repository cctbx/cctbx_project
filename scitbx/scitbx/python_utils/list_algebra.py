def plus(lhs, rhs):
  return [l + r for l, r in zip(lhs, rhs)]

def minus(lhs, rhs):
  return [l - r for l, r in zip(lhs, rhs)]

def multiplies(lhs, rhs):
  return [l * r for l, r in zip(lhs, rhs)]

def divides(lhs, rhs):
  return [l / r for l, r in zip(lhs, rhs)]

def modulus(lhs, rhs):
  return [l % r for l, r in zip(lhs, rhs)]

def dot_product(lhs, rhs=0):
  if (rhs == 0): rhs = lhs
  result = 0
  for l, r in zip(lhs, rhs): result += l * r
  return result
