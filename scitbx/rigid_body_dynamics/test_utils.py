from scitbx import matrix

def potential_energy(sites, wells, A, J, AJA_tree=None):
  result = 0
  AJA = A.T_inv * J.T * A.T
  if (AJA_tree is not None): AJA = AJA * AJA_tree
  for s, w in zip(sites, wells):
    result += (AJA * s - w).dot()
  return result

def potential_f_ext_pivot_at_origin(sites, wells, A, J, AJA_tree=None):
  AJA = A.T_inv * J.T * A.T
  if (AJA_tree is not None): AJA = AJA * AJA_tree
  f_cart_ff = [2 * (AJA * s - w) for s, w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force_ff in zip(sites, f_cart_ff):
    f += force_ff
    nc += (AJA * s).cross(force_ff)
  return matrix.col((nc, f)).resolve_partitions()

def potential_energy_no_align(sites, wells, J):
  result = 0
  for s, w in zip(sites, wells):
    result += (J.T * s - w).dot()
  return result

def potential_f_ext_no_align_pivot_at_origin(sites, wells, J):
  f_cart_ff = [2 * (J.T * s - w) for s, w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force_ff in zip(sites, f_cart_ff):
    f += force_ff
    nc += (J.T * s).cross(force_ff)
  return matrix.col((nc, f)).resolve_partitions()
