from scitbx import matrix

def potential_energy(sites, wells, A, J):
  result = 0
  for s, w in zip(sites, wells):
    result += (A.T_inv * J.T * A.T * s - w).dot()
  return result

def potential_f_ext_pivot_at_origin(sites, wells, A, J):
  AJ = A.T_inv * J.T
  AJA = AJ * A.T
  f_cart_ff = [2 * (AJA * s - w) for s, w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force_ff in zip(sites, f_cart_ff):
    force_bf = (AJ.r).transpose() * force_ff
    f += force_bf
    nc += (A.T * s).cross(force_bf)
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
    force_bf = J.T.r.transpose() * force_ff
    f += force_bf
    nc += s.cross(force_bf)
  return matrix.col((nc, f)).resolve_partitions()
