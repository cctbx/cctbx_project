from scitbx import matrix

def potential_energy(sites, wells, A_T, J_T_inv):
  result = 0
  for s, w in zip(sites, wells):
    result += (A_T * s - J_T_inv * A_T * w).dot()
  return result

def potential_f_ext_pivot_at_origin(sites, wells, A_T, J_T_inv):
  f_cart = [2 * (A_T * s - J_T_inv * A_T * w) for s, w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force in zip(sites, f_cart):
    f += force
    nc += (A_T * s).cross(force)
  return matrix.col((nc, f)).resolve_partitions()

def potential_energy_no_align(sites, wells, J_T_inv):
  result = 0
  for s, w in zip(sites, wells):
    result += (s - J_T_inv * w).dot()
  return result

def potential_f_ext_no_align_pivot_at_origin(sites, wells, J_T_inv):
  f_cart = [2 * (s - J_T_inv * w) for s, w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force in zip(sites, f_cart):
    f += force
    nc += s.cross(force)
  return matrix.col((nc, f)).resolve_partitions()
