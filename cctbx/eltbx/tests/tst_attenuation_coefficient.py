from __future__ import division

def tst_for_z(z):
  from cctbx.eltbx import attenuation_coefficient

  eps = 1e-7

  # Get the table
  table = attenuation_coefficient.get_table(z)

  # Get list of energies and coefficients
  energy = table.energy()[1:-1]
  mu_rho = table.mu_rho()[1:-1]
  mu_en_rho = table.mu_en_rho()[1:-1]

  # Check values at measured energies match
  for i in range(len(energy)):
    e = energy[i]
    mr = mu_rho[i]
    mer = mu_en_rho[i]

    if i < len(energy) - 1:
      if abs(e - energy[i+1]) < eps:
        continue

    mr2 = table.mu_rho_at_ev(e * 1000000.0)
    mer2 = table.mu_en_rho_at_ev(e * 1000000.0)
    assert(abs(mr - mr2) <= eps)
    assert(abs(mer - mer2) <= eps)

  print 'OK'

def run():
  from cctbx.eltbx import attenuation_coefficient

  for z in attenuation_coefficient.nist_elements().atomic_number_list():
    tst_for_z(z)

if __name__ == '__main__':
  run()
