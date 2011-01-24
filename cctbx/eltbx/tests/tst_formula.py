from cctbx.eltbx.formula import formula

def exercise():
  f = formula({'C': 10, 'H': 30, 'N': 2, 'O': 1, 'B': 3, 'Cl': 2, 'Rh': 1}
              ).sorted_as_c_h_then_by_increasing_atomic_number()
  assert str(f) == 'C10 H30 B3 N2 O Cl2 Rh'

  f = formula({'Li': 1, 'Br': 2, 'O': 3}
              ).sorted_as_c_h_then_by_increasing_atomic_number()
  assert str(f) == 'Li O3 Br2'

  f = formula({'Li': 1, 'Br': 2, 'C': 3}
              ).sorted_as_c_h_then_by_increasing_atomic_number()
  assert str(f) == 'C3 Li Br2'

  f = formula({'C': 4.33333, 'B': 2.66666, 'O': 0.123454, 'Li': 3}
              ).sorted_as_c_h_then_by_increasing_atomic_number()
  assert str(f) == 'C13/3 Li3 B8/3 O0.12345'

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
