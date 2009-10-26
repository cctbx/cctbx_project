from scitbx.array_family import flex
from scitbx import lstbx

def exercise_normal_equations_separating_scale_factor():
  eqs = lstbx.normal_equations_separating_scale_factor(3)
  eqs.add_equation(y_calc=1,
                   grad_y_calc=flex.double((1, 2, 3)),
                   y_obs=1,
                   weight=1)
  eqs.finalise()
  a, b = eqs.reduced_equations()
  assert a.size() == 6
  assert list(b) == [0, 0, 0]

def run():
  exercise_normal_equations_separating_scale_factor()
  print 'OK'

if __name__ == '__main__':
  run()
