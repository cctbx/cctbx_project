from smtbx import refinement
from iotbx.shelx.from_ins import from_ins
from iotbx.reflection_file_reader import any_reflection_file
from cctbx import xray

tests_data_dir = '../smtbx_tests_data/refinement'

def exercise_disorder(name):
  xs0 = from_ins('%s/%s.res' % (tests_data_dir, name))
  data_file = any_reflection_file('hklf4=%s/%s.hkl' % (tests_data_dir, name))
  f_obs_sqr = data_file.as_miller_arrays(xs0.crystal_symmetry())[0]
  xs = refinement.tests.shaked_structure(xs0,
                                         thermal_shift=0.01, # in %
                                         site_shift=0.02)
  for a in xs.scatterers():
    a.flags.set_grad_site(True)
    if a.flags.use_u_iso(): a.flags.set_grad_u_iso(True)
    if a.flags.use_u_aniso(): a.flags.set_grad_u_aniso(True)
    assert not( a.flags.use_u_iso() and a.flags.grad_u_aniso() )
    assert not( a.flags.use_u_aniso() and a.flags.grad_u_iso() )
  ls = xray.unified_least_squares_residual(f_obs_sqr)
  minimisation = refinement.minimization.lbfgs(
    target_functor=ls,
    xray_structure=xs,
    cos_sin_table=True,
  )
  
  

def run():
  exercise_disorder(name='6-ring-twist-disorder')
  print 'OK'

if __name__ == '__main__':
  run()




