import scitbx.random
from libtbx.test_utils import approx_equal

def exercise_distributions():
  n = scitbx.random.ext.normal_distribution()
  assert (n.mean(), n.sigma()) == (0, 1)
  n = scitbx.random.ext.normal_distribution(mean=5, sigma=10)
  assert (n.mean(), n.sigma()) == (5, 10)

def exercise_variate_generators():
  for i in range(10):
    scitbx.random.set_random_seed(0)
    g = scitbx.random.normal_variate_generator()
    assert approx_equal(g(), -1.2780081289048213)
    assert approx_equal(g(10),
      (-0.40474189234755492, -0.41845505596083288,
       -1.8825790263067721, -1.5779112018107659,
       -1.1888174422378859, -1.8619619179878537,
       -0.53946818661388318, -1.2400941724410812,
       0.64511959841907285, -0.59934120033270688))

def run():
  exercise_distributions()
  exercise_variate_generators()
  print "OK"

if __name__ == '__main__':
  run()
