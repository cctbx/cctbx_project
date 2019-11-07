from __future__ import absolute_import, division, print_function
from scitbx.math import clustering
from scitbx.array_family import flex

def exercise_two_means_clustering():
  data = flex.double((0.21947204867684716, 0.21947204867684716,
                      0.21947204867684716, 0.21947204867684714))
  c = clustering.two_means(data)
  assert c.cut == 4 # precision limit

  data = flex.double((1, 1, 1, 1))
  c = clustering.two_means(data)
  assert c.cut == 4

  data = flex.double((10, 4, 3, 2, 2, 1))
  c = clustering.two_means(data)
  assert c.cut == 1

  data = flex.double((10, 8, 5, 4, 4, 3))
  assert clustering.two_means(data).cut == 2

  data = flex.double((10, 9, 9, 8, 8, 8))
  assert clustering.two_means(data).cut == 3

  data = flex.double((1000, 900, 800, 500, 200, 200, 100))
  assert clustering.two_means(data).cut == 3

  data = flex.double((0.91, 0.29, 0.29, 0.27, 0.27, 0.27, 0.27))
  assert clustering.two_means(data).cut == 1

  data = flex.double((0.998, 0.449, 0.152, 0.152, 0.151))
  clusters = clustering.two_means(data)
  assert clusters.cut == 2

  data = flex.double((0.998, 0.449, 0.152))
  clusters = clustering.two_means(data)
  assert clusters.cut == 1

def exercise_two_medians_clustering():
  data = flex.double((0.998, 0.449, 0.152, 0.152, 0.151))
  clusters = clustering.two_medians(data)
  assert clusters.cut == 2

def run():
  exercise_two_means_clustering()
  exercise_two_medians_clustering()
  print('OK')

if __name__ == '__main__':
  run()
