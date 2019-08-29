from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import time
from six.moves import range

def block(data_size, n_repeats):
  data = flex.random_double(size=data_size)
  permutation = flex.random_permutation(data.size())
  for use_pointers, use_iterators_range in [(False, 2), (True, 3)]:
    print("  use_pointers =", use_pointers)
    for use_iterators in range(use_iterators_range):
      t0 = time.time()
      result = flex.integer_offsets_vs_pointers(
        data, permutation, 0, use_pointers, use_iterators)
      t_overhead = time.time() - t0
      t0 = time.time()
      result = flex.integer_offsets_vs_pointers(
        data,
        permutation,
        n_repeats,
        use_pointers,
        use_iterators)
      print("    use_iterators =", use_iterators, \
            "time = %.2f s" % (time.time() - t0 - t_overhead), \
            "overhead = %.2f s" % t_overhead, \
            "(%-7s %.6g)" % result)

def run():
  for data_size in [1000, 10000, 100000, 1000000, 10000000]:
    n_repeats = 30000000 // data_size
    print("data_size =", data_size, "n_repeats =", n_repeats)
    block(data_size, n_repeats)
  print("OK")

if (__name__ == "__main__"):
  run()
