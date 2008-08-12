import sys, math, re, cStringIO
from libtbx.utils import show_times_at_exit
from libtbx.test_utils import approx_equal

class cosine_integrals(object):

  def __init__(self, i_range, n):
    self.i_range = xrange(*i_range)
    def func(i):
      s = 0
      n = 1e5/(i*i % 5 + 1)
      h = math.pi/2/n
      for j in xrange(int(n)):
        s += math.cos(i*j*h)
      s = abs(s*h)
      print "%i: S=%.2f" % (i,s)
      return s
    self.func = func

  def exercise(self, n_workers):
    import libtbx.symmetric_multi_processing as smp
    if (not smp.is_available):
      print "Skipped!"
      return
    show_times_at_exit()
    result_pat = re.compile("(\d+): S=(\d\.\d\d)")
    ref_results = [ abs(math.sin(i*math.pi/2)/i) for i in self.i_range ]
    ref_printout = zip(self.i_range, [ "%.2f" % s for s in ref_results ])
    if 0:
      debug = True
      timeout = 0.1
    else:
      debug = False
      timeout = 0.001
    f = smp.parallelized_function(self.func, n_workers=n_workers,
                                  stdout=cStringIO.StringIO(),
                                  timeout=timeout,
                                  debug=debug)
    results = list(f(self.i_range))
    if 0:
      print results
      print f.printout.getvalue()
    assert approx_equal(results, ref_results, eps=1e-3)
    matches = [ result_pat.search(li)
                for li in f.stdout.getvalue().strip().split('\n') ]
    printout = [ (int(m.group(1)), m.group(2)) for m in matches ]
    printout.sort()
    assert printout == ref_printout
    sys.stdout.flush()

def run():
  from libtbx import option_parser
  parser = (option_parser.option_parser()
            .option('--test_case',
                    action='store',
                    default='cosine_integrals(i_range=(1,33), n=1e5)')
            .option('--workers',
                    action='store',
                    type='int',
                    default=2)
            )
  co = parser.process(sys.argv[1:], max_nargs=0).options
  test_case = eval(co.test_case)
  test_case.exercise(n_workers=co.workers)
  print "OK"

if __name__ == '__main__':
  run()
