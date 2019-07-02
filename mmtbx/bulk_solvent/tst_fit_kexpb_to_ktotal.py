from __future__ import absolute_import, division, print_function
import mmtbx.bulk_solvent
from scitbx.array_family import flex
import math, time
from six.moves import range

def run(nref, k, b, k_start, b_start):
  data = flex.double()
  ss = flex.double()
  for ss_ in range(1, nref):
    ss_ = ss_/nref*1.
    y = k*math.exp(-b*ss_)
    data.append(y)
    ss.append(ss_)
  t0 = time.time()
  r = mmtbx.bulk_solvent.fit_k_exp_b_to_k_total(data, ss, k_start, b_start)
  return list(r), "Time: %6.4f"%(time.time()-t0)

if (__name__ == "__main__"):
  for i in [(10000, 1, 10, 1, 10),
            (10000,10,100,10,100),
            (10000, 5, 10, 5, 10),
            (10000, 0,  0, 0,  0),
            (10000, 90, -10, -10, 10) # outside convergence well
            ]:
    r = run(nref=i[0], k=i[1], b=i[2], k_start=i[3], b_start=i[4])
    print(r)
