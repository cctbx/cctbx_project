from cctbx.eltbx import tiny_pse
from scitbx.math import continued_fraction
from boost import rational

class formula(object):

  def __init__(self, count_of_element):
    self.count_of_element = count_of_element

  def sorted_as_c_h_then_by_increasing_atomic_number(self):
    head = []
    for elt in ('C', 'H'):
      n = self.count_of_element.get(elt)
      if n: head.append((elt, n))
    tail = [ (tiny_pse.table(elt).atomic_number(), (elt, n))
             for elt, n in self.count_of_element.iteritems()
             if elt not in ('C', 'H') ]
    tail.sort()
    self.element_count_pairs = head + [ item[-1] for item in tail ]
    return self

  def __str__(self):
    result = []
    for elt, n in self.element_count_pairs:
      r = continued_fraction.from_real(n, eps=1e-5).as_rational()
      if r.denominator() == 1:
        if r.numerator() > 1: m = "%i" % r.numerator()
        else: m = ""
      elif r.denominator() < 10: m = "%i/%i" % (r.numerator(), r.denominator())
      else: m = "%.5f" % n
      result.append("%s%s" % (elt, m))
    return ' '.join(result)
