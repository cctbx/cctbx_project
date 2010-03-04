from scitbx.array_family import flex

class two_stats(object):
  """ The clustering of data sorted in reverse order as data[:n] U data[n:].
  Instances feature:
    - self.cut: the index n
    - self.highest_mean: the mean of the strongest cluster data[:n]
    - self.lowest_mean: the mean of the weakest cluster data[n:]
  """

  def __init__(self, data):
    self.n_data = n = len(data)
    if n == 0: return 0
    if n == 1: return 1
    cut = None
    new_cut = n//2
    while new_cut != cut:
      cut = new_cut
      hi, lo = self.statistics(data[:cut]), self.statistics(data[cut:])
      for i, x in enumerate(data):
        if x >= hi: continue
        if abs(x - lo) <= abs(x - hi):
          new_cut = i
          break
    self.cut = cut
    self.highest_mean = hi
    self.lowest_mean = lo

  def __str__(self):
    h, l = [ "%.3g" % x for x in (self.highest_mean, self.lowest_mean) ]
    dash_h, dash_l = [ "-"*len(x) for x in (h, l) ]
    c = "%i" % self.cut
    buf_c = " "*len(c)
    return (  "|-%s-|%s|-%s-|%i\n" % (dash_h, c, dash_l, self.n_data)
            + "  %s <%s> %s" % (h, buf_c, l))


class two_means(two_stats):
  statistics = staticmethod(flex.mean)

class two_medians(two_stats):
  statistics = staticmethod(flex.median)
