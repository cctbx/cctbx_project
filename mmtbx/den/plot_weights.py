
def _plot_weights (figure, gamma, weight, z) :
  from matplotlib import cm
  assert (len(z[0]) == len(gamma)) and (len(z) == len(weight))
  p = figure.add_subplot(111)
  p.set_position([0.1,0.1,0.85,0.85])
  cset = p.contourf(z, 20, cmap=cm.YlGnBu, interpolation='bilinear')
  p.contour(z, 20, colors=[(0.5,0.5,0.5)], linewidth=1)
  figure.colorbar(cset, ax=p)
  p.set_xticks(range(len(gamma)))
  p.set_yticks(range(len(weight)))
  p.set_xticklabels([ "%g" % x for x in gamma ])
  p.set_yticklabels([ "%g" % x for x in weight ])
  p.set_xlabel("gamma")
  p.set_ylabel("weight")
  p.set_title("DEN optimization (R-free)")
  return p

def plot_weights_pyplot (gamma, weight, z) :
  from matplotlib import pyplot as plt
  figure = plt.figure()
  _plot_weights(figure, gamma, weight, z)
  plt.show()

def exercise () :
  grid_results = [
    (0.0, 0.5, 0.2330),
    (0.0, 1.0, 0.2313),
    (0.0, 3.0, 0.2259),
    (0.0, 5.0, 0.2295),
    (0.0, 10.0, 0.2252),
    (0.0, 25.0, 0.2334),
    (0.0, 50.0, 0.2329),
    (0.0, 100.0, 0.2434),
    (0.0, 200.0, 0.2612),
    (0.2, 0.5, 0.2343),
    (0.2, 1.0, 0.2346),
    (0.2, 3.0, 0.2331),
    (0.2, 5.0, 0.2328),
    (0.2, 10.0, 0.2320),
    (0.2, 25.0, 0.2280),
    (0.2, 50.0, 0.2313),
    (0.2, 100.0, 0.2379),
    (0.2, 200.0, 0.2555),
    (0.4, 0.5, 0.2352),
    (0.4, 1.0, 0.2353),
    (0.4, 3.0, 0.2328),
    (0.4, 5.0, 0.2307),
    (0.4, 10.0, 0.2272),
    (0.4, 25.0, 0.2308),
    (0.4, 50.0, 0.2300),
    (0.4, 100.0, 0.2359),
    (0.4, 200.0, 0.2427),
    (0.6, 0.5, 0.2357),
    (0.6, 1.0, 0.2300),
    (0.6, 3.0, 0.2352),
    (0.6, 5.0, 0.2356),
    (0.6, 10.0, 0.2347),
    (0.6, 25.0, 0.2249),
    (0.6, 50.0, 0.2344),
    (0.6, 100.0, 0.2332),
    (0.6, 200.0, 0.2407),
    (0.8, 0.5, 0.2311),
    (0.8, 1.0, 0.2362),
    (0.8, 3.0, 0.2288),
    (0.8, 5.0, 0.2335),
    (0.8, 10.0, 0.2337),
    (0.8, 25.0, 0.2295),
    (0.8, 50.0, 0.2283),
    (0.8, 100.0, 0.2276),
    (0.8, 200.0, 0.2310),
    (1.0, 0.5, 0.2367),
    (1.0, 1.0, 0.2293),
    (1.0, 3.0, 0.2289),
    (1.0, 5.0, 0.2373),
    (1.0, 10.0, 0.2343),
    (1.0, 25.0, 0.2322),
    (1.0, 50.0, 0.2326),
    (1.0, 100.0, 0.2269),
    (1.0, 200.0, 0.2328),
  ]
  gamma = sorted(list(set([ x[0] for x in grid_results ])))
  weight = sorted(list(set([ x[1] for x in grid_results ])))
  values = [ [] for x in weight ]
  n = len(weight)
  for i in range(len(weight)) :
    for j in range(len(gamma)) :
      k = i + j * n
      values[i].append(grid_results[k][2])
  plot_weights_pyplot(gamma, weight, values)

if (__name__ == "__main__") :
  exercise()
