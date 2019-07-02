from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex

def I_and_G_base_estimate(data,params=None):
  """Estimate I and G, based on simple weighted averages, nothing fancy"""
  Nframes = flex.max(data.frame) + 1
  # worth some explanation.  For a half_data_set on 1000 images, one of the
  # halves will only have a data.frame max of 999. An apparent discrepancy,
  # not functionally important.
  Nmiller = flex.max(data.miller) + 1

  # figure out if d_max, d_min were requested
  inv_d_max_sq = 0.0
  inv_d_min_sq = 0.0
  if params is not None:
    if params.d_max is not None: inv_d_max_sq = (1./params.d_max)*(1./params.d_max)
    inv_d_min_sq = (1./params.d_min)*(1./params.d_min)

  G,G_visited = data.estimate_G(Nframes,inv_d_max_sq,inv_d_min_sq)
  I,I_visited = data.estimate_I(Nmiller,inv_d_max_sq,inv_d_min_sq)
  return I,I_visited,G,G_visited

def plot_it(fit, sim, mode=None):
  sub_fit = []
  sub_sim = []
  for i in range(len(fit)):
    if sim[i] > 0:
      sub_fit.append(fit[i])
      sub_sim.append(sim[i])
  from matplotlib import pyplot as plt
  if mode=="I":
    plt.loglog(sub_sim, sub_fit, "b,")
    plt.title("Log plot of fitted intensity vs. simulated intensity")
  else:
    plt.plot(sub_sim, sub_fit, "b.")
    if max(sub_sim)>1E8:
      plt.axes().set_xlim(0,5E7)
    if max(sub_fit)>1E8:
      plt.axes().set_ylim(0,5E7)
  if mode=="G":
    plt.title("Fitted scale-factor G vs. simulated scale-factor")
  elif mode=="B":
    plt.axes().set_ylim(-12,12)
    plt.title("Fitted B-factor vs. simulated B-factor")

  plt.show()

def show_correlation(A, B, selection, message):
  data_A = A.select(selection==1)
  data_B = B.select(selection==1)
  LC = flex.linear_correlation(data_A,data_B)
  print(message,LC.coefficient(),"on %d values"%LC.n())

def show_histogram(data,title):
  from matplotlib import pyplot as plt

  nbins = 40
  n,bins,patches = plt.hist(data,
    nbins, normed=0, facecolor="orange", alpha=0.75)

  plt.xlabel("Rotational angle (degrees)")
  plt.title(title)
  #plt.axis([0,maxangle,0,1450])
  #plt.plot([median],[1375],"b|")
  #plt.plot([rmsd],[1375],"b|")
  plt.show()
