from __future__ import division
from scitbx.array_family import flex

def I_and_G_base_estimate(data):
  """Estimate I and G, based on simple weighted averages, nothing fancy"""
  Nframes = flex.max(data.frame) + 1
  Nmiller = flex.max(data.miller) + 1
  print "%d frames, %d Miller indices"%(Nframes, Nmiller)

  G,G_visited = data.estimate_G(Nframes)

  I,I_visited = data.estimate_I(Nmiller)
  return I,I_visited,G,G_visited

def plot_it(fit, sim, mode=None):
  sub_fit = []
  sub_sim = []
  for i in xrange(len(fit)):
    if sim[i] > 0:
      sub_fit.append(fit[i])
      sub_sim.append(sim[i])
  from matplotlib import pyplot as plt
  if mode=="I":
    plt.loglog(sub_sim, sub_fit, "b.")
    plt.title("Log plot of fitted intensity vs. simulated intensity")
  else:
    plt.plot(sub_sim, sub_fit, "b.")
    if max(sub_sim)>1E8:
      plt.axes().set_xlim(0,5E7)
    if max(sub_fit)>1E8:
      plt.axes().set_ylim(0,5E7)
  if mode=="B":
    plt.axes().set_ylim(-12,12)
    plt.title("Fitted B-factor vs. simulated B-factor")

  plt.show()

def show_correlation(A, B, selection, message):
  data_A = A.select(selection==1)
  data_B = B.select(selection==1)
  LC = flex.linear_correlation(data_A,data_B)
  print message,LC.coefficient(),"on %d values"%LC.n()
