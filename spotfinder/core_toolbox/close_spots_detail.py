from __future__ import absolute_import, division, print_function
from six.moves import range
import math
from spotfinder.array_family import flex

'''Separate module for special treatment of close spots.
'''

#extraordinary procedure, when many spots are close.  Include closely
#  spaced spots in autoindexing, if it appears that they truly
#  reflect lattice spacing.

# a poor man's Fourier transform of the image, to determine
#   if close spots form a regular pattern

class NearNeighborVectors:

  def __init__(self, ave_area, query_centers, fstats, ann_adaptor):
    self.plot = {}
    self.pmax = []

    plot_to_spot_ratio = 20.
    significant_vector_fraction_of_highest = 0.2
    significant_vector_fraction_of_recent_child = 0.05

    # make plot with area 20*average spot area (an arbitrary cutoff)
    self.deltaxylimit = int(math.sqrt(plot_to_spot_ratio*ave_area)/2.)+1
    for x in range(-self.deltaxylimit,self.deltaxylimit+1):
      self.plot[x]=flex.int(2*self.deltaxylimit+1)

    for d in range(len(query_centers)//2):
      vector = (int( fstats.master[ann_adaptor.nn[d]].max_pxl_x()-query_centers[2*d] ),
                int( fstats.master[ann_adaptor.nn[d]].max_pxl_y()-query_centers[2*d+1]))
      if abs(vector[0])<=self.deltaxylimit and \
         abs(vector[1])<=self.deltaxylimit:
        self.plot[vector[0]][self.deltaxylimit+vector[1]]+=1
        self.plot[-vector[0]][self.deltaxylimit-vector[1]]+=1

    # get the maximum value of this "Fourier" plot
    allx = ()
    for x in range(-self.deltaxylimit+1,self.deltaxylimit):
      allx=allx+tuple(self.plot[x])
    maxp = max(allx) #"the largest peak is",maxp

    # Search for all local maxima in the plot. Only look in positive hemisphere.
    for locx in range(-self.deltaxylimit+1,self.deltaxylimit):
      for locy in range(-self.deltaxylimit+1,self.deltaxylimit):
        if locy>0 or (locy==0 and locx>0):  #primary hemisphere
          tlocy = self.deltaxylimit+locy
          refval = self.plot[locx][tlocy]
          #(Two important parameters subject to adjustment)
          if refval > significant_vector_fraction_of_highest*maxp and \
             refval > significant_vector_fraction_of_recent_child*\
                      fstats['N_%s'%fstats.most_recent_child()] and \
             refval >= self.plot[locx][tlocy+1] and \
             refval >= self.plot[locx+1][tlocy+1] and \
             refval >= self.plot[locx+1][tlocy] and \
             refval >= self.plot[locx+1][tlocy-1] and \
             refval > self.plot[locx][tlocy-1] and \
             refval > self.plot[locx-1][tlocy-1] and \
             refval > self.plot[locx-1][tlocy] and \
             refval > self.plot[locx-1][tlocy+1]:
             self.pmax.append((locx,locy,refval))

  def show_maxima(self):
    print("all maxima:",self.pmax)

  def show_vector_map(self):
     for x in range(-self.deltaxylimit,self.deltaxylimit+1):
      for i,y in enumerate(self.plot[x]):
        if x==0 and i==len(self.plot[x])//2:
          print("   X", end=' ')
        else:
          print("%4d"%y, end=' ')
      print()

  def vectors(self):
    return self.pmax
