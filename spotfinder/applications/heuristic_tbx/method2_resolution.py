import math
from libtbx.table_utils import Spreadsheet,Formula
from spotfinder.core_toolbox import bin_populations

'''Support for determining the resolution limit from a group of
Bragg spots.  This module provides tools for organizing resolution
bin information into a spreadsheet; the caller is responsible for use
of the information'''

class ResolutionShells(Spreadsheet):
  def __init__(self,spots,targetResolution,fractionCalculator):
    #Total rows by formula
    Total_rows = int(math.pow(
       targetResolution/spots[len(spots)-1],3.))+1
    Spreadsheet.__init__(self,rows=Total_rows)
    self.addColumn('Limit')
    self.addColumn('Fract')
    self.addColumn('Population',0)
    self.addColumn('adjustPop',
                   Formula('self.Population[%row] / self.Fract[%row]'))

    # the outer resolution limit of the lowest-resolution shell
    self.Limit[0] = targetResolution
    self.Fract[0] = 1.0

    def BinRes(shellnumber):
      return self.Limit[0]/math.pow(shellnumber+1,1.0/3.0)

    #this little loop consumes 1.2 seconds of CPU time for a 1400-row spreadsheet
    for xrow in xrange(1,Total_rows):
      self.Limit[xrow] = BinRes(xrow)
      self.Fract[xrow] = fractionCalculator(self.Limit[xrow])

    bp = bin_populations(spots,self.Limit[0])
    for c in xrange(min(Total_rows,len(bp))): #reconcile inconsistent bin counts
        self.Population[c]=bp[c]

    self.Limit.format = "%.2f"
    self.Fract.format = "%.2f"
    self.Population.format = "%7d"
    self.adjustPop.format = "%7.1f"

  def show(self,default=['Limit','Population','Fract']):
    self.printTable(default)

  def vetter(self,cutoff,last):
    redcutoff = 10
    """The problem:  we are given a list of shell populations, pop.
       We are told we are only interested from index 0 thru last.
       We have to find a new slice from 0 to newlast such that there is
       no contiguous stretch of populations < cutoff, which is longer than 10.
    """
    pop = []
    for x in xrange(0,last+1):
      if self.Population[x]>cutoff: pop.append(1)
      else: pop.append(0)

    requalify=[0,]
    for x in xrange(1,last+1):
      if pop[x]==0: requalify.append(requalify[x-1]+1)
      else:
         requalify.append(0)
      if requalify[x] == redcutoff:
        break

    pop = pop[0:x+1]
    for x in xrange(len(pop)-1,-1,-1):
      if pop[x]==1:
        #print "In Vetter with input",last,"output",x
        return x
