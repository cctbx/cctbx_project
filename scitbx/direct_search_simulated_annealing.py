from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.stdlib import math, random
from libtbx.utils import Sorry
from copy import deepcopy
from scitbx import simplex
from six.moves import range

class dssa(object):
  """
  Directed Simplex Simulated Annealing
  http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/go_files/DSSA.pdf
  Hybrid simulated annealing and direct search method for nonlinear unconstrained global optimization
  """

  def __init__(self,
               dimension,
               matrix, # ndm * (ndm+1)
               evaluator,
               further_opt=False,
               n_candidate=None,
               tolerance=1e-8,
               max_iter=1e9,
               coolfactor = 0.6,
               T_ratio = 1e4,
               simplex_scale=0.2,
               monitor_cycle=11):

    self.max_iter = max_iter
    self.dimension=dimension
    self.tolerance=tolerance
    self.evaluator = evaluator
    self.coolfactor = coolfactor
    self.T_ratio = T_ratio
    self.simplex_scale = simplex_scale
    if(n_candidate is None):
      self.n_candidate=dimension+1
    else:
      self.n_candidate=n_candidate

    if((len(matrix) != self.dimension+1) or (matrix[0].size() != self.dimension)):
       raise Sorry("The initial simplex matrix does not match dimensions specified")
    for vector in matrix[1:]:
      if (vector.size() !=  matrix[0].size()):
        raise Sorry("Vector length in intial simplex do not match up" )

    self.monitor_cycle=monitor_cycle
    self.initialize(matrix)
    self.candidates = []
    self.optimize()
    if(further_opt):
      self.optimize_further()

  def initialize(self,matrix):
    self.end=False
    self.found=False
    self.simplexValue=flex.double()
    self.matrix=[]
    for vector in matrix:
     self.matrix.append(vector.deep_copy())
     self.simplexValue.append(self.function(vector))
    self.centroid=flex.double()
    self.reflectionPt=flex.double()

  def optimize(self):
    found = False
    end = False
    self.Nstep=self.dimension * 2
    monitor_score=0
    self.sort()

    for point in self.matrix:
      self.candidates.append(point.deep_copy() )
    self.candi_value = self.simplexValue.deep_copy()

    rd = 0
    self.count = 0
    while ((not found ) and (not end)):
      self.explore()
      self.update_candi()

      self.count += 1
      self.min_score=self.simplexValue[0]

      if self.count%self.monitor_cycle==0:
        rd = abs(self.simplexValue[self.dimension]-self.simplexValue[0])
        rd = rd/(abs(self.min_score)+self.tolerance*self.tolerance)
        if rd < self.tolerance:
          found = True
        else:
          monitor_score = self.min_score

      if self.count>=self.max_iter or self.T < self.min_T:
        end =True

  def update_candi(self):
    for ii in range(2):
      for jj in range(ii, self.n_candidate):
        if self.simplexValue[ii] < self.candi_value[jj]:
          if(self.simplexValue[ii] > self.candi_value[jj-1]):
            self.candi_value.insert( jj, self.simplexValue[ii])
            self.candidates.insert( jj, self.matrix[ii].deep_copy() )
            if(self.candi_value.size() > self.n_candidate):
              self.candi_value.pop_back()
              self.candidates.pop()
          break

  def optimize_further(self):
    self.solutions=[]
    self.scores= flex.double()
    for candi in self.candidates:
      starting_simplex = []
      for ii in range(self.dimension+1):
        starting_simplex.append(flex.random_double(self.dimension)*self.simplex_scale + candi)
      optimizer = simplex.simplex_opt(     dimension=self.dimension,
                                           matrix = starting_simplex,
                                           evaluator = self.evaluator,
                                           tolerance=self.tolerance
                                     )
      self.solutions.append( optimizer.get_solution() )
      self.scores.append( optimizer.get_score() )

    min_index = flex.min_index( self.scores )
    self.best_solution = self.solutions[ min_index ]
    self.best_score = self.scores[ min_index ]


  def explore(self):
     if(self.count == 0):
       self.T=(flex.mean(flex.pow2(self.simplexValue - flex.mean(self.simplexValue) )))**0.5 * 10.0
       self.min_T = self.T /self.T_ratio
     elif(self.count%self.Nstep == 0):
       self.T = self.T*self.coolfactor

     for kk in range(1,self.dimension+1):
       self.FindCentroidPt(self.dimension+1-kk)
       self.FindReflectionPt(kk)
     self.sort()

     return # end of this explore step

  def sort(self):
    tmp_matrix=deepcopy(self.matrix)
    tmp_value=list(self.simplexValue)
    sort_value=list(self.simplexValue)
    sort_value.sort()
    indx_array=[]
    for value in sort_value:
      indx=tmp_value.index(value)
      indx_array.append(indx)
    for ii in range(self.dimension+1):
      self.matrix[ii]=tmp_matrix[indx_array[ii]]
      self.simplexValue[ii]=tmp_value[indx_array[ii]]
    return

  def FindCentroidPt(self,kk):
   self.centroid=self.matrix[0]*0
   for ii in range (kk):
     self.centroid += self.matrix[ii]
   self.centroid /= kk

  def FindReflectionPt(self,kk):
    reflect_matrix=[]
    reflect_value=[]
    self.alpha=random.random()*0.2+0.9
    for ii in range(self.dimension+1-kk, self.dimension+1):
      reflectionPt=(self.centroid*(1.0+self.alpha) - self.alpha*self.matrix[ii])
      reflect_matrix.append(reflectionPt)
      reflect_value.append(self.function(reflectionPt))
    self.reflectionPtValue=min(reflect_value)
    if(self.reflectionPtValue > self.simplexValue[0]):
      p=math.exp(-(self.reflectionPtValue-self.simplexValue[0])/self.T)
      #print p
      if(p >= random.random()):
        self.ReplacePt(kk, reflect_matrix, reflect_value)
    else:
      self.ReplacePt(kk, reflect_matrix, reflect_value)



  def ReplacePt(self,kk,reflect_matrix,reflect_value):
    for ii in range(self.dimension+1-kk,self.dimension+1):
      self.matrix[ii] = reflect_matrix[ii-(self.dimension+1-kk)]
      self.simplexValue[ii]=reflect_value[ii-(self.dimension+1-kk)]
    self.sort()
    #self.update_candi()

  def get_solution(self):
    return self.best_solution

  def get_candi(self):
    return self.candidates

  def function(self,point):
    return self.evaluator.target( point )
