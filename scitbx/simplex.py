from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.utils import Sorry
from six.moves import range

class simplex_opt(object):
  """
  Python implementaion of Nelder-Mead Simplex method
  The first publication is By JA Nelder and R. Mead:
    A simplex method for function minimization, computer journal 7(1965), 308-313
  The routine is following the description by Lagarias et al:
    Convergence properties of the nelder-mead simplex method in low dimensions SIAM J. Optim. Vol 9, No. 1, pp. 112-147
  And a previous implementation in (c++) by Adam Gurson and Virginia Torczon can be found at:
    http://www.cs.wm.edu/~va/software/SimplexSearch/AGEssay.html

  User need to specifiy the dimension and an initial 'guessed' simplex (n*(n+1)) matrix
    where n is the dimension

  A minimum is desired for simplex algorithm to work well, as the searching is not bounded.

  !! Target function should be modified to reflect specific purpose

  Problems should be referred to haiguang.liu@gmail.com
  """

  def __init__(self,
               dimension,
               matrix, # ndm * (ndm+1)
               evaluator,
               tolerance=1e-6,
               max_iter=1e9,
               alpha=1.0,
               beta=0.5,
               gamma=2.0,
               sigma=0.5,
               monitor_cycle=10):

    self.max_iter = max_iter
    self.dimension=dimension
    self.tolerance=tolerance
    self.evaluator = evaluator
    if((len(matrix) != self.dimension+1) or (matrix[0].size() != self.dimension)):
       raise Sorry("The initial simplex matrix does not match dimensions specified")
    for vector in matrix[1:]:
      if (vector.size() !=  matrix[0].size()):
        raise Sorry("Vector length in intial simplex do not match up" )

    self.alpha=alpha
    self.beta=beta
    self.gamma=gamma
    self.sigma=sigma
    self.monitor_cycle=monitor_cycle
    self.initialize(matrix)
    self.optimize()

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
    self.expansionPt=flex.double()
    self.contractionPt=flex.double()

    self.min_indx=0
    self.second_indx=0
    self.max_indx=0
    self.maxPrimePtId=0
    self.secondHigh=0.0

    self.max=0.0
    self.min=0.0

  def optimize(self):
    found = False
    end = False
    self.count = 0
    monitor_score=0
    while ((not found ) and (not end)):
      self.explore()
      self.count += 1
      self.min_score=self.simplexValue[self.min_indx]
      if self.count%self.monitor_cycle==0:
        rd = abs(monitor_score-self.min_score)
        #rd = abs(self.simplexValue[self.max_indx]-self.min_score)
        rd = rd/(abs(self.min_score)+self.tolerance*self.tolerance)
        if rd < self.tolerance:
          found = True
        else:
          monitor_score = self.min_score

      if self.count>=self.max_iter:
        end =True

  def explore(self):
     self.FindMinMaxIndices()
     self.FindCentroidPt()
     self.FindReflectionPt()
     self.secondHigh=self.simplexValue[self.second_indx]
     if self.simplexValue[self.min_indx] > self.reflectionPtValue:
        self.FindExpansionPt()
        if self.reflectionPtValue > self.expansionPtValue:
          self.ReplaceSimplexPoint(self.expansionPt)
          self.simplexValue[self.max_indx]=self.expansionPtValue
        else:
          self.ReplaceSimplexPoint(self.reflectionPt)
          self.simplexValue[self.max_indx]=self.reflectionPtValue
     elif (self.secondHigh > self.reflectionPtValue) and (self.reflectionPtValue >= self.simplexValue[self.min_indx]):
        self.ReplaceSimplexPoint(self.reflectionPt)
        self.simplexValue[self.max_indx]=self.reflectionPtValue
     elif self.reflectionPtValue >= self.secondHigh:
        self.FindContractionPt()
        if(self.maxPrimePtId == 0):
          if(self.contractionPtValue>self.maxPrimePtValue):
             self.ShrinkSimplex()
          else:
             self.ReplaceSimplexPoint(self.contractionPt)
             self.simplexValue[self.max_indx] = self.contractionPtValue
        elif (self.maxPrimePtId == 1):
          if(self.contractionPtValue >= self.maxPrimePtValue):
             self.ShrinkSimplex()
          else:
             self.ReplaceSimplexPoint(self.contractionPt)
             self.simplexValue[self.max_indx] = self.contractionPtValue
     return # end of this explore step

  def FindMinMaxIndices(self):
    self.max=self.min=second_max=self.simplexValue[0]
    self.max_indx=0
    self.min_indx=0
    self.second_indx=0
    for ii in range(1,self.dimension+1):
      if(self.simplexValue[ii] > self.max):
        second_max=self.max
        self.max=self.simplexValue[ii]
        self.second_indx=self.max_indx
        self.max_indx=ii
      elif(self.simplexValue[ii] > second_max):
        second_max=self.simplexValue[ii]
        self.second_indx=ii
      elif(self.simplexValue[ii] < self.min):
        self.min=self.simplexValue[ii]
        self.min_indx=ii
    return

  def FindCentroidPt(self):
   self.centroid=self.matrix[0]*0
   for ii in range (self.dimension+1):
      if(ii != self.max_indx):
        self.centroid += self.matrix[ii]
   self.centroid /= self.dimension

  def FindReflectionPt(self):
    self.reflectionPt=self.centroid*(1.0+self.alpha) -self.alpha*self.matrix[self.max_indx]
    self.reflectionPtValue=self.function(self.reflectionPt)

  def FindExpansionPt(self):
    self.expansionPt=self.centroid*(1.0-self.gamma)+self.gamma*self.reflectionPt
    self.expansionPtValue=self.function(self.expansionPt)

  def FindContractionPt(self):
    if(self.max <= self.reflectionPtValue):
      self.maxPrimePtId=1
      self.maxPrimePtValue=self.max
      self.contractionPt=self.centroid*(1.0-self.beta)+self.beta*self.matrix[self.max_indx]
    else:
      self.maxPrimePtId = 0
      self.maxPrimePtValue = self.reflectionPtValue
      self.contractionPt = self.centroid*(1.0-self.beta)+self.beta*self.reflectionPt
    self.contractionPtValue=self.function(self.contractionPt)

  def ShrinkSimplex(self):
    for ii in range(self.dimension+1):
      if(ii != self.min_indx):
        self.matrix[ii] = self.matrix[self.min_indx]+self.sigma*(self.matrix[ii] - self.matrix[self.min_indx])
        self.simplexValue[ii]=self.function(self.matrix[ii])

  def ReplaceSimplexPoint(self,vector):
    #self.FindCentroidPt()
    self.matrix[self.max_indx] = vector

  def get_solution(self):
    return self.matrix[self.min_indx]

  def get_score(self):
    return self.simplexValue[ self.min_indx ]

  def function(self,point):
    return self.evaluator.target( point )
