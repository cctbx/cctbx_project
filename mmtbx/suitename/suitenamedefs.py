"Class definitions used throughout suitename"

#        Copyright 2021  Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function
import numpy as np
from numpy import array
from enum import Enum


class Holder(object):
  "these objects exist only to have case specific attributes added to them"
  pass

globals = Holder()
# globals.options will be created in the main program


# reasons why a suite may fail to be classified:
Issue = Enum('Issue', 'DELTA_M EPSILON_M ZETA_M ALPHA BETA GAMMA DELTA')
reasons = {
  Issue.DELTA_M:    "delta-1",
  Issue.EPSILON_M:  "epsilon-1",
  Issue.ZETA_M:     "zeta-1",
  Issue.ALPHA:      "alpha",
  Issue.BETA:       "beta",
  Issue.GAMMA:      "gamma",
  Issue.DELTA:      "delta"
}

failMessages = {
  Issue.DELTA_M:    "bad deltam",
  Issue.GAMMA:      "g out",
  Issue.DELTA:      "bad delta"
}


class Bin(object):
  """ primary (coarse grained) classification of suites """
  # permanence properties
  name = ""
  ordinal = 0
  cluster = ()
  # a tuple of cluster objects
  dominant = -1
  # statistics gathered during the run
  count = 0
  active = False

  def __init__(self, ordinal, name, clusters=()):
    self.ordinal = ordinal
    self.name = name
    self.cluster = clusters
    self.dominant = -1
    self.active = False

    for i, c in enumerate(clusters):
      if c.dominance == "dom":
        self.dominant = i
        break


class Cluster(object):
  """secondary (fine grained) classification of suites"""
  # intrinsic data:
  ordinal = 0       # its place in the bin
  name = ""         # the name of the cluster
  status = ""       # certain, wannabe, triaged, outlier, nothing, incomplete
  clusterColor = "" # kinemage color names
  dominance = ""    # dom, sat, ord, out, tri, inc
  satelliteInfo = None    # present only if this cluster is a satellite
  angle = ()        # tuple of 9 angles:
  # chiMinus, deltaMinus, epsilon, zeta, alpha, beta, gamma, delta, chi

  # gathered statistics:
  count = 0  # number of data points found in this cluster
  suitenessSum = 0
  suitenessCounts = None

  def __init__(self, ordinal, name, status, color, dominance, angles):
    self.ordinal = ordinal
    self.name = name
    self.LOK = (name != "!!")
    self.status = status
    self.clusterColor = color
    self.dominance = dominance
    self.angle = array(angles)
    self.satelliteInfo = None
    # sometimes modified in suiteninit.buildBin
    self.suitenessCounts = np.zeros(12)
    self.suitenessSum = 0


class SatelliteInfo(object):
  # numbers used when suite is between satellite and dominant centers
  name = ""
  satelliteWidths = ()   # vector of 9 angles
  dominantWidths = ()    # vector of 9 angles

  def __init__(self, name, satelliteWidths, dominantWidths):
    self.name = name
    self.satelliteWidths = satelliteWidths
    self.dominantWidths = dominantWidths


class Residue(object):
  '''
  # A residue as normally read in, consisting of its six dihedral angles
  Used only briefly as input.
  '''
  sequence = -1  # sequence number in PDB file
  pointIDs = []
  base = " "    # A, C, G, U, ...
  angle = np.empty(0)  # will have 6 or 7 elements:
  # alpha, beta, gamma, delta, epsilon, zeta [, chi]

  def __init__(self, ID, base, angles):
    self.pointIDs = ID
    self.base = base
    self.angle = angles

  def is_dead(self):
    dead = [a == 9999 for a in self.angle]
    return all(dead)


  # def __str__(self):
  #     string = "%s:%s:%s:%4s:%s:%s:%s:%s:%s:%s:%s:%s:%s" \
  #        % (" ",
  #           "1",
  #           chainID,
  #           resnum,
  #           i_code,
  #           altloc,
  #           resname,
  #           alpha,
  #           beta,
  #           gamma,
  #           delta,
  #           epsilon,
  #           zeta)
  #     return string


  # nicknames: for ease of reading the code, each angle is given
  # a meaningful alias. Here they are:
  # 0   alpha
  # 1   beta
  # 2   gamma
  # 3   delta
  # 4   epsilon
  # 5   zeta
  # 7   chi
  @property
  def alpha(self):
    return self.angle[0]

  @alpha.setter
  def alpha(self, value):
    self.angle[0] = value

  @property
  def beta(self):
    return self.angle[1]

  @beta.setter
  def beta(self, value):
    self.angle[1] = value

  @property
  def gamma(self):
    return self.angle[2]

  @gamma.setter
  def gamma(self, value):
    self.angle[2] = value

  @property
  def delta(self):
    return self.angle[3]

  @delta.setter
  def delta(self, value):
    self.angle[3] = value

  @property
  def epsilon(self):
    return self.angle[4]

  @epsilon.setter
  def epsilon(self, value):
    self.angle[4] = value

  @property
  def zeta(self):
    return self.angle[5]

  @zeta.setter
  def zeta(self, value):
    self.angle[5] = value

  @property
  def chi(self):
    if len(self.angle) > 6:
      return self.angle[6]

  @chi.setter
  def chi(self, value):
    if len(self.angle) > 6:
      self.angle[6] = value


class Suite(object):
  '''
  The set of angles forming the linkage BETWEEN residues.
  This is the core data structure used in most operations of the program.
  '''
  pointID = ()
  base = " "    # A, C, G, U, ...
  angle = np.empty(0) # will become an np.array of 9 angles:
  # chiMinus, deltaMinus, epsilon, zeta, alpha, beta, gamma, delta, chi

  # fields computed during analysis:
  valid = False  # False means an incomplete, malformed suite
  cluster = None  # The cluster to which it is assigned
  suiteness = 0.0
  distance = 0.0
  situation = ""  # by what logical path this cluster was assigned
  pointMaster = ""
  pointColor = ""

  def __init__(self, ID, base, angles=None):
    self.pointID = ID
    self.base = base
    if angles is None:
      self.angle = np.full(9, 0.0)
    else:
      self.angle = angles
    self.suiteness = 0.0
    self.distance = 0.0
    self.notes = ""
    self.dbflag = False

  def validate(self):
    # make sure that angles deltaMinus through delta are reasonable
    self.valid = True
    for i in range(1, 8):
      if self.angle[i] < 0 or self.angle[i] > 360:
        self.valid = False
    return self.valid

  # nicknames: for ease of reading the code, each angle is given
  # a meaningful alias. Here they are:
  # 0   chiMinus
  # 1   deltaMinus
  # 2   epsilon
  # 3   zeta
  # 4   alpha
  # 5   beta
  # 6   gamma
  # 7   delta
  # 8   chi

  @property
  def chiMinus(self):
    return self.angle[0]

  @chiMinus.setter
  def chiMinus(self, value):
    self.angle[0] = value

  @property
  def deltaMinus(self):
    return self.angle[1]

  @deltaMinus.setter
  def deltaMinus(self, value):
    self.angle[1] = value

  @property
  def epsilon(self):
    return self.angle[2]

  @epsilon.setter
  def epsilon(self, value):
    self.angle[2] = value

  @property
  def zeta(self):
    return self.angle[3]

  @zeta.setter
  def zeta(self, value):
    self.angle[3] = value

  @property
  def alpha(self):
    return self.angle[4]

  @alpha.setter
  def alpha(self, value):
    self.angle[4] = value

  @property
  def beta(self):
    return self.angle[5]

  @beta.setter
  def beta(self, value):
    self.angle[5] = value

  @property
  def gamma(self):
    return self.angle[6]

  @gamma.setter
  def gamma(self, value):
    self.angle[6] = value

  @property
  def delta(self):
    return self.angle[7]

  @delta.setter
  def delta(self, value):
    self.angle[7] = value

  @property
  def chi(self):
    return self.angle[8]

  @chi.setter
  def chi(self, value):
    self.angle[8] = value


# The great variety of codes that may represent each base in the input file
NAListA = ":ADE:  A:A  : Ar:ATP:ADP:AMP:T6A:1MA:RIA:  I:I  :"
NAListG = ":GUA:  G:G  : Gr:GTP:GDP:GMP:GSP:1MG:2MG:M2G:OMG: 7MG:"
NAListC = ":CYT:  C:C  : Cr:CTP:CDP:CMP:5MC:OMC:"
NAListU = ":URA:URI:  U: Ur:U  :UTP:UDP:UMP:5MU:H2U:PSU:4SU:"
NAListY = ": YG:YG :  Y:Y  :"
#NAListT = ":THY:  T:T  : Tr:TTP:TDP:TMP:"
IgnoreDNAList = ": DA: DG: DC: DT:THY:  T:T  : Tr:TTP:TDP:TMP:"

# out of the noise, determine the base
def findBase(baseCode):
  if IgnoreDNAList.find(baseCode) >= 0:
    return None  # we ignore DNA residues
  elif len(baseCode) != 3:
    return "Z"

  if NAListA.find(baseCode) >= 0:
    base = "A"
  elif NAListG.find(baseCode) >= 0:
    base = "G"
  elif NAListC.find(baseCode) >= 0:
    base = "C"
  elif NAListU.find(baseCode) >= 0:
    base = "U"
  elif NAListY.find(baseCode) >= 0:
    base = "Y"
  else:
    base = "?"
  return base
