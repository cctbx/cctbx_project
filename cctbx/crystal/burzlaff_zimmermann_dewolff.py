'''
Reference implementation of BZW lattice characters, from Table 9.3.1 of
the International Tables.  Use of this method is not recommended.  See
counterexamples given below.

Literature:  Crystal lattices. H. Burzlaff, H. Zimmermann and P.M. de Wolff.
In International Tables for Cyrstallography, Volume A: Space-group symmetry
ed. Theo Hahn.  Dordrecht: Kluwer Academic Publishers (1996).
Note correction in the lattice_text table as given by W. Kabsch,
J. Appl. Cryst. 26: 795-800 (1993).  The 1996 edition of International
Tables is apparently uncorrected from the 1989 edition.

Author: N.K. Sauter
'''
from __future__ import absolute_import, division, print_function

from cctbx import crystal
from cctbx.uctbx import unit_cell
from scitbx import matrix
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
import math
import sys

class MetricCriteria(object):
  def equal(self,A,B):
    # A == B within stated fractional tolerance
    if self.int_tol < A/B < 1./self.int_tol or \
       self.int_tol < B/A < 1./self.int_tol: return True
    else: return False

  def equal_penalty(self,A,B):
    ratio = A/B
    # assume ratio is positive; it could only be negative for those lattice
    #  characters that have failed the tolerance tests, and penalty scores
    #  will be irrelevant in those cases.
    if ratio > 1.0: ratio = 1./ratio
    return 1.-ratio

  def zero(self,Product,B):
    """
    Product == 0 within fractional tolerance referred to by the quantities B.
    Usage: zero(D,(B,C)) | zero(E,(A,C)) | zero(F,(A,B))
    Only these combinations can be used.  This translates exactly into
    a tolerance on right angle, based on formula D=b.dot(c)=|b||c|cos(theta)
    internal_tolerance=0.97 corresponds to 1.72 degrees."""
    ldot = B[0]
    rdot = B[1]
    if abs(Product/math.sqrt(ldot*rdot)) < (1.-self.int_tol):  return True
    else: return False

  def zero_penalty(self,Product,B):
    ldot = B[0]
    rdot = B[1]
    return abs(Product/math.sqrt(ldot*rdot))

  def __init__(self,unit_cell,int_tol = 0.97):
    def criteria_encoding():
      return {

      1:[AB, AC, Type(1), equal(A/2.,D), equal(A/2.,E),equal(A/2.,F)],
      2:[AB, AC, Type(1), equal(E,D), equal(F,D)],
      3:[AB, AC, Type(2), D0, E0, F0],
      5:[AB, AC, Type(2), equal(-A/3.,D), equal(-A/3.,E),equal(-A/3.,F)],
      4:[AB, AC, Type(2), equal(E,D), equal(F,D)],
      6:[AB, AC, Type(2), Dstar, equal(E,D)],
      7:[AB, AC, Type(2), Dstar, equal(F,E)],
      8:[AB, AC, Type(2), Dstar],

       9:[AB, Type(1), equal(A/2.,D), equal(A/2.,E), equal(A/2.,F)],
      10:[AB, Type(1), equal(E,D)],
      11:[AB, Type(2), D0, E0, F0],
      12:[AB, Type(2), D0, E0, equal(-A/2.,F)],
      13:[AB, Type(2), D0, E0],
      15:[AB, Type(2), equal(-A/2.,D), equal(-A/2.,E), F0],
      16:[AB, Type(2), Dstar,equal(E,D)],
      14:[AB, Type(2), equal(E,D)],
      17:[AB, Type(2), Dstar],

      18:[BC, Type(1), equal(A/4.,D), equal(A/2.,E), equal(A/2.,F)],
      19:[BC, Type(1), equal(A/2.,E), equal(A/2.,F)],
      20:[BC, Type(1), equal(F,E)],
      21:[BC, Type(2), D0, E0, F0],
      22:[BC, Type(2), equal(-B/2.,D), E0, F0],
      23:[BC, Type(2), E0, F0],
      24:[BC, Type(2), Dstar, equal(-A/3.,E), equal(-A/3.,F)],
      25:[BC, Type(2), equal(F,E)],

      26:[Type(1), equal(A/4.,D), equal(A/2.,E), equal(A/2.,F)],
      27:[Type(1), equal(A/2.,E), equal(A/2.,F)],
      28:[Type(1), equal(A/2.,E), equal(2.*D,F)],
      29:[Type(1), equal(2.*D,E), equal(A/2.,F)],
      30:[Type(1), equal(B/2.,D), equal(2.*E,F)],
      31:[Type(1)],

      32:[Type(2), D0, E0, F0],
      40:[Type(2), equal(-B/2.,D), E0, F0],
      35:[Type(2), E0, F0],
      36:[Type(2), D0, equal(-A/2.,E), F0],
      33:[Type(2), D0, F0],
      38:[Type(2), D0, E0, equal(-A/2.,F)],
      34:[Type(2), D0, E0],
      42:[Type(2), equal(-B/2.,D), equal(-A/2.,E), F0],
      41:[Type(2), equal(-B/2.,D), F0],
      37:[Type(2), equal(-A/2.,E), F0],
      39:[Type(2), E0, equal(-A/2.,F)],
      43:[Type(2), Dstar, equal(abs(2*D+F),B)],
      44:[Type(2)],

    }
    self.int_tol = int_tol

    unit_cell = self.apply_sign_correction(unit_cell)

    mm = unit_cell.metrical_matrix()
    A=mm[0]; B=mm[1]; C=mm[2]; D=mm[5]; E=mm[4]; F=mm[3];
    #print "%.0f %.0f %.0f %.0f %.0f %.0f"%(A,B,C,D,E,F)

    #Part 1.  Encoded criteria give True / False answers
    equal = self.equal
    zero = self.zero
    D0 = zero(D,(B,C))
    E0 = zero(E,(A,C))
    F0 = zero(F,(A,B))
    if D*E*F <= 0: fType=2
    elif D0 or E0 or F0:
      # should never get here if sign correction has been applied
      # cell must be treated as Type 2 even though numerically Type 1
      fType=2
    else: fType=1
    def Type(flag):
      return flag==fType

    AB = equal(A,B)
    BC = equal(B,C)
    AC = equal(A,C)

    Dstar = equal(2*abs(D+E+F),A+B)

    self.metric_tests = criteria_encoding()

    #Part 2.  Encoded criteria give numerical penalty scores
    equal = self.equal_penalty
    zero = self.zero_penalty
    D0 = zero(D,(B,C))
    E0 = zero(E,(A,C))
    F0 = zero(F,(A,B))
    if D*E*F <= 0: type_penalty = (0.0,0.0,0.0) #type 2
    elif Type(2): type_penalty = (0.0,0.0,min(D0,E0,F0)) #type 2
    else: type_penalty = (0.0,0.0,0.0) # type 1
    def Type(flag):
      return type_penalty[flag]

    AB = equal(A,B)
    BC = equal(B,C)
    AC = equal(A,C)

    Dstar = equal(2*abs(D+E+F),A+B)

    self.penalties = criteria_encoding()

    if (False): # for debugging
      for key in self.penalties:
        print(key, self.metric_tests[key], [
          "%4.2f%%"%(100.*x) for x in self.penalties[key]])

  def apply_sign_correction(self,cell):
    # Function helps accomodate the tolerance.  Cell reduction may
    # produce a mathematically type=1 cell (DEF>0), which however has D~0, or
    # E~0 or F~0.  Therefore for practical purposes it is considered type=2
    # The basis must then be inverted to make it behave similarly to
    # a mathematically type=2 basis.
    # Several of the counterexamples produce correct results if this cor-
    # rection is applied; but wrong results if the correction is commented out

    mm = cell.metrical_matrix()
    A=mm[0]; B=mm[1]; C=mm[2]; D=mm[5]; E=mm[4]; F=mm[3];
    D0 = self.zero(D,(B,C))
    E0 = self.zero(E,(A,C))
    F0 = self.zero(F,(A,B))
    if D*E*F > 0 and (D0 or E0 or F0):
      #print "%.0f %.0f %.0f %.0f %.0f %.0f"%(A,B,C,D,E,F)
      # Program only gets here if either 0 or 2 of D,E,F are negative,
      # and if D~0 or E~0 or F~0
      orth = matrix.sqr(cell.orthogonalization_matrix())
      if D0 and E0:
        if F>0: inversion = matrix.sqr((-1.,0.,0.,0.,1.,0.,0.,0.,-1.))
      elif E0 and F0:
        if D>0: inversion = matrix.sqr((-1.,0.,0.,0.,-1.,0.,0.,0.,1.))
      elif F0 and D0:
        if E>0: inversion = matrix.sqr((1.,0.,0.,0.,-1.,0.,0.,0.,-1.))
      elif E0 and F>0: inversion = matrix.sqr((-1.,0.,0.,0.,1.,0.,0.,0.,-1.))
      elif F0 and D>0: inversion = matrix.sqr((-1.,0.,0.,0.,-1.,0.,0.,0.,1.))
      elif D0 and E>0: inversion = matrix.sqr((1.,0.,0.,0.,-1.,0.,0.,0.,-1.))
      uc = unit_cell(orthogonalization_matrix=orth*inversion)
      print("Adjusted cell:",uc)
      return uc
    return cell

  def character_tests(self,):
    return self.metric_tests

  def penalty_scores(self,):
    return self.penalties

lattice_text = '''
1         cubic  cF  1,-1,1,1,1,-1,-1,1,1
2  rhombohedral  hR  1,-1,0,-1,0,1,-1,-1,-1
3         cubic  cP  1,0,0,0,1,0,0,0,1
5         cubic  cI  1,0,1,1,1,0,0,1,1
4  rhombohedral  hR  1,-1,0,-1,0,1,-1,-1,-1
6    tetragonal  tI  0,1,1,1,0,1,1,1,0
7    tetragonal  tI  1,0,1,1,1,0,0,1,1
8  orthorhombic  oI  -1,-1,0,-1,0,-1,0,-1,-1
9  rhombohedral  hR  1,0,0,-1,1,0,-1,-1,3
10   monoclinic  mC  1,1,0,1,-1,0,0,0,-1
11   tetragonal  tP  1,0,0,0,1,0,0,0,1
12    hexagonal  hP  1,0,0,0,1,0,0,0,1
13 orthorhombic  oC  1,1,0,-1,1,0,0,0,1
15   tetragonal  tI  1,0,0,0,1,0,1,1,2
16 orthorhombic  oF  -1,-1,0,1,-1,0,1,1,2
14   monoclinic  mC  1,1,0,-1,1,0,0,0,1
17   monoclinic  mC  1,-1,0,-1,-1,0,-1,0,-1
18   tetragonal  tI  0,-1,1,1,-1,-1,1,0,0
19 orthorhombic  oI  -1,0,0,0,-1,1,-1,1,1
20   monoclinic  mC  0,1,1,0,1,-1,-1,0,0
21   tetragonal  tP  0,1,0,0,0,1,1,0,0
22    hexagonal  hP  0,1,0,0,0,1,1,0,0
23 orthorhombic  oC  0,1,1,0,-1,1,1,0,0
24 rhombohedral  hR  1,2,1,0,-1,1,1,0,0
25   monoclinic  mC  0,1,1,0,-1,1,1,0,0
26 orthorhombic  oF  1,0,0,-1,2,0,-1,0,2
27   monoclinic  mC  -1,2,0,-1,0,0,0,-1,1
28   monoclinic  mC  -1,0,0,-1,0,2,0,1,0
29   monoclinic  mC  1,0,0,1,-2,0,0,0,-1
30   monoclinic  mC  0,1,0,0,1,-2,-1,0,0
31    triclinic  aP  1,0,0,0,1,0,0,0,1
32 orthorhombic  oP  1,0,0,0,1,0,0,0,1
40 orthorhombic  oC  0,-1,0,0,1,2,-1,0,0
35   monoclinic  mP  0,-1,0,-1,0,0,0,0,-1
36 orthorhombic  oC  1,0,0,-1,0,-2,0,1,0
33   monoclinic  mP  1,0,0,0,1,0,0,0,1
38 orthorhombic  oC  -1,0,0,1,2,0,0,0,-1
34   monoclinic  mP  -1,0,0,0,0,-1,0,-1,0
42 orthorhombic  oI  -1,0,0,0,-1,0,1,1,2
41   monoclinic  mC  0,-1,-2,0,-1,0,-1,0,0
37   monoclinic  mC  1,0,2,1,0,0,0,1,0
39   monoclinic  mC  -1,-2,0,-1,0,0,0,0,-1
43   monoclinic  mI  -1,0,0,-1,-1,-2,0,-1,0
44    triclinic  aP  1,0,0,0,1,0,0,0,1
'''

lattices = []

for line in lattice_text.split('\n')[1:-1]:
  line = line.strip()
  items = line.split()
  d = {'number':int(items[0]),
       'system':items[1],
       'bravais':items[2],
       'matrix':tuple([int(x) for x in items[3].split(',')])
       }
  lattices.append(d)

class LatticeCharacter(object):
  #assumes that unit_cell is already in primitive setting and NIGGLI reduced
  # tolerance is given as a percentage
  def __init__(self,unit_cell,tolerance=3.):
    self.unit_cell = unit_cell
    self.int_tol = 1.0 - tolerance/100.
    self.possible_characters = self.best_bravais__()

  def best(self):
    return self.possible_charcters[0]

  def all(self):
    return self.possible_characters

  def show_summary(self):
    for item in self.possible_characters:
      print("%2d: %2s %14s, matrix"%(item['number'],item['bravais'],
            item['system']), end=' ')
      print("%30s"%str(item['matrix']), end=' ')
      print("max_penalty=%4.2f%%"%(100.*item['penalty']))

  def show_summary_labelit_format(self):
    for item in self.possible_characters:
      print("    %2d"%(item['number'],), end=' ')
      print("     %4.2f%%"%(100.*item['penalty'],), end=' ')
      print("               %14s %2s"%(item['system'],item['bravais'],), end=' ')
      print("%30s"%str(item['matrix']))

  def best_bravais__(self):
    MC = MetricCriteria(self.unit_cell,self.int_tol)
    answers = []
    for character in lattices:
      if not False in MC.character_tests()[character['number']]:
        character['penalty'] = max(MC.penalty_scores()[character['number']])
        answers.append(character)
    return answers

def counterexamples():

  for example in [
    # Stated example known to be hR, Table lookup gives mC
    (143.11386252141878, 143.60212158695211, 191.65636462983394,
     90.01331794194698, 111.85072371897915, 119.88503851099796),

    # Stated examples known to be oC, Table lookup gives mP
    # ...but applying sign correction produces oC.
    (121.48, 122.45, 144.06,  89.94,  65.09,  89.99),
    (121.39, 122.16, 143.98,  89.94,  65.24,  89.97),

    # Stated examples known to be hP, Table lookup gives oC
    # ...but applying sign correction produces hP
    (81.003144781582421, 81.130355147032134, 169.50959209067071,
     89.896542797449115, 89.999041351249176, 60.10397864241261),
    (81.88,  81.92, 170.38,  89.95,  89.98,  60.04),
    (81.84,  81.85, 169.28,  89.94,  89.97,  60.02),
    (149.74, 149.89, 154.90, 89.96,  89.77,  60.10),

    # Stated example known to be tP, Table lookup gives oP
    (64.5259, 65.5211, 140.646,  90.0599,  90.0254,  90.0023),

    # Stated example known to be tI, Table lookup gives oI
    (100.66, 100.78, 101.00, 116.96,  95.54, 116.63), ]:

    uc = unit_cell(example)
    assert uc.is_niggli_cell()
    print("Input Niggli cell:",uc)
    print("BZW Table lookup:")
    LatticeCharacter(uc,3.0).show_summary()

    input_symmetry = crystal.symmetry(
      unit_cell=uc,space_group_symbol="P 1")
    Groups = metric_subgroups(input_symmetry, 3.0)
    Groups.show()
    print();print()

def run():
  if len(sys.argv) < 2:
    counterexamples()
    sys.exit()
  # first six arguments are unit cell parameters a,b,c,alpha,beta,gamma
  # (degrees)
  params = [float(i) for i in sys.argv[1:7]]
  uc = unit_cell(params)
  niggli = uc.is_niggli_cell()
  print("Input cell:",uc,['is not Niggli reduced','is Niggli reduced'][niggli])
  if not niggli:
    uc = uc.niggli_cell()
    print("Niggli cell:",uc)
  # last argument is percent tolerance; default 3%
  tolerance = float(sys.argv[7])
  LatticeCharacter(uc,tolerance).show_summary()

if __name__=='__main__':
  run()
