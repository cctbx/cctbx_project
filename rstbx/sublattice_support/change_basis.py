from __future__ import division, print_function
from scitbx import matrix
from six.moves import range

class sublattice_change_of_basis:

  def __init__(self,max_modulus=3):
    from cctbx.sgtbx.sub_lattice_tools import generate_matrix
    self.maxmod = max_modulus
    self._reindex_N = {}

    for matS in generate_matrix(1):
        self._reindex_N[1]=[sublattice_relation(matS,0)]

    for mod in range(2,self.maxmod+1):
      #if mod!=5:
        self._reindex_N[mod] = []
        for idx,matS in enumerate(generate_matrix(mod)):
          self._reindex_N[mod].append(sublattice_relation(matS,idx))

  def yield_transformations_ascending_modulus(self):
    for x in range(2,self.maxmod+1):
      for relation in self._reindex_N[x]:
        yield relation

  def yield_transformations_descending_modulus(self):
    for x in range(self.maxmod,1,-1):
      for relation in self._reindex_N[x]:
        yield relation

  def yield_identity_and_transformations_descending(self):
    yield self._reindex_N[1][0]
    for x in range(self.maxmod,1,-1):
      for relation in self._reindex_N[x]:
        yield relation

  def show_summary(self):
    for keymod in self._reindex_N:
      for relation in self._reindex_N[keymod]:
        relation.show_summary()

  def set_size(self):
    size=0
    for keymod in self._reindex_N: size+=len(self._reindex_N[keymod])
    return size

  def get_element(self,number):
    icount=0
    for item in self.yield_identity_and_transformations_descending():
      if icount==number:
        return item
      icount+=1

  def get_relation(self,code):
    found = [
        item
        for item in self.yield_identity_and_transformations_descending()
        if item.lookup_code()==code]
    if len(found)==0:  raise Exception("No index code=%s found"%(str(code)))
    return found[0]

class sublattice_relation:
  def __init__(self, transformation,idx):
    self._reindex_N = transformation # scitbx matrix sqr of rational ints
    self.printout_order = idx

  def matS(self):
    return self._reindex_N

  def lookup_code(self):
    return "%d.%d"%(self.index_N(),self.printout_order)

  def show_summary(self):
    self.show_M()
    # sublattice basis a',b',c'
    print (self.as_abc(mat=self._reindex_N,vec=["a","b","c"]))
    # index basis a,b,c
    print (self.as_abc(mat=self._reindex_N.inverse(),vec=["a'","b'","c'"]))
    print (self.sublattice_cosets())
    print ()

  def print_any_matrix(self, matrix):
    for i in range(3):
      for j in range(3):
        if j<2:
          print ("%4s"%(str(matrix[3*i+j])), end=" ")
        else:
          print ("%4s"%(str(matrix[3*i+j])))

  def show_M(self):
    # submitted manuscript, not consistent with Billiet or Rutherford
    # self.print_any_matrix(self._reindex_N.transpose())
    # revised manuscript
    self.print_any_matrix(self._reindex_N)

  def as_abc(self,mat,vec):
    basis_comp=[]
    for irow in [0,1,2]:
      component_sum=[]
      for idx,icol in enumerate([0,3,6]):
        fac = mat[irow+icol]
        if fac!=0:
          if fac==1:
            str_fac=''
          elif fac==-1:
            str_fac='-'
          else:
            str_fac="%s*"%(str(fac))
          component_sum.append("%s%s"%(str_fac,vec[idx]))
      all_plus_expression = "+".join(component_sum)
      basis_comp.append(all_plus_expression.replace("+-","-"))
    return ",".join(basis_comp)

  def index_N(self):
    return self._reindex_N.determinant()

  def sublattice_cosets(self):
    N = int(self.index_N()) # for Python 3 cast bost rational int to int
    cosets = {}
    inv_transform = self._reindex_N.inverse()

    ui_rows = [(i * inv_transform[0], i*inv_transform[1], i*inv_transform[2])
               for i in range(0,N+1)]
    uj_rows = [(j*inv_transform[3], j*inv_transform[4], j*inv_transform[5])
               for j in range(0,N+1)]
    uk_rows = [(k*inv_transform[6], k*inv_transform[7], k*inv_transform[8])
               for k in range(0,N+1)]
    for ui_row in ui_rows:
      for uj_row in uj_rows:
        for uk_row in uk_rows:
          alt_u_list = []
          for idx in range(3):
            #numerous Boost.python calls; could potentially be pushed to C++
            idx_suppl = ui_row[idx] + uj_row[idx] + uk_row[idx]
            den = idx_suppl.denominator()
            num = idx_suppl.numerator()
            alt_u_list.append( idx_suppl - (num//den))

          if alt_u_list != [0,0,0]:
            cosets[tuple(alt_u_list)]='present'

    assert len(list(cosets.keys()))==N-1
    return list(cosets.keys())
